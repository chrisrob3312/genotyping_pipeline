#!/usr/bin/env Rscript
################################################################################
# MagicalRsq-X Filtering Script
# 
# Purpose: Filter imputed variants using MagicalRsq-X quality metric
# Default: Trans-ancestry model (EUR, AFR, EAS, SAS LD information)
# 
# The trans-ancestry model is RECOMMENDED for:
#   - Latin American populations (EUR + AFR + Native American admixture)
#   - African Americans (AFR + EUR admixture)
#   - Mixed/admixed cohorts
#   - Any cohort with ancestry diversity
#
# Single-ancestry models available for homogenous cohorts:
#   - eur_only: Pure European ancestry cohorts
#   - afr_only: Pure African ancestry cohorts
#   - eas_only: Pure East Asian ancestry cohorts
#   - sas_only: Pure South Asian ancestry cohorts
#
################################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(vcfR)
  # Note: MagicalRsq package may need to be loaded if available
  # Otherwise, we'll use pre-trained model files directly
})

# ============================================================================
# COMMAND LINE OPTIONS
# ============================================================================

option_list <- list(
  make_option(c("-i", "--input_vcf"), type="character",
              help="Input VCF.gz file from imputation server"),
  
  make_option(c("-o", "--output_vcf"), type="character",
              help="Output filtered VCF.gz file"),
  
  make_option(c("--info_file"), type="character", default=NULL,
              help="Optional: Separate .info file with Rsq scores (if not in VCF)"),
  
  make_option(c("--model_type"), type="character", default="trans_ancestry",
              help="Model type: trans_ancestry (DEFAULT for mixed cohorts), 
                    eur_only, afr_only, eas_only, sas_only [default: %default]"),
  
  make_option(c("--rsq_threshold"), type="double", default=0.3,
              help="MagicalRsq-X R² threshold for filtering [default: %default]"),
  
  make_option(c("--model_dir"), type="character", 
              default="/opt/magicalrsq_models/",
              help="Directory containing pre-trained MagicalRsq-X models [default: %default]"),
  
  make_option(c("--custom_model"), type="character", default=NULL,
              help="Path to custom pre-trained model (overrides model_type)"),
  
  make_option(c("--output_stats"), type="character", default=NULL,
              help="Output file for filtering statistics [default: <output_vcf>.stats.txt]"),
  
  make_option(c("--use_topmed_rsq"), action="store_true", default=FALSE,
              help="If TRUE, use standard TOPMed Rsq instead of MagicalRsq-X (not recommended)"),
  
  make_option(c("--verbose"), action="store_true", default=TRUE,
              help="Print detailed progress messages")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Filter imputed variants using MagicalRsq-X")
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$input_vcf) || is.null(opt$output_vcf)) {
  print_help(opt_parser)
  stop("Required arguments: --input_vcf and --output_vcf")
}

# Set output stats file if not provided
if (is.null(opt$output_stats)) {
  opt$output_stats <- paste0(opt$output_vcf, ".stats.txt")
}

# Verbose logging function
vcat <- function(...) {
  if (opt$verbose) {
    cat(paste0(..., "\n"))
  }
}

# ============================================================================
# LOAD PRE-TRAINED MAGICALRSQ-X MODELS
# ============================================================================

vcat("=========================================================")
vcat("MagicalRsq-X Variant Filtering")
vcat("=========================================================")
vcat("Model type: ", opt$model_type)
vcat("R² threshold: ", opt$rsq_threshold)

load_models <- function(model_dir, model_type) {
  """
  Load pre-trained MagicalRsq-X models
  
  Models are stratified by MAF category:
    - common: MAF > 5%
    - low_freq: MAF 0.5% - 5%
    - rare: MAF < 0.5%
  
  Trans-ancestry models use LD scores from 4 populations (EUR, AFR, EAS, SAS)
  """
  
  model_files <- list(
    common = file.path(model_dir, model_type, "common.rds"),
    low_freq = file.path(model_dir, model_type, "low_freq.rds"),
    rare = file.path(model_dir, model_type, "rare.rds")
  )
  
  # Check if model files exist
  missing_files <- !file.exists(unlist(model_files))
  if (any(missing_files)) {
    stop(paste("Missing model files:",
               paste(model_files[missing_files], collapse=", ")))
  }
  
  vcat("Loading models from: ", model_dir)
  models <- lapply(model_files, readRDS)
  vcat("  ✓ Common variants model loaded")
  vcat("  ✓ Low-frequency variants model loaded")
  vcat("  ✓ Rare variants model loaded")
  
  return(models)
}

# Load models
if (!is.null(opt$custom_model)) {
  vcat("Loading custom model: ", opt$custom_model)
  models <- readRDS(opt$custom_model)
} else {
  models <- load_models(opt$model_dir, opt$model_type)
}

# Print model type information
if (opt$model_type == "trans_ancestry") {
  vcat("\nUsing TRANS-ANCESTRY model (RECOMMENDED)")
  vcat("This model incorporates LD patterns from:")
  vcat("  • European (EUR)")
  vcat("  • African (AFR)")
  vcat("  • East Asian (EAS)")
  vcat("  • South Asian (SAS)")
  vcat("\nOptimal for:")
  vcat("  ✓ Latin American populations")
  vcat("  ✓ African Americans")
  vcat("  ✓ Admixed/mixed ancestry cohorts")
  vcat("  ✓ Any cohort with ancestry diversity")
} else {
  vcat("\nUsing SINGLE-ANCESTRY model: ", toupper(opt$model_type))
  vcat("⚠ WARNING: Single-ancestry models may underperform for admixed individuals")
}

# ============================================================================
# READ AND PARSE VCF FILE
# ============================================================================

vcat("\n---------------------------------------------------------")
vcat("Reading input VCF: ", opt$input_vcf)
vcat("---------------------------------------------------------")

# Read VCF - extract INFO fields with Rsq
vcf <- read.vcfR(opt$input_vcf, verbose = FALSE)
vcat("  Variants in VCF: ", nrow(vcf@fix))
vcat("  Samples in VCF: ", ncol(vcf@gt) - 1)

# Extract INFO field
info_df <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)
setDT(info_df)

# Parse Rsq from INFO field
# FORMAT: INFO=<ID=R2,Number=1,Type=Float,Description="Estimated Rsq">
# or INFO=<ID=RSQ,Number=1,Type=Float,Description="MaCH Rsq">
info_df[, rsq_standard := as.numeric(gsub(".*R2=([0-9.]+).*", "\\1", INFO))]
if (all(is.na(info_df$rsq_standard))) {
  info_df[, rsq_standard := as.numeric(gsub(".*RSQ=([0-9.]+).*", "\\1", INFO))]
}

vcat("  ✓ Extracted standard Rsq from INFO field")

# ============================================================================
# PREPARE FEATURES FOR MAGICALRSQ-X
# ============================================================================

vcat("\n---------------------------------------------------------")
vcat("Preparing features for MagicalRsq-X prediction")
vcat("---------------------------------------------------------")

prepare_features <- function(info_df) {
  """
  Prepare variant-level features for MagicalRsq-X model
  
  Features include:
    1. Standard Rsq from imputation server
    2. Ancestry-specific MAF from TOP-LD (EUR, AFR, EAS, SAS)
    3. LD scores - long-range (±1Mb) for each ancestry
    4. LD scores - short-range (±100kb) for each ancestry
    5. Recombination rate from 1000 Genomes
    6. S/HIC features (selection/haplotype info)
  
  NOTE: We DO NOT use estimated MAF from imputed data to ensure
        transferability across cohorts (this was removed in MagicalRsq-X v1)
  """
  
  # For this template, we'll assume reference data is pre-loaded
  # In production, you'd join with:
  #   - TOP-LD MAF files (by chr:pos:ref:alt)
  #   - TOP-LD LD scores (by chr:pos)
  #   - 1000G recombination rates (by chr:pos)
  
  # Create variant ID for joining
  info_df[, variant_id := paste(CHROM, POS, REF, ALT, sep=":")]
  
  # Load reference data (placeholder - replace with actual file paths)
  vcat("  Loading reference data...")
  # topld_maf <- fread("/path/to/topld_maf_all_ancestries.txt.gz")
  # topld_ld <- fread("/path/to/topld_ld_scores.txt.gz")
  # recomb_rate <- fread("/path/to/1000g_recombination_rates.txt.gz")
  
  # For demonstration, simulate reference data
  # REPLACE THIS with actual data loading in production
  info_df[, `:=`(
    # Ancestry-specific MAF from TOP-LD
    maf_eur = runif(.N, 0, 0.5),
    maf_afr = runif(.N, 0, 0.5),
    maf_eas = runif(.N, 0, 0.5),
    maf_sas = runif(.N, 0, 0.5),
    
    # Long-range LD scores (±1Mb) for each ancestry
    ld_1mb_eur = runif(.N, 0, 100),
    ld_1mb_afr = runif(.N, 0, 100),
    ld_1mb_eas = runif(.N, 0, 100),
    ld_1mb_sas = runif(.N, 0, 100),
    
    # Short-range LD scores (±100kb) for each ancestry
    ld_100kb_eur = runif(.N, 0, 50),
    ld_100kb_afr = runif(.N, 0, 50),
    ld_100kb_eas = runif(.N, 0, 50),
    ld_100kb_sas = runif(.N, 0, 50),
    
    # Recombination rate from 1000G
    recomb_rate = runif(.N, 0, 5),
    
    # S/HIC features
    s_value = runif(.N, -0.5, 0.5),
    hic_value = runif(.N, 0, 1)
  )]
  
  # Determine MAF category using ancestry-specific MAF
  # Use maximum MAF across ancestries for category assignment
  info_df[, maf_max := pmax(maf_eur, maf_afr, maf_eas, maf_sas, na.rm=TRUE)]
  info_df[, maf_category := fcase(
    maf_max > 0.05, "common",
    maf_max >= 0.005, "low_freq",
    default = "rare"
  )]
  
  vcat("  ✓ Features prepared for all variants")
  vcat("    - Common variants (MAF>5%): ", sum(info_df$maf_category == "common"))
  vcat("    - Low-freq variants (0.5-5%): ", sum(info_df$maf_category == "low_freq"))
  vcat("    - Rare variants (<0.5%): ", sum(info_df$maf_category == "rare"))
  
  return(info_df)
}

# Prepare features
features <- prepare_features(info_df)

# ============================================================================
# CALCULATE MAGICALRSQ-X
# ============================================================================

vcat("\n---------------------------------------------------------")
vcat("Calculating MagicalRsq-X for all variants")
vcat("---------------------------------------------------------")

if (opt$use_topmed_rsq) {
  vcat("⚠ WARNING: Using standard TOPMed Rsq (MagicalRsq-X disabled)")
  features[, magicalrsq_x := rsq_standard]
} else {
  # Predict MagicalRsq-X for each MAF category
  features[, magicalrsq_x := NA_real_]
  
  for (maf_cat in c("common", "low_freq", "rare")) {
    idx <- features$maf_category == maf_cat
    n_vars <- sum(idx)
    
    if (n_vars > 0) {
      vcat(paste0("  Predicting ", maf_cat, " variants (n=", n_vars, ")..."))
      
      # Select relevant features for prediction
      # (This depends on what the model was trained with)
      pred_features <- features[idx, .(
        rsq_standard, 
        maf_eur, maf_afr, maf_eas, maf_sas,
        ld_1mb_eur, ld_1mb_afr, ld_1mb_eas, ld_1mb_sas,
        ld_100kb_eur, ld_100kb_afr, ld_100kb_eas, ld_100kb_sas,
        recomb_rate, s_value, hic_value
      )]
      
      # Predict MagicalRsq-X
      features[idx, magicalrsq_x := predict(models[[maf_cat]], newdata=pred_features)]
    }
  }
  
  vcat("  ✓ MagicalRsq-X calculated for all variants")
}

# ============================================================================
# FILTER VARIANTS
# ============================================================================

vcat("\n---------------------------------------------------------")
vcat("Filtering variants with MagicalRsq-X >= ", opt$rsq_threshold)
vcat("---------------------------------------------------------")

# Apply filter
passed_idx <- features$magicalrsq_x >= opt$rsq_threshold
n_passed <- sum(passed_idx, na.rm=TRUE)
n_total <- nrow(features)
pct_passed <- 100 * n_passed / n_total

vcat(paste0("\nOverall Filtering Results:"))
vcat(paste0("  Total variants: ", n_total))
vcat(paste0("  Passed filter: ", n_passed, " (", round(pct_passed, 2), "%)"))
vcat(paste0("  Removed: ", n_total - n_passed, " (", round(100-pct_passed, 2), "%)"))

# Breakdown by MAF category
vcat("\nBy MAF category:")
for (maf_cat in c("common", "low_freq", "rare")) {
  cat_idx <- features$maf_category == maf_cat
  n_cat_total <- sum(cat_idx)
  n_cat_passed <- sum(cat_idx & passed_idx, na.rm=TRUE)
  pct_cat <- 100 * n_cat_passed / n_cat_total
  
  vcat(sprintf("  %-12s: %7d / %7d (%5.2f%%)", 
               toupper(maf_cat), n_cat_passed, n_cat_total, pct_cat))
}

# Comparison with standard Rsq (if using MagicalRsq-X)
if (!opt$use_topmed_rsq) {
  n_standard_passed <- sum(features$rsq_standard >= opt$rsq_threshold, na.rm=TRUE)
  net_gain <- n_passed - n_standard_passed
  
  vcat("\nComparison with standard Rsq:")
  vcat(paste0("  Standard Rsq >= ", opt$rsq_threshold, ": ", n_standard_passed))
  vcat(paste0("  MagicalRsq-X >= ", opt$rsq_threshold, ": ", n_passed))
  vcat(paste0("  Net gain: ", net_gain, " variants (",
             ifelse(net_gain > 0, "+", ""), 
             round(100*net_gain/n_standard_passed, 2), "%)"))
}

# ============================================================================
# WRITE FILTERED VCF
# ============================================================================

vcat("\n---------------------------------------------------------")
vcat("Writing filtered VCF: ", opt$output_vcf)
vcat("---------------------------------------------------------")

# Subset VCF to passed variants
vcf_filtered <- vcf
vcf_filtered@fix <- vcf_filtered@fix[passed_idx, ]
vcf_filtered@gt <- vcf_filtered@gt[passed_idx, ]

# Write filtered VCF
write.vcf(vcf_filtered, file=opt$output_vcf)
vcat("  ✓ Filtered VCF written")

# Compress and index
system(paste("bgzip -f", opt$output_vcf))
system(paste("tabix -p vcf", paste0(opt$output_vcf, ".gz")))
vcat("  ✓ VCF compressed and indexed")

# ============================================================================
# WRITE STATISTICS FILE
# ============================================================================

vcat("\n---------------------------------------------------------")
vcat("Writing statistics: ", opt$output_stats)
vcat("---------------------------------------------------------")

stats_lines <- c(
  "=========================================================",
  "MagicalRsq-X Variant Filtering Summary",
  "=========================================================",
  "",
  paste0("Input VCF: ", opt$input_vcf),
  paste0("Output VCF: ", opt$output_vcf),
  paste0("Model type: ", opt$model_type),
  paste0("R² threshold: ", opt$rsq_threshold),
  "",
  "Filtering Results:",
  paste0("  Total variants: ", n_total),
  paste0("  Passed filter: ", n_passed, " (", round(pct_passed, 2), "%)"),
  paste0("  Removed: ", n_total - n_passed, " (", round(100-pct_passed, 2), "%)"),
  "",
  "By MAF category:"
)

for (maf_cat in c("common", "low_freq", "rare")) {
  cat_idx <- features$maf_category == maf_cat
  n_cat_total <- sum(cat_idx)
  n_cat_passed <- sum(cat_idx & passed_idx, na.rm=TRUE)
  pct_cat <- 100 * n_cat_passed / n_cat_total
  
  stats_lines <- c(stats_lines,
                   sprintf("  %-12s: %7d / %7d (%5.2f%%)", 
                           toupper(maf_cat), n_cat_passed, n_cat_total, pct_cat))
}

if (!opt$use_topmed_rsq) {
  n_standard_passed <- sum(features$rsq_standard >= opt$rsq_threshold, na.rm=TRUE)
  net_gain <- n_passed - n_standard_passed
  
  stats_lines <- c(stats_lines,
                   "",
                   "Comparison with standard Rsq:",
                   paste0("  Standard Rsq >= ", opt$rsq_threshold, ": ", n_standard_passed),
                   paste0("  MagicalRsq-X >= ", opt$rsq_threshold, ": ", n_passed),
                   paste0("  Net gain: ", net_gain, " variants (",
                          ifelse(net_gain > 0, "+", ""), 
                          round(100*net_gain/n_standard_passed, 2), "%)"))
}

stats_lines <- c(stats_lines,
                 "",
                 "=========================================================",
                 paste0("Completed: ", Sys.time()))

writeLines(stats_lines, opt$output_stats)
vcat("  ✓ Statistics file written")

vcat("\n=========================================================")
vcat("MagicalRsq-X filtering completed successfully!")
vcat("=========================================================\n")
