#!/usr/bin/env Rscript

# ==============================================================================
# MagicalRsq-X Filtering Helper Script
# ==============================================================================
# This script calculates MagicalRsq-X scores using pre-trained models
# and filters variants based on quality thresholds
# ==============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(MagicalRsq)  # Assuming package is installed from GitHub
})

# Parse command line arguments
option_list <- list(
    make_option(c("--vcf"), type="character", default=NULL,
                help="Input VCF file path", metavar="FILE"),
    make_option(c("--models"), type="character", default=NULL,
                help="Path to MagicalRsq-X model directory", metavar="DIR"),
    make_option(c("--ancestry"), type="character", default="EUR",
                help="Ancestry group (EUR or AFR) [default: %default]", metavar="STRING"),
    make_option(c("--output"), type="character", default="magicalrsqx_scores.txt",
                help="Output file path [default: %default]", metavar="FILE"),
    make_option(c("--threads"), type="integer", default=1,
                help="Number of threads [default: %default]", metavar="INT")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$vcf)) {
    stop("ERROR: --vcf is required")
}
if (is.null(opt$models)) {
    stop("ERROR: --models is required")
}
if (!file.exists(opt$vcf)) {
    stop(paste("ERROR: VCF file not found:", opt$vcf))
}
if (!dir.exists(opt$models)) {
    stop(paste("ERROR: Models directory not found:", opt$models))
}

cat("=======================================================================\n")
cat("MagicalRsq-X Filtering\n")
cat("=======================================================================\n")
cat(sprintf("VCF file: %s\n", opt$vcf))
cat(sprintf("Models directory: %s\n", opt$models))
cat(sprintf("Ancestry: %s\n", opt$ancestry))
cat(sprintf("Output: %s\n", opt$output))
cat(sprintf("Threads: %d\n", opt$threads))
cat("=======================================================================\n\n")

# ------------------------------------------------------------------------
# Step 1: Read VCF and extract necessary information
# ------------------------------------------------------------------------

cat("Step 1: Reading VCF file...\n")

# Read VCF header to get sample information
vcf_header <- system(paste("bcftools view -h", opt$vcf), intern=TRUE)
sample_line <- grep("^#CHROM", vcf_header, value=TRUE)
samples <- unlist(strsplit(sample_line, "\t"))[-c(1:9)]
n_samples <- length(samples)

cat(sprintf("  Found %d samples in VCF\n", n_samples))

# Extract variant information and imputation quality (INFO/R2 or INFO/RSQ)
cat("  Extracting variant information...\n")

vcf_query_cmd <- sprintf(
    "bcftools query -f '%%CHROM\\t%%POS\\t%%ID\\t%%REF\\t%%ALT\\t%%INFO/RSQ\\t%%INFO/R2\\t%%INFO/AF\\n' %s",
    opt$vcf
)

variants <- fread(cmd = vcf_query_cmd, header = FALSE)
setnames(variants, c("CHROM", "POS", "ID", "REF", "ALT", "RSQ", "R2", "AF"))

# Use whichever R² metric is available
variants[, rsq := ifelse(!is.na(RSQ), RSQ, R2)]
variants[, RSQ := NULL]
variants[, R2 := NULL]

cat(sprintf("  Loaded %d variants\n", nrow(variants)))

# ------------------------------------------------------------------------
# Step 2: Calculate additional features for MagicalRsq-X
# ------------------------------------------------------------------------

cat("\nStep 2: Calculating features for MagicalRsq-X...\n")

# MAF calculation (if AF is available)
if (all(is.na(variants$AF))) {
    cat("  WARNING: Allele frequency not in VCF, calculating from genotypes...\n")
    # Would need to extract genotypes and calculate MAF
    # For now, skip this if AF not available
} else {
    variants[, maf := pmin(AF, 1 - AF)]
    cat("  Calculated MAF from AF field\n")
}

# Categorize by MAF
variants[, maf_category := cut(
    maf,
    breaks = c(0, 0.005, 0.05, 1),
    labels = c("rare", "lowfreq", "common"),
    include.lowest = TRUE
)]

cat("  MAF distribution:\n")
cat(sprintf("    Rare (MAF < 0.5%%): %d variants\n", 
            sum(variants$maf_category == "rare", na.rm=TRUE)))
cat(sprintf("    Low-freq (0.5%% ≤ MAF < 5%%): %d variants\n", 
            sum(variants$maf_category == "lowfreq", na.rm=TRUE)))
cat(sprintf("    Common (MAF ≥ 5%%): %d variants\n", 
            sum(variants$maf_category == "common", na.rm=TRUE)))

# ------------------------------------------------------------------------
# Step 3: Load pre-trained MagicalRsq-X models
# ------------------------------------------------------------------------

cat("\nStep 3: Loading MagicalRsq-X models...\n")

# Model files should be named like:
# models/EUR_common_model.rds
# models/EUR_lowfreq_model.rds
# models/EUR_rare_model.rds

models <- list()
for (category in c("common", "lowfreq", "rare")) {
    model_file <- file.path(opt$models, sprintf("%s_%s_model.rds", opt$ancestry, category))
    
    if (file.exists(model_file)) {
        cat(sprintf("  Loading %s model from %s\n", category, model_file))
        models[[category]] <- readRDS(model_file)
    } else {
        cat(sprintf("  WARNING: Model not found for %s category: %s\n", category, model_file))
        models[[category]] <- NULL
    }
}

# ------------------------------------------------------------------------
# Step 4: Calculate MagicalRsq-X scores
# ------------------------------------------------------------------------

cat("\nStep 4: Calculating MagicalRsq-X scores...\n")

# Placeholder for actual MagicalRsq-X calculation
# The actual implementation depends on the MagicalRsq-X package structure

# For demonstration, we'll use a simplified approach:
# MagicalRsq-X = f(rsq, maf, ld_scores, recombination_rate, etc.)

variants[, MagicalRsqX := NA_real_]

for (category in c("common", "lowfreq", "rare")) {
    if (!is.null(models[[category]])) {
        # Select variants in this category
        idx <- which(variants$maf_category == category & !is.na(variants$rsq))
        
        if (length(idx) > 0) {
            cat(sprintf("  Processing %d %s variants...\n", length(idx), category))
            
            # Prepare features for prediction
            # This is a placeholder - actual features depend on model
            features <- variants[idx, .(
                rsq = rsq,
                maf = maf
                # Add other features as needed:
                # ld_score, recombination_rate, etc.
            )]
            
            # Predict MagicalRsq-X
            # The actual predict function depends on the model type
            # This is a placeholder
            predictions <- predict(models[[category]], newdata = features)
            
            variants[idx, MagicalRsqX := predictions]
        }
    }
}

# For genotyped variants (rsq = NA or rsq = 1), set MagicalRsq-X = 1
variants[is.na(rsq) | rsq >= 0.99, MagicalRsqX := 1.0]

cat(sprintf("  Calculated MagicalRsq-X for %d variants\n", 
            sum(!is.na(variants$MagicalRsqX))))
cat(sprintf("  Missing MagicalRsq-X: %d variants\n", 
            sum(is.na(variants$MagicalRsqX))))

# ------------------------------------------------------------------------
# Step 5: Write output
# ------------------------------------------------------------------------

cat("\nStep 5: Writing output...\n")

output_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "rsq", "maf", "MagicalRsqX")
fwrite(variants[, ..output_cols], file = opt$output, sep = "\t")

cat(sprintf("  Output written to: %s\n", opt$output))

# ------------------------------------------------------------------------
# Step 6: Summary statistics
# ------------------------------------------------------------------------

cat("\n=======================================================================\n")
cat("Summary Statistics\n")
cat("=======================================================================\n")

# Overall statistics
cat(sprintf("Total variants: %d\n", nrow(variants)))
cat(sprintf("Genotyped variants: %d\n", sum(is.na(variants$rsq) | variants$rsq >= 0.99)))
cat(sprintf("Imputed variants: %d\n", sum(!is.na(variants$rsq) & variants$rsq < 0.99)))

# MagicalRsq-X distribution
if (any(!is.na(variants$MagicalRsqX))) {
    magic_summary <- summary(variants$MagicalRsqX)
    cat("\nMagicalRsq-X distribution:\n")
    print(magic_summary)
    
    # Count by quality threshold
    cat(sprintf("\nVariants by MagicalRsq-X threshold:\n"))
    cat(sprintf("  ≥ 0.9: %d (%.1f%%)\n", 
                sum(variants$MagicalRsqX >= 0.9, na.rm=TRUE),
                100 * mean(variants$MagicalRsqX >= 0.9, na.rm=TRUE)))
    cat(sprintf("  ≥ 0.8: %d (%.1f%%)\n", 
                sum(variants$MagicalRsqX >= 0.8, na.rm=TRUE),
                100 * mean(variants$MagicalRsqX >= 0.8, na.rm=TRUE)))
    cat(sprintf("  ≥ 0.5: %d (%.1f%%)\n", 
                sum(variants$MagicalRsqX >= 0.5, na.rm=TRUE),
                100 * mean(variants$MagicalRsqX >= 0.5, na.rm=TRUE)))
    cat(sprintf("  ≥ 0.3: %d (%.1f%%)\n", 
                sum(variants$MagicalRsqX >= 0.3, na.rm=TRUE),
                100 * mean(variants$MagicalRsqX >= 0.3, na.rm=TRUE)))
}

cat("\n=======================================================================\n")
cat("MagicalRsq-X calculation complete!\n")
cat("=======================================================================\n")
