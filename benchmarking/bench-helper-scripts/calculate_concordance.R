#!/usr/bin/env Rscript
################################################################################
# calculate_concordance.R
#
# Calculate concordance metrics between imputed genotypes and WGS truth data.
#
# Metrics calculated:
#   - Genotype concordance (exact GT match)
#   - Dosage R² (correlation of imputed dosage vs truth)
#   - Non-reference concordance (NRC)
#   - Imputation accuracy by MAF bin
#
# Stratified by:
#   - MAF bins (0-0.5%, 0.5-1%, 1-5%, 5-50%)
#   - Ancestry group (if provided)
#   - Chromosome
#
# Usage:
#   Rscript calculate_concordance.R \
#       --imputed imputed.vcf.gz \
#       --truth wgs_truth.vcf.gz \
#       --output concordance_results.txt \
#       --ancestry sample_ancestry.txt \
#       --maf-bins "0,0.005,0.01,0.05,0.5"
#
################################################################################

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
})

# =============================================================================
# Parse Arguments
# =============================================================================

option_list <- list(
    make_option(c("-i", "--imputed"), type = "character", default = NULL,
                help = "Imputed VCF/PLINK file"),
    make_option(c("-t", "--truth"), type = "character", default = NULL,
                help = "WGS truth VCF/PLINK file"),
    make_option(c("-o", "--output"), type = "character", default = "concordance_results",
                help = "Output prefix [default: %default]"),
    make_option(c("-a", "--ancestry"), type = "character", default = NULL,
                help = "Sample ancestry file (tab-separated: IID, Ancestry)"),
    make_option(c("--maf-bins"), type = "character", default = "0,0.005,0.01,0.05,0.5",
                help = "MAF bin boundaries [default: %default]"),
    make_option(c("--samples"), type = "character", default = NULL,
                help = "File with samples to include (one per line)"),
    make_option(c("--chr"), type = "character", default = NULL,
                help = "Chromosome to process (for parallel processing)"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Verbose output")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$imputed) || is.null(opt$truth)) {
    print_help(opt_parser)
    stop("Both --imputed and --truth are required")
}

# Parse MAF bins
maf_bins <- as.numeric(strsplit(opt$`maf-bins`, ",")[[1]])

log_msg <- function(msg) {
    if (opt$verbose) {
        cat(sprintf("[%s] %s\n", Sys.time(), msg))
    }
}

# =============================================================================
# Helper Functions
# =============================================================================

#' Calculate genotype concordance
#' @param imputed Vector of imputed genotypes (0, 1, 2)
#' @param truth Vector of truth genotypes (0, 1, 2)
#' @return Concordance rate (0-1)
calc_concordance <- function(imputed, truth) {
    valid <- !is.na(imputed) & !is.na(truth)
    if (sum(valid) == 0) return(NA)
    sum(imputed[valid] == truth[valid]) / sum(valid)
}

#' Calculate dosage R²
#' @param imputed_dosage Vector of imputed dosages (0-2)
#' @param truth_dosage Vector of truth dosages (0-2)
#' @return R² (squared correlation)
calc_dosage_r2 <- function(imputed_dosage, truth_dosage) {
    valid <- !is.na(imputed_dosage) & !is.na(truth_dosage)
    if (sum(valid) < 10) return(NA)

    # Need variance in both
    if (var(imputed_dosage[valid]) == 0 || var(truth_dosage[valid]) == 0) return(NA)

    cor(imputed_dosage[valid], truth_dosage[valid])^2
}

#' Calculate non-reference concordance
#' @param imputed Vector of imputed genotypes
#' @param truth Vector of truth genotypes
#' @return NRC rate
calc_nrc <- function(imputed, truth) {
    # Only consider sites where truth has at least one alt allele
    has_alt <- truth > 0 & !is.na(truth)
    if (sum(has_alt) == 0) return(NA)

    # Concordance among non-reference sites
    sum(imputed[has_alt] == truth[has_alt], na.rm = TRUE) / sum(has_alt)
}

#' Calculate MAF from genotypes
#' @param geno Matrix of genotypes (samples x variants)
#' @return Vector of MAFs
calc_maf <- function(geno) {
    af <- colMeans(geno, na.rm = TRUE) / 2
    pmin(af, 1 - af)
}

#' Read VCF genotypes using bcftools
#' @param vcf_file Path to VCF file
#' @param samples Vector of sample IDs to include
#' @param region Genomic region (e.g., "chr22")
#' @return data.table with genotypes
read_vcf_genotypes <- function(vcf_file, samples = NULL, region = NULL) {
    # Build bcftools command
    cmd <- sprintf("bcftools query -f '%%CHROM\\t%%POS\\t%%ID\\t%%REF\\t%%ALT[\\t%%GT]\\n' %s", vcf_file)

    if (!is.null(region)) {
        cmd <- paste(cmd, "-r", region)
    }

    if (!is.null(samples)) {
        sample_file <- tempfile()
        writeLines(samples, sample_file)
        cmd <- paste(cmd, "-S", sample_file)
    }

    # Get sample names
    sample_cmd <- sprintf("bcftools query -l %s", vcf_file)
    if (!is.null(samples)) {
        sample_cmd <- paste(sample_cmd, "-S", sample_file)
    }
    sample_names <- system(sample_cmd, intern = TRUE)

    # Read genotypes
    log_msg(sprintf("Reading genotypes from %s...", basename(vcf_file)))
    geno_raw <- fread(cmd = cmd, header = FALSE, sep = "\t")

    if (nrow(geno_raw) == 0) {
        warning("No variants found in ", vcf_file)
        return(NULL)
    }

    # Set column names
    setnames(geno_raw, 1:5, c("CHROM", "POS", "ID", "REF", "ALT"))
    if (ncol(geno_raw) > 5) {
        setnames(geno_raw, 6:ncol(geno_raw), sample_names)
    }

    # Convert GT to dosage
    gt_to_dosage <- function(gt) {
        gt <- as.character(gt)
        ifelse(gt %in% c("0/0", "0|0"), 0,
               ifelse(gt %in% c("0/1", "1/0", "0|1", "1|0"), 1,
                      ifelse(gt %in% c("1/1", "1|1"), 2, NA)))
    }

    for (col in sample_names) {
        geno_raw[, (col) := gt_to_dosage(get(col))]
    }

    # Create variant ID
    geno_raw[, VAR_ID := paste(CHROM, POS, REF, ALT, sep = ":")]

    return(geno_raw)
}

# =============================================================================
# Main Analysis
# =============================================================================

log_msg("Starting concordance analysis")
log_msg(sprintf("Imputed: %s", opt$imputed))
log_msg(sprintf("Truth: %s", opt$truth))

# Read ancestry if provided
ancestry <- NULL
if (!is.null(opt$ancestry)) {
    log_msg("Reading ancestry assignments...")
    ancestry <- fread(opt$ancestry)
    setnames(ancestry, 1:2, c("IID", "Ancestry"))
}

# Read sample list if provided
samples <- NULL
if (!is.null(opt$samples)) {
    samples <- readLines(opt$samples)
    log_msg(sprintf("Restricting to %d samples", length(samples)))
}

# Read imputed genotypes
imputed <- read_vcf_genotypes(opt$imputed, samples = samples, region = opt$chr)
if (is.null(imputed)) {
    stop("Failed to read imputed genotypes")
}

# Read truth genotypes
truth <- read_vcf_genotypes(opt$truth, samples = samples, region = opt$chr)
if (is.null(truth)) {
    stop("Failed to read truth genotypes")
}

# Find overlapping variants and samples
log_msg("Finding overlapping variants...")
common_vars <- intersect(imputed$VAR_ID, truth$VAR_ID)
log_msg(sprintf("  Imputed variants: %d", nrow(imputed)))
log_msg(sprintf("  Truth variants: %d", nrow(truth)))
log_msg(sprintf("  Overlapping variants: %d", length(common_vars)))

if (length(common_vars) == 0) {
    stop("No overlapping variants found!")
}

# Get sample columns
sample_cols_imp <- setdiff(names(imputed), c("CHROM", "POS", "ID", "REF", "ALT", "VAR_ID"))
sample_cols_truth <- setdiff(names(truth), c("CHROM", "POS", "ID", "REF", "ALT", "VAR_ID"))
common_samples <- intersect(sample_cols_imp, sample_cols_truth)
log_msg(sprintf("Overlapping samples: %d", length(common_samples)))

if (length(common_samples) == 0) {
    stop("No overlapping samples found!")
}

# Subset to common variants
setkey(imputed, VAR_ID)
setkey(truth, VAR_ID)
imputed_sub <- imputed[common_vars, c("VAR_ID", "CHROM", "POS", common_samples), with = FALSE]
truth_sub <- truth[common_vars, c("VAR_ID", common_samples), with = FALSE]

# Ensure same order
setkey(imputed_sub, VAR_ID)
setkey(truth_sub, VAR_ID)

# Calculate MAF from truth data
log_msg("Calculating MAF from truth data...")
truth_geno_mat <- as.matrix(truth_sub[, common_samples, with = FALSE])
maf_truth <- calc_maf(truth_geno_mat)
imputed_sub[, MAF := maf_truth]

# Assign MAF bins
imputed_sub[, MAF_BIN := cut(MAF, breaks = maf_bins, include.lowest = TRUE,
                              labels = paste0(maf_bins[-length(maf_bins)] * 100, "-",
                                             maf_bins[-1] * 100, "%"))]

# =============================================================================
# Calculate Metrics
# =============================================================================

log_msg("Calculating concordance metrics...")

results <- list()

# Overall metrics
overall_concordance <- numeric(length(common_vars))
overall_nrc <- numeric(length(common_vars))
overall_dosage_r2 <- numeric(length(common_vars))

for (i in seq_len(length(common_vars))) {
    imp_geno <- as.numeric(imputed_sub[i, common_samples, with = FALSE])
    truth_geno <- as.numeric(truth_sub[i, common_samples, with = FALSE])

    overall_concordance[i] <- calc_concordance(imp_geno, truth_geno)
    overall_nrc[i] <- calc_nrc(imp_geno, truth_geno)
    overall_dosage_r2[i] <- calc_dosage_r2(imp_geno, truth_geno)
}

imputed_sub[, `:=`(
    Concordance = overall_concordance,
    NRC = overall_nrc,
    Dosage_R2 = overall_dosage_r2
)]

# Summary by MAF bin
log_msg("Summarizing by MAF bin...")
maf_summary <- imputed_sub[, .(
    N_variants = .N,
    Mean_Concordance = mean(Concordance, na.rm = TRUE),
    SD_Concordance = sd(Concordance, na.rm = TRUE),
    Mean_NRC = mean(NRC, na.rm = TRUE),
    SD_NRC = sd(NRC, na.rm = TRUE),
    Mean_Dosage_R2 = mean(Dosage_R2, na.rm = TRUE),
    SD_Dosage_R2 = sd(Dosage_R2, na.rm = TRUE)
), by = MAF_BIN]

results$maf_summary <- maf_summary

# Summary by ancestry (if provided)
if (!is.null(ancestry)) {
    log_msg("Calculating metrics by ancestry group...")

    ancestry_results <- list()

    for (anc in unique(ancestry$Ancestry)) {
        anc_samples <- ancestry[Ancestry == anc, IID]
        anc_samples <- intersect(anc_samples, common_samples)

        if (length(anc_samples) < 10) {
            log_msg(sprintf("  Skipping %s (only %d samples)", anc, length(anc_samples)))
            next
        }

        log_msg(sprintf("  Processing %s (%d samples)...", anc, length(anc_samples)))

        anc_concordance <- numeric(length(common_vars))
        anc_nrc <- numeric(length(common_vars))
        anc_dosage_r2 <- numeric(length(common_vars))

        for (i in seq_len(length(common_vars))) {
            imp_geno <- as.numeric(imputed_sub[i, anc_samples, with = FALSE])
            truth_geno <- as.numeric(truth_sub[i, anc_samples, with = FALSE])

            anc_concordance[i] <- calc_concordance(imp_geno, truth_geno)
            anc_nrc[i] <- calc_nrc(imp_geno, truth_geno)
            anc_dosage_r2[i] <- calc_dosage_r2(imp_geno, truth_geno)
        }

        # Summary by MAF for this ancestry
        anc_dt <- data.table(
            MAF_BIN = imputed_sub$MAF_BIN,
            Concordance = anc_concordance,
            NRC = anc_nrc,
            Dosage_R2 = anc_dosage_r2
        )

        anc_summary <- anc_dt[, .(
            N_variants = .N,
            N_samples = length(anc_samples),
            Mean_Concordance = mean(Concordance, na.rm = TRUE),
            SD_Concordance = sd(Concordance, na.rm = TRUE),
            Mean_NRC = mean(NRC, na.rm = TRUE),
            SD_NRC = sd(NRC, na.rm = TRUE),
            Mean_Dosage_R2 = mean(Dosage_R2, na.rm = TRUE),
            SD_Dosage_R2 = sd(Dosage_R2, na.rm = TRUE)
        ), by = MAF_BIN]

        anc_summary[, Ancestry := anc]
        ancestry_results[[anc]] <- anc_summary
    }

    results$ancestry_summary <- rbindlist(ancestry_results)
}

# Overall summary
overall_summary <- data.table(
    N_variants = length(common_vars),
    N_samples = length(common_samples),
    Mean_Concordance = mean(overall_concordance, na.rm = TRUE),
    SD_Concordance = sd(overall_concordance, na.rm = TRUE),
    Mean_NRC = mean(overall_nrc, na.rm = TRUE),
    SD_NRC = sd(overall_nrc, na.rm = TRUE),
    Mean_Dosage_R2 = mean(overall_dosage_r2, na.rm = TRUE),
    SD_Dosage_R2 = sd(overall_dosage_r2, na.rm = TRUE)
)
results$overall <- overall_summary

# =============================================================================
# Write Results
# =============================================================================

log_msg("Writing results...")

# Overall summary
fwrite(overall_summary, paste0(opt$output, "_overall.txt"), sep = "\t")
log_msg(sprintf("  Wrote %s_overall.txt", opt$output))

# MAF summary
fwrite(maf_summary, paste0(opt$output, "_by_maf.txt"), sep = "\t")
log_msg(sprintf("  Wrote %s_by_maf.txt", opt$output))

# Ancestry summary (if available)
if (!is.null(results$ancestry_summary)) {
    fwrite(results$ancestry_summary, paste0(opt$output, "_by_ancestry.txt"), sep = "\t")
    log_msg(sprintf("  Wrote %s_by_ancestry.txt", opt$output))
}

# Per-variant metrics
fwrite(imputed_sub[, .(VAR_ID, CHROM, POS, MAF, MAF_BIN, Concordance, NRC, Dosage_R2)],
       paste0(opt$output, "_per_variant.txt.gz"), sep = "\t")
log_msg(sprintf("  Wrote %s_per_variant.txt.gz", opt$output))

# =============================================================================
# Print Summary
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("CONCORDANCE SUMMARY\n")
cat("================================================================================\n")
cat(sprintf("Variants analyzed: %d\n", length(common_vars)))
cat(sprintf("Samples analyzed: %d\n", length(common_samples)))
cat("\n")
cat("Overall Metrics:\n")
cat(sprintf("  Genotype Concordance: %.4f (SD: %.4f)\n",
            overall_summary$Mean_Concordance, overall_summary$SD_Concordance))
cat(sprintf("  Non-Reference Concordance: %.4f (SD: %.4f)\n",
            overall_summary$Mean_NRC, overall_summary$SD_NRC))
cat(sprintf("  Dosage R²: %.4f (SD: %.4f)\n",
            overall_summary$Mean_Dosage_R2, overall_summary$SD_Dosage_R2))
cat("\n")
cat("By MAF Bin:\n")
print(maf_summary)
cat("\n")

if (!is.null(results$ancestry_summary)) {
    cat("By Ancestry:\n")
    print(results$ancestry_summary[, .(Ancestry, N_samples, Mean_Concordance, Mean_NRC, Mean_Dosage_R2)])
    cat("\n")
}

cat("================================================================================\n")
cat(sprintf("Results written to: %s_*.txt\n", opt$output))
cat("================================================================================\n")

log_msg("Done!")
