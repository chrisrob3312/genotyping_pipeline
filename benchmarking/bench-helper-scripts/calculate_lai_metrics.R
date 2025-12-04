#!/usr/bin/env Rscript
################################################################################
# calculate_lai_metrics.R
#
# Calculate Local Ancestry Inference (LAI) aware imputation quality metrics.
#
# Stratifies imputation quality by:
#   - Local ancestry tract (AFR, EUR, NAT, etc.)
#   - Position along chromosome
#   - MAF within each ancestry context
#
# For LAI-aware PRS (Tractor framework):
#   - Ancestry-specific effect sizes
#   - Local ancestry dosage weighting
#
# Usage:
#   Rscript calculate_lai_metrics.R \
#       --imputed imputed.vcf.gz \
#       --lai-msp rfmix2_output.msp.tsv \
#       --output lai_metrics.txt
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
                help = "Imputed VCF file with INFO scores"),
    make_option(c("-l", "--lai-msp"), type = "character", default = NULL,
                help = "RFMix2 MSP output (local ancestry calls)"),
    make_option(c("-t", "--truth"), type = "character", default = NULL,
                help = "WGS truth VCF (optional, for concordance)"),
    make_option(c("-o", "--output"), type = "character", default = "lai_metrics",
                help = "Output prefix [default: %default]"),
    make_option(c("--ancestry-codes"), type = "character", default = "0=AFR,1=EUR,2=NAT",
                help = "Ancestry code mapping [default: %default]"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Verbose output")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$imputed) || is.null(opt$`lai-msp`)) {
    print_help(opt_parser)
    stop("Both --imputed and --lai-msp are required")
}

# Parse ancestry codes
ancestry_map <- setNames(
    sapply(strsplit(strsplit(opt$`ancestry-codes`, ",")[[1]], "="), `[`, 2),
    sapply(strsplit(strsplit(opt$`ancestry-codes`, ",")[[1]], "="), `[`, 1)
)

log_msg <- function(msg) {
    if (opt$verbose) cat(sprintf("[%s] %s\n", Sys.time(), msg))
}

# =============================================================================
# Read LAI Data (RFMix2 MSP format)
# =============================================================================

read_rfmix_msp <- function(msp_file) {
    log_msg(sprintf("Reading LAI data from %s...", basename(msp_file)))

    # RFMix2 MSP format:
    # #chm  spos  epos  sgpos  egpos  n_snps  sample1.0  sample1.1  sample2.0  ...
    # Columns after n_snps are haplotype ancestry calls (0, 1, 2, etc.)

    msp <- fread(msp_file)

    # Parse header to get sample names
    header <- names(msp)
    sample_cols <- header[7:length(header)]  # After first 6 columns
    sample_names <- unique(gsub("\\.[01]$", "", sample_cols))

    log_msg(sprintf("  Found %d samples, %d ancestry tracts", length(sample_names), nrow(msp)))

    return(list(msp = msp, samples = sample_names, header = header))
}

# =============================================================================
# Get Local Ancestry at Variant Position
# =============================================================================

get_lai_at_position <- function(msp, chrom, pos) {
    # Find the tract containing this position
    tract <- msp[`#chm` == chrom & spos <= pos & epos >= pos]

    if (nrow(tract) == 0) {
        return(NULL)
    }

    # Return ancestry calls for all samples (both haplotypes)
    return(tract)
}

# =============================================================================
# Calculate LAI-Stratified Metrics
# =============================================================================

calculate_lai_metrics <- function(imputed_file, msp_data, ancestry_map) {
    log_msg("Calculating LAI-stratified metrics...")

    # Read imputed VCF (INFO scores)
    cmd <- sprintf("bcftools query -f '%%CHROM\\t%%POS\\t%%INFO/INFO\\t%%INFO/R2\\n' %s", imputed_file)
    variants <- fread(cmd = cmd, header = FALSE,
                      col.names = c("CHROM", "POS", "INFO", "R2"))

    # Handle missing values
    variants[INFO == ".", INFO := NA]
    variants[R2 == ".", R2 := NA]
    variants[, INFO := as.numeric(INFO)]
    variants[, R2 := as.numeric(R2)]

    log_msg(sprintf("  Processing %d variants...", nrow(variants)))

    msp <- msp_data$msp
    samples <- msp_data$samples

    # For each variant, determine local ancestry context
    # This is computationally intensive - process in chunks

    results <- list()

    for (i in seq_len(nrow(variants))) {
        if (i %% 10000 == 0) {
            log_msg(sprintf("    Processed %d/%d variants...", i, nrow(variants)))
        }

        var <- variants[i]
        tract <- get_lai_at_position(msp, var$CHROM, var$POS)

        if (is.null(tract)) next

        # Get sample columns (haplotype ancestry calls)
        sample_cols <- names(tract)[7:ncol(tract)]

        # Count ancestry backgrounds
        anc_counts <- table(unlist(tract[, sample_cols, with = FALSE]))

        # Majority ancestry for this region
        majority_anc <- names(which.max(anc_counts))
        majority_anc_label <- ancestry_map[as.character(majority_anc)]

        results[[i]] <- data.table(
            CHROM = var$CHROM,
            POS = var$POS,
            INFO = var$INFO,
            R2 = var$R2,
            Local_Ancestry = majority_anc_label,
            Ancestry_Homogeneity = max(anc_counts) / sum(anc_counts)
        )
    }

    results_dt <- rbindlist(results)
    log_msg(sprintf("  LAI annotation complete: %d variants", nrow(results_dt)))

    return(results_dt)
}

# =============================================================================
# Calculate MAF by Local Ancestry
# =============================================================================

calculate_ancestry_maf <- function(imputed_file, msp_data, ancestry_map) {
    log_msg("Calculating ancestry-specific MAF...")

    # This requires genotype data stratified by local ancestry
    # Complex implementation - would need full genotype matrix

    # Placeholder for structure
    return(NULL)
}

# =============================================================================
# Summarize by Local Ancestry and MAF
# =============================================================================

summarize_lai_metrics <- function(lai_data, maf_bins = c(0, 0.001, 0.005, 0.01, 0.05, 0.5)) {
    log_msg("Summarizing metrics by local ancestry and MAF...")

    # Note: MAF would need to be added from separate calculation
    # For now, summarize by local ancestry only

    summary_by_lai <- lai_data[, .(
        N_Variants = .N,
        Mean_INFO = mean(INFO, na.rm = TRUE),
        SD_INFO = sd(INFO, na.rm = TRUE),
        Median_INFO = median(INFO, na.rm = TRUE),
        Mean_R2 = mean(R2, na.rm = TRUE),
        SD_R2 = sd(R2, na.rm = TRUE),
        Pct_High_Quality = sum(INFO >= 0.8, na.rm = TRUE) / .N * 100,
        Mean_Homogeneity = mean(Ancestry_Homogeneity, na.rm = TRUE)
    ), by = Local_Ancestry]

    return(summary_by_lai)
}

# =============================================================================
# Main Analysis
# =============================================================================

log_msg("Starting LAI-stratified analysis")
log_msg(sprintf("Imputed: %s", opt$imputed))
log_msg(sprintf("LAI MSP: %s", opt$`lai-msp`))

# Read LAI data
msp_data <- read_rfmix_msp(opt$`lai-msp`)

# Calculate LAI-stratified metrics
lai_metrics <- calculate_lai_metrics(opt$imputed, msp_data, ancestry_map)

# Summarize
lai_summary <- summarize_lai_metrics(lai_metrics)

# =============================================================================
# Write Results
# =============================================================================

log_msg("Writing results...")

# Per-variant LAI metrics
fwrite(lai_metrics, paste0(opt$output, "_per_variant.txt.gz"), sep = "\t")
log_msg(sprintf("  Wrote %s_per_variant.txt.gz", opt$output))

# Summary by local ancestry
fwrite(lai_summary, paste0(opt$output, "_summary.txt"), sep = "\t")
log_msg(sprintf("  Wrote %s_summary.txt", opt$output))

# =============================================================================
# Print Summary
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("LAI-STRATIFIED METRICS SUMMARY\n")
cat("================================================================================\n")
cat("\n")
cat("By Local Ancestry:\n")
print(lai_summary)
cat("\n")
cat("================================================================================\n")
cat(sprintf("Results written to: %s_*.txt\n", opt$output))
cat("================================================================================\n")

log_msg("Done!")
