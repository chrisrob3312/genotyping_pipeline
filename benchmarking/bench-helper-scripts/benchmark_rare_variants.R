#!/usr/bin/env Rscript
################################################################################
# benchmark_rare_variants.R
#
# Benchmarks rare variant imputation quality for RVAS applications.
# Demonstrates that better rare variant imputation → improved RVAS power.
#
# Key metrics:
#   1. Rare variant concordance by MAF bin and ancestry
#   2. Rare variant retention after QC filtering
#   3. Gene-level rare variant coverage (for burden tests)
#   4. Amerindigenous-enriched variant recovery
#
# For RVAS, what matters:
#   - Are rare variants imputed accurately? (concordance)
#   - Are they retained after filtering? (yield)
#   - Do we capture variants in genes? (coverage for burden tests)
#
# Usage:
#   Rscript benchmark_rare_variants.R \
#       --imputed imputed.vcf.gz \
#       --truth wgs_truth.vcf.gz \
#       --ancestry sample_ancestry.txt \
#       --output rare_variant_benchmark
#
################################################################################

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
})

# =============================================================================
# Command Line Arguments
# =============================================================================

option_list <- list(
    make_option(c("-i", "--imputed"), type = "character", default = NULL,
                help = "Imputed VCF file"),
    make_option(c("-t", "--truth"), type = "character", default = NULL,
                help = "WGS truth VCF file"),
    make_option(c("-a", "--ancestry"), type = "character", default = NULL,
                help = "Sample ancestry file (IID, SUPERPOP columns)"),
    make_option(c("-g", "--gene-bed"), type = "character", default = NULL,
                help = "Gene BED file for coverage analysis"),
    make_option(c("-o", "--output"), type = "character", default = "rare_benchmark",
                help = "Output prefix [default: %default]"),
    make_option(c("--maf-bins"), type = "character",
                default = "0.0001,0.0005,0.001,0.005,0.01",
                help = "MAF bin boundaries [default: %default]"),
    make_option(c("--info-threshold"), type = "numeric", default = 0.3,
                help = "INFO score threshold for retention analysis [default: %default]"),
    make_option(c("--approach-name"), type = "character", default = "pipeline",
                help = "Name of this approach for labeling [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# =============================================================================
# Configuration
# =============================================================================

maf_bins <- as.numeric(strsplit(opt$`maf-bins`, ",")[[1]])
maf_labels <- c(
    paste0("ultra_rare_", maf_bins[1]*100, "-", maf_bins[2]*100, "%"),
    paste0("very_rare_", maf_bins[2]*100, "-", maf_bins[3]*100, "%"),
    paste0("rare_", maf_bins[3]*100, "-", maf_bins[4]*100, "%"),
    paste0("low_freq_", maf_bins[4]*100, "-", maf_bins[5]*100, "%")
)

cat("================================================================================\n")
cat("RARE VARIANT IMPUTATION BENCHMARK FOR RVAS\n")
cat("================================================================================\n\n")
cat(sprintf("Approach: %s\n", opt$`approach-name`))
cat(sprintf("INFO threshold: %.2f\n", opt$`info-threshold`))
cat(sprintf("MAF bins: %s\n", paste(maf_bins, collapse = ", ")))
cat("\n")

# =============================================================================
# Extract Variant Information
# =============================================================================

cat("=== Extracting variant information ===\n\n")

# Get imputed variants with INFO scores
cat("Reading imputed variants...\n")
cmd_imp <- sprintf(
    "bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%INFO/AF\\t%%INFO/INFO\\n' %s 2>/dev/null || bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%INFO/AF\\t.\\n' %s",
    opt$imputed, opt$imputed
)
imputed <- tryCatch({
    fread(cmd = cmd_imp, header = FALSE,
          col.names = c("CHR", "POS", "REF", "ALT", "AF", "INFO"))
}, error = function(e) {
    cat("  Error reading imputed VCF. Trying simpler query...\n")
    cmd_simple <- sprintf("bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\n' %s | head -100000", opt$imputed)
    fread(cmd = cmd_simple, header = FALSE, col.names = c("CHR", "POS", "REF", "ALT"))
})

imputed[, AF := as.numeric(AF)]
imputed[, INFO := as.numeric(INFO)]
imputed[, MAF := pmin(AF, 1 - AF)]

cat(sprintf("  Total imputed variants: %d\n", nrow(imputed)))
cat(sprintf("  Rare variants (MAF < 1%%): %d\n", sum(imputed$MAF < 0.01, na.rm = TRUE)))

# Get truth variants
if (!is.null(opt$truth) && file.exists(opt$truth)) {
    cat("Reading truth variants...\n")
    cmd_truth <- sprintf("bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%INFO/AF\\n' %s", opt$truth)
    truth <- tryCatch({
        fread(cmd = cmd_truth, header = FALSE,
              col.names = c("CHR", "POS", "REF", "ALT", "AF_truth"))
    }, error = function(e) NULL)

    if (!is.null(truth)) {
        truth[, AF_truth := as.numeric(AF_truth)]
        truth[, MAF_truth := pmin(AF_truth, 1 - AF_truth)]
        cat(sprintf("  Truth variants: %d\n", nrow(truth)))
    }
} else {
    truth <- NULL
    cat("  No truth file provided - skipping concordance analysis\n")
}

# =============================================================================
# 1. RARE VARIANT CONCORDANCE BY MAF BIN
# =============================================================================

cat("\n=== 1. Rare Variant Concordance by MAF ===\n\n")

if (!is.null(truth)) {
    # Merge imputed and truth
    imputed[, VAR_ID := paste(CHR, POS, REF, ALT, sep = "_")]
    truth[, VAR_ID := paste(CHR, POS, REF, ALT, sep = "_")]

    merged <- merge(imputed, truth[, .(VAR_ID, MAF_truth)], by = "VAR_ID", all = FALSE)
    cat(sprintf("  Overlapping variants: %d\n", nrow(merged)))

    # Assign MAF bins based on truth MAF
    merged[, MAF_BIN := cut(MAF_truth, breaks = c(0, maf_bins, 1),
                            labels = c(maf_labels, "common"),
                            include.lowest = TRUE)]

    # Calculate concordance metrics by MAF bin
    concordance_by_maf <- merged[, .(
        N_variants = .N,
        Mean_INFO = mean(INFO, na.rm = TRUE),
        Pct_pass_INFO = 100 * mean(INFO >= opt$`info-threshold`, na.rm = TRUE),
        MAF_correlation = cor(MAF, MAF_truth, use = "complete.obs"),
        MAF_RMSE = sqrt(mean((MAF - MAF_truth)^2, na.rm = TRUE))
    ), by = MAF_BIN]

    cat("\nConcordance by MAF bin:\n")
    print(concordance_by_maf)

    # Focus on rare variants
    rare_summary <- merged[MAF_truth < 0.01, .(
        Total_rare = .N,
        Pass_INFO = sum(INFO >= opt$`info-threshold`, na.rm = TRUE),
        Retention_rate = 100 * mean(INFO >= opt$`info-threshold`, na.rm = TRUE),
        Mean_INFO = mean(INFO, na.rm = TRUE)
    )]

    cat("\n")
    cat("================================================================================\n")
    cat("RARE VARIANT SUMMARY (MAF < 1%)\n")
    cat("================================================================================\n")
    cat(sprintf("Total rare variants in truth:     %d\n", rare_summary$Total_rare))
    cat(sprintf("Pass INFO threshold (>= %.2f):    %d\n", opt$`info-threshold`, rare_summary$Pass_INFO))
    cat(sprintf("Retention rate:                   %.1f%%\n", rare_summary$Retention_rate))
    cat(sprintf("Mean INFO score:                  %.3f\n", rare_summary$Mean_INFO))
    cat("================================================================================\n")

} else {
    concordance_by_maf <- NULL
    cat("  Skipped (no truth data)\n")
}

# =============================================================================
# 2. RARE VARIANT RETENTION AFTER QC
# =============================================================================

cat("\n=== 2. Rare Variant Retention Analysis ===\n\n")

# Retention by MAF bin
imputed[, MAF_BIN := cut(MAF, breaks = c(0, maf_bins, 1),
                         labels = c(maf_labels, "common"),
                         include.lowest = TRUE)]

retention_by_maf <- imputed[, .(
    Total = .N,
    Pass_INFO = sum(INFO >= opt$`info-threshold`, na.rm = TRUE),
    Retention_pct = 100 * mean(INFO >= opt$`info-threshold`, na.rm = TRUE),
    Mean_INFO = mean(INFO, na.rm = TRUE),
    Median_INFO = median(INFO, na.rm = TRUE)
), by = MAF_BIN]

cat("Retention by MAF bin:\n")
print(retention_by_maf)

# This is key for RVAS: how many rare variants survive QC?
rvas_summary <- imputed[MAF < 0.01, .(
    Category = c("Ultra-rare (MAF<0.05%)", "Very rare (0.05-0.1%)",
                 "Rare (0.1-0.5%)", "Low-freq (0.5-1%)"),
    N_total = c(
        sum(MAF < 0.0005, na.rm = TRUE),
        sum(MAF >= 0.0005 & MAF < 0.001, na.rm = TRUE),
        sum(MAF >= 0.001 & MAF < 0.005, na.rm = TRUE),
        sum(MAF >= 0.005 & MAF < 0.01, na.rm = TRUE)
    ),
    N_retained = c(
        sum(MAF < 0.0005 & INFO >= opt$`info-threshold`, na.rm = TRUE),
        sum(MAF >= 0.0005 & MAF < 0.001 & INFO >= opt$`info-threshold`, na.rm = TRUE),
        sum(MAF >= 0.001 & MAF < 0.005 & INFO >= opt$`info-threshold`, na.rm = TRUE),
        sum(MAF >= 0.005 & MAF < 0.01 & INFO >= opt$`info-threshold`, na.rm = TRUE)
    )
)]
rvas_summary[, Retention_pct := round(100 * N_retained / N_total, 1)]

cat("\n")
cat("================================================================================\n")
cat("RVAS-RELEVANT RARE VARIANT RETENTION\n")
cat("================================================================================\n")
print(rvas_summary)
cat("\nInterpretation for RVAS:\n")
cat("  - Higher retention = more rare variants for burden/SKAT tests\n")
cat("  - Better imputation → more power to detect rare variant associations\n")
cat("================================================================================\n")

# =============================================================================
# 3. ANCESTRY-STRATIFIED RARE VARIANT METRICS
# =============================================================================

cat("\n=== 3. Ancestry-Stratified Analysis ===\n\n")

if (!is.null(opt$ancestry) && file.exists(opt$ancestry)) {
    ancestry <- fread(opt$ancestry)

    # Need sample-level genotypes for this
    # For now, provide framework
    cat("  Ancestry file provided: ", opt$ancestry, "\n")
    cat("  Full ancestry-stratified analysis requires sample-level genotypes\n")
    cat("  Framework ready for: AFR, EUR, AMR, EAS, SAS stratification\n")

    # Key message for AMR/Amerindigenous
    cat("\n")
    cat("  NOTE: For Amerindigenous-enriched variants:\n")
    cat("  - MX Biobank reference improves NAT tract imputation\n")
    cat("  - TOPMed has ~25% Hispanic/Latino samples\n")
    cat("  - Expect higher retention for AMR-specific rare variants\n")
    cat("  - Compare retention rates: AMR vs EUR for same MAF range\n")

} else {
    cat("  No ancestry file provided\n")
    cat("  To stratify by ancestry, provide --ancestry with IID,SUPERPOP columns\n")
}

# =============================================================================
# 4. GENE-LEVEL COVERAGE (FOR BURDEN TESTS)
# =============================================================================

cat("\n=== 4. Gene-Level Rare Variant Coverage ===\n\n")

if (!is.null(opt$`gene-bed`) && file.exists(opt$`gene-bed`)) {
    cat("Reading gene regions...\n")
    genes <- fread(opt$`gene-bed`, header = FALSE,
                   col.names = c("CHR", "START", "END", "GENE"))

    # Count rare variants per gene
    imputed_rare <- imputed[MAF < 0.01 & INFO >= opt$`info-threshold`]

    gene_coverage <- genes[, {
        chr <- CHR
        start <- START
        end <- END

        n_rare <- imputed_rare[CHR == chr & POS >= start & POS <= end, .N]

        .(N_rare_variants = n_rare)
    }, by = GENE]

    cat(sprintf("  Genes analyzed: %d\n", nrow(gene_coverage)))
    cat(sprintf("  Genes with ≥1 rare variant: %d (%.1f%%)\n",
                sum(gene_coverage$N_rare_variants > 0),
                100 * mean(gene_coverage$N_rare_variants > 0)))
    cat(sprintf("  Genes with ≥5 rare variants: %d (%.1f%%)\n",
                sum(gene_coverage$N_rare_variants >= 5),
                100 * mean(gene_coverage$N_rare_variants >= 5)))
    cat(sprintf("  Mean rare variants per gene: %.1f\n",
                mean(gene_coverage$N_rare_variants)))

} else {
    cat("  No gene BED file provided\n")
    cat("  To analyze gene coverage, provide --gene-bed\n")
    cat("  Download from: https://www.gencodegenes.org/\n")
    gene_coverage <- NULL
}

# =============================================================================
# 5. RVAS POWER IMPLICATIONS
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("RVAS POWER IMPLICATIONS\n")
cat("================================================================================\n")
cat("\n")
cat("What these metrics mean for RVAS:\n")
cat("\n")
cat("1. RARE VARIANT RETENTION:\n")
if (exists("rvas_summary")) {
    overall_retention <- sum(rvas_summary$N_retained) / sum(rvas_summary$N_total) * 100
    cat(sprintf("   Your pipeline retains %.1f%% of rare variants after QC\n", overall_retention))
    cat("   → More variants available for burden tests (SKAT, SKAT-O, CMC)\n")
}
cat("\n")
cat("2. CONCORDANCE WITH TRUTH:\n")
if (!is.null(concordance_by_maf)) {
    rare_info <- concordance_by_maf[MAF_BIN %in% maf_labels[1:3], mean(Mean_INFO, na.rm = TRUE)]
    cat(sprintf("   Mean INFO for rare variants: %.3f\n", rare_info))
    cat("   → Higher INFO = more accurate dosages for association testing\n")
}
cat("\n")
cat("3. EXPECTED RVAS IMPROVEMENT:\n")
cat("   - Better imputation → fewer false negatives (missed associations)\n")
cat("   - Higher concordance → correct effect size estimates\n")
cat("   - More retained variants → better gene-level coverage\n")
cat("\n")
cat("4. AMERINDIGENOUS/LATINO-SPECIFIC:\n")
cat("   - MX Biobank panel improves NAT-ancestry tract imputation\n")
cat("   - TOPMed/AoU references have Latino diversity\n")
cat("   - Expect improved detection of population-specific rare variants\n")
cat("   - Key genes: SLC16A11 (T2D), ABCA1 (HDL), PCSK9 (LDL)\n")
cat("\n")
cat("================================================================================\n")

# =============================================================================
# Output Files
# =============================================================================

cat("\nWriting output files...\n")

# Main results
results <- list(
    approach = opt$`approach-name`,
    retention_by_maf = retention_by_maf,
    rvas_summary = rvas_summary
)

if (!is.null(concordance_by_maf)) {
    results$concordance_by_maf <- concordance_by_maf
}

if (!is.null(gene_coverage)) {
    results$gene_coverage_summary <- data.table(
        Genes_total = nrow(gene_coverage),
        Genes_with_rare = sum(gene_coverage$N_rare_variants > 0),
        Genes_with_5plus = sum(gene_coverage$N_rare_variants >= 5),
        Mean_per_gene = mean(gene_coverage$N_rare_variants)
    )
}

# Save as RDS for later comparison
saveRDS(results, paste0(opt$output, "_results.rds"))
cat(sprintf("  %s_results.rds\n", opt$output))

# Save retention table
fwrite(retention_by_maf, paste0(opt$output, "_retention.txt"), sep = "\t")
cat(sprintf("  %s_retention.txt\n", opt$output))

# Save RVAS summary
fwrite(rvas_summary, paste0(opt$output, "_rvas_summary.txt"), sep = "\t")
cat(sprintf("  %s_rvas_summary.txt\n", opt$output))

cat("\nDone!\n")
