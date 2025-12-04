#!/usr/bin/env Rscript
################################################################################
# check_hit_recovery.R
#
# Checks whether known causal/GWAS variants are recovered after imputation.
# This is a key metric for comparing imputation approaches.
#
# Outputs:
#   - Recovery rate at different p-value thresholds
#   - Effect size correlation (observed vs expected)
#   - Power by MAF and ancestry
#
# Usage:
#   Rscript check_hit_recovery.R \
#       --gwas gwas_results.glm.linear \
#       --known-hits known_hits.txt \
#       --output hit_recovery_results
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
    make_option(c("-g", "--gwas"), type = "character", default = NULL,
                help = "GWAS results file (PLINK2 .glm.linear or .glm.logistic)"),
    make_option(c("-k", "--known-hits"), type = "character", default = NULL,
                help = "Known hits file (CHR, POS, RSID, BETA columns)"),
    make_option(c("-c", "--causal"), type = "character", default = NULL,
                help = "Causal variants from simulation (alternative to --known-hits)"),
    make_option(c("-o", "--output"), type = "character", default = "hit_recovery",
                help = "Output prefix [default: %default]"),
    make_option(c("--pvalue-thresholds"), type = "character",
                default = "5e-8,1e-6,1e-4,0.01,0.05",
                help = "P-value thresholds to test [default: %default]"),
    make_option(c("--window"), type = "integer", default = 500000,
                help = "Window around known hit to search (bp) [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# =============================================================================
# Read Data
# =============================================================================

cat("================================================================================\n")
cat("GWAS HIT RECOVERY ANALYSIS\n")
cat("================================================================================\n\n")

if (is.null(opt$gwas)) {
    stop("--gwas is required")
}

if (is.null(opt$`known-hits`) && is.null(opt$causal)) {
    stop("Either --known-hits or --causal is required")
}

# Read GWAS results
cat("Reading GWAS results...\n")
gwas <- fread(opt$gwas)

# Detect format (PLINK2 vs other)
if ("#CHROM" %in% names(gwas)) {
    # PLINK2 format
    setnames(gwas, "#CHROM", "CHR")
    gwas[, CHR := gsub("chr", "", CHR)]
} else if ("CHROM" %in% names(gwas)) {
    setnames(gwas, "CHROM", "CHR")
}

# Standardize column names
if ("P" %in% names(gwas)) {
    # Already named P
} else if ("P_VALUE" %in% names(gwas)) {
    setnames(gwas, "P_VALUE", "P")
}

cat(sprintf("  %d variants in GWAS results\n", nrow(gwas)))

# Read known hits
cat("Reading known hits...\n")
known_file <- ifelse(!is.null(opt$`known-hits`), opt$`known-hits`, opt$causal)
known <- fread(known_file)

# Standardize
if (!"CHR" %in% names(known) && "#CHR" %in% names(known)) {
    setnames(known, "#CHR", "CHR")
}
known[, CHR := as.character(CHR)]
known[, CHR := gsub("chr", "", CHR)]
gwas[, CHR := as.character(CHR)]

cat(sprintf("  %d known hits to check\n", nrow(known)))

# =============================================================================
# Check Recovery
# =============================================================================

cat("\nChecking hit recovery...\n")

pval_thresholds <- as.numeric(strsplit(opt$`pvalue-thresholds`, ",")[[1]])
window <- opt$window

# For each known hit, find best match in GWAS results
recovery_results <- lapply(1:nrow(known), function(i) {
    hit <- known[i]

    # Find variants in window
    matches <- gwas[CHR == hit$CHR &
                    POS >= (hit$POS - window) &
                    POS <= (hit$POS + window)]

    if (nrow(matches) == 0) {
        return(data.table(
            CHR = hit$CHR,
            POS = hit$POS,
            RSID = hit$RSID,
            EXPECTED_BETA = hit$BETA,
            FOUND = FALSE,
            BEST_P = NA,
            BEST_BETA = NA,
            DISTANCE = NA
        ))
    }

    # Find best match (lowest p-value)
    best_idx <- which.min(matches$P)
    best <- matches[best_idx]

    # Get observed beta
    obs_beta <- NA
    if ("BETA" %in% names(best)) {
        obs_beta <- best$BETA
    } else if ("OR" %in% names(best)) {
        obs_beta <- log(best$OR)
    }

    return(data.table(
        CHR = hit$CHR,
        POS = hit$POS,
        RSID = hit$RSID,
        EXPECTED_BETA = hit$BETA,
        FOUND = TRUE,
        BEST_P = best$P,
        BEST_BETA = obs_beta,
        DISTANCE = abs(best$POS - hit$POS),
        MATCH_POS = best$POS,
        MATCH_ID = if ("ID" %in% names(best)) best$ID else NA
    ))
})

results <- rbindlist(recovery_results, fill = TRUE)

# =============================================================================
# Calculate Recovery Rates
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("RECOVERY RATES BY P-VALUE THRESHOLD\n")
cat("================================================================================\n\n")

recovery_summary <- lapply(pval_thresholds, function(thresh) {
    n_found <- sum(results$FOUND, na.rm = TRUE)
    n_sig <- sum(results$BEST_P < thresh, na.rm = TRUE)
    n_total <- nrow(results)

    data.table(
        Threshold = thresh,
        N_Total = n_total,
        N_Found = n_found,
        N_Significant = n_sig,
        Pct_Found = 100 * n_found / n_total,
        Pct_Significant = 100 * n_sig / n_total
    )
})

recovery_summary <- rbindlist(recovery_summary)
print(recovery_summary)

# =============================================================================
# Effect Size Correlation
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("EFFECT SIZE CORRELATION\n")
cat("================================================================================\n\n")

# Only for variants with both expected and observed betas
complete <- results[FOUND == TRUE & !is.na(EXPECTED_BETA) & !is.na(BEST_BETA)]

if (nrow(complete) >= 3) {
    cor_test <- cor.test(complete$EXPECTED_BETA, complete$BEST_BETA)
    cat(sprintf("Correlation:     r = %.3f (p = %.2e)\n", cor_test$estimate, cor_test$p.value))
    cat(sprintf("R-squared:       %.3f\n", cor_test$estimate^2))

    # Regression slope (should be ~1 if unbiased)
    lm_fit <- lm(BEST_BETA ~ EXPECTED_BETA, data = complete)
    cat(sprintf("Slope:           %.3f (expected: 1.0)\n", coef(lm_fit)[2]))
    cat(sprintf("Intercept:       %.3f (expected: 0.0)\n", coef(lm_fit)[1]))
} else {
    cat("Too few variants with effect sizes for correlation analysis\n")
}

# =============================================================================
# Per-Hit Details
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("PER-HIT DETAILS\n")
cat("================================================================================\n\n")

# Sort by p-value
results_sorted <- results[order(BEST_P)]

# Print top hits
print(results_sorted[, .(RSID, CHR, POS, EXPECTED_BETA, BEST_P, BEST_BETA, DISTANCE)])

# =============================================================================
# Write Output
# =============================================================================

cat("\n")
cat("Writing output files...\n")

# Full results
fwrite(results, paste0(opt$output, "_per_variant.txt"), sep = "\t")
cat(sprintf("  %s_per_variant.txt\n", opt$output))

# Summary
fwrite(recovery_summary, paste0(opt$output, "_summary.txt"), sep = "\t")
cat(sprintf("  %s_summary.txt\n", opt$output))

# =============================================================================
# Quick Assessment
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("QUICK ASSESSMENT\n")
cat("================================================================================\n")

# Genome-wide significance recovery
gw_recovery <- 100 * sum(results$BEST_P < 5e-8, na.rm = TRUE) / nrow(results)
nominal_recovery <- 100 * sum(results$BEST_P < 0.05, na.rm = TRUE) / nrow(results)

if (gw_recovery >= 50) {
    cat(sprintf("✓ GOOD: %.0f%% of known hits reach genome-wide significance\n", gw_recovery))
} else if (gw_recovery >= 25) {
    cat(sprintf("○ MODERATE: %.0f%% of known hits reach genome-wide significance\n", gw_recovery))
} else {
    cat(sprintf("✗ LOW: Only %.0f%% of known hits reach genome-wide significance\n", gw_recovery))
}

cat(sprintf("  (%.0f%% reach nominal significance p < 0.05)\n", nominal_recovery))

if (exists("cor_test") && !is.na(cor_test$estimate)) {
    if (cor_test$estimate > 0.8) {
        cat(sprintf("✓ GOOD: Effect size correlation r = %.2f\n", cor_test$estimate))
    } else if (cor_test$estimate > 0.5) {
        cat(sprintf("○ MODERATE: Effect size correlation r = %.2f\n", cor_test$estimate))
    } else {
        cat(sprintf("✗ LOW: Effect size correlation r = %.2f\n", cor_test$estimate))
    }
}

cat("================================================================================\n")
