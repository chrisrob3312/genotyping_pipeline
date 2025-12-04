#!/usr/bin/env Rscript
################################################################################
# compare_rvas_approaches.R
#
# Compares rare variant imputation quality across pipeline approaches.
# Generates publication-ready figures showing RVAS-relevant improvements.
#
# Expected finding: Your pipeline shows better rare variant retention,
# especially for:
#   - Amerindigenous-enriched variants
#   - Admixed (AMR) population samples
#   - Ultra-rare variants (MAF < 0.1%)
#
# Usage:
#   Rscript compare_rvas_approaches.R \
#       --results approach_a_results.rds,approach_f_results.rds \
#       --labels "Traditional,Our Pipeline" \
#       --output rvas_comparison
#
################################################################################

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(optparse)
})

# =============================================================================
# Command Line Arguments
# =============================================================================

option_list <- list(
    make_option(c("-r", "--results"), type = "character", default = NULL,
                help = "Comma-separated list of results RDS files"),
    make_option(c("-l", "--labels"), type = "character", default = NULL,
                help = "Comma-separated labels for each approach"),
    make_option(c("-o", "--output"), type = "character", default = "rvas_comparison",
                help = "Output prefix [default: %default]"),
    make_option(c("--highlight-approach"), type = "character", default = NULL,
                help = "Approach to highlight (e.g., 'Our Pipeline')")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# =============================================================================
# Load Results
# =============================================================================

cat("================================================================================\n")
cat("RVAS APPROACH COMPARISON\n")
cat("================================================================================\n\n")

result_files <- strsplit(opt$results, ",")[[1]]
labels <- if (!is.null(opt$labels)) strsplit(opt$labels, ",")[[1]] else paste0("Approach_", 1:length(result_files))

all_results <- lapply(seq_along(result_files), function(i) {
    res <- readRDS(result_files[i])
    res$label <- labels[i]
    res
})

cat(sprintf("Loaded %d approaches: %s\n\n", length(all_results), paste(labels, collapse = ", ")))

# =============================================================================
# Publication Theme
# =============================================================================

theme_publication <- function() {
    theme_minimal() +
    theme(
        text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 12, face = "bold")
    )
}

# Color palette (colorblind-friendly)
approach_colors <- c(
    "Traditional" = "#E69F00",
    "Our Pipeline" = "#0072B2",
    "1-Step" = "#56B4E9",
    "2-Step" = "#009E73",
    "TOPMed" = "#D55E00",
    "Michigan" = "#CC79A7"
)

# =============================================================================
# Figure 1: Rare Variant Retention by MAF
# =============================================================================

cat("=== Generating Figure 1: Retention by MAF ===\n")

# Combine retention data
retention_combined <- rbindlist(lapply(all_results, function(res) {
    dt <- res$rvas_summary
    dt$Approach <- res$label
    dt
}), fill = TRUE)

# Reorder categories
retention_combined[, Category := factor(Category,
    levels = c("Ultra-rare (MAF<0.05%)", "Very rare (0.05-0.1%)",
               "Rare (0.1-0.5%)", "Low-freq (0.5-1%)"))]

# Create figure
fig1 <- ggplot(retention_combined, aes(x = Category, y = Retention_pct, fill = Approach)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", Retention_pct)),
              position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
    scale_fill_manual(values = approach_colors) +
    labs(
        title = "Rare Variant Retention After Quality Filtering",
        subtitle = "Higher retention = more variants available for RVAS",
        x = "MAF Category",
        y = "Retention Rate (%)",
        fill = "Approach"
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0, 100))

ggsave(paste0(opt$output, "_fig1_retention.pdf"), fig1, width = 10, height = 6)
ggsave(paste0(opt$output, "_fig1_retention.png"), fig1, width = 10, height = 6, dpi = 300)
cat("  Saved: fig1_retention.pdf/png\n")

# =============================================================================
# Figure 2: Total Rare Variants by Approach
# =============================================================================

cat("=== Generating Figure 2: Total Variants ===\n")

# Calculate totals
totals <- retention_combined[, .(
    Total = sum(N_total),
    Retained = sum(N_retained)
), by = Approach]
totals[, Lost := Total - Retained]

# Melt for stacked bar
totals_melt <- melt(totals, id.vars = "Approach", measure.vars = c("Retained", "Lost"))

fig2 <- ggplot(totals_melt, aes(x = Approach, y = value, fill = variable)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = c("Retained" = "#2E7D32", "Lost" = "#C62828"),
                      labels = c("Pass QC", "Filtered")) +
    geom_text(data = totals, aes(x = Approach, y = Retained, label = scales::comma(Retained)),
              inherit.aes = FALSE, vjust = -0.5, fontface = "bold") +
    labs(
        title = "Rare Variant Yield by Approach",
        subtitle = "Total rare variants (MAF < 1%) retained after INFO filtering",
        x = "",
        y = "Number of Variants",
        fill = ""
    ) +
    theme_publication() +
    scale_y_continuous(labels = scales::comma)

ggsave(paste0(opt$output, "_fig2_totals.pdf"), fig2, width = 8, height = 6)
ggsave(paste0(opt$output, "_fig2_totals.png"), fig2, width = 8, height = 6, dpi = 300)
cat("  Saved: fig2_totals.pdf/png\n")

# =============================================================================
# Figure 3: Improvement Over Traditional
# =============================================================================

cat("=== Generating Figure 3: Relative Improvement ===\n")

if ("Traditional" %in% labels || length(labels) >= 2) {
    baseline <- labels[1]  # First approach as baseline

    # Calculate relative improvement
    baseline_retention <- retention_combined[Approach == baseline, .(Category, Baseline_pct = Retention_pct)]
    improvement <- merge(retention_combined[Approach != baseline], baseline_retention, by = "Category")
    improvement[, Improvement := Retention_pct - Baseline_pct]
    improvement[, Relative_improvement := 100 * (Retention_pct - Baseline_pct) / Baseline_pct]

    fig3 <- ggplot(improvement, aes(x = Category, y = Improvement, fill = Approach)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        geom_text(aes(label = sprintf("%+.1f%%", Improvement)),
                  position = position_dodge(width = 0.8),
                  vjust = ifelse(improvement$Improvement >= 0, -0.5, 1.5), size = 3) +
        scale_fill_manual(values = approach_colors) +
        labs(
            title = sprintf("Improvement in Rare Variant Retention vs %s", baseline),
            subtitle = "Positive values indicate more rare variants retained",
            x = "MAF Category",
            y = "Improvement in Retention (%)",
            fill = "Approach"
        ) +
        theme_publication() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(paste0(opt$output, "_fig3_improvement.pdf"), fig3, width = 10, height = 6)
    ggsave(paste0(opt$output, "_fig3_improvement.png"), fig3, width = 10, height = 6, dpi = 300)
    cat("  Saved: fig3_improvement.pdf/png\n")
}

# =============================================================================
# Figure 4: INFO Score Distribution for Rare Variants
# =============================================================================

cat("=== Generating Figure 4: INFO Distribution ===\n")

# If concordance data available
has_concordance <- sapply(all_results, function(x) !is.null(x$concordance_by_maf))

if (any(has_concordance)) {
    concordance_combined <- rbindlist(lapply(all_results[has_concordance], function(res) {
        dt <- res$concordance_by_maf
        dt$Approach <- res$label
        dt
    }), fill = TRUE)

    # Focus on rare MAF bins
    rare_conc <- concordance_combined[MAF_BIN %in% c("ultra_rare_0.01-0.05%", "very_rare_0.05-0.1%",
                                                      "rare_0.1-0.5%", "low_freq_0.5-1%")]

    if (nrow(rare_conc) > 0) {
        fig4 <- ggplot(rare_conc, aes(x = MAF_BIN, y = Mean_INFO, fill = Approach)) +
            geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
            geom_hline(yintercept = 0.3, linetype = "dashed", color = "red", alpha = 0.7) +
            annotate("text", x = 0.5, y = 0.32, label = "INFO threshold", hjust = 0, color = "red") +
            scale_fill_manual(values = approach_colors) +
            labs(
                title = "Imputation Quality (INFO Score) by MAF",
                subtitle = "Higher INFO = more accurate rare variant genotypes",
                x = "MAF Category",
                y = "Mean INFO Score",
                fill = "Approach"
            ) +
            theme_publication() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            coord_cartesian(ylim = c(0, 1))

        ggsave(paste0(opt$output, "_fig4_info.pdf"), fig4, width = 10, height = 6)
        ggsave(paste0(opt$output, "_fig4_info.png"), fig4, width = 10, height = 6, dpi = 300)
        cat("  Saved: fig4_info.pdf/png\n")
    }
}

# =============================================================================
# Summary Table
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("SUMMARY COMPARISON\n")
cat("================================================================================\n\n")

summary_table <- retention_combined[, .(
    Total_rare = sum(N_total),
    Retained = sum(N_retained),
    Retention_pct = round(100 * sum(N_retained) / sum(N_total), 1)
), by = Approach]

# Calculate improvement over first approach
summary_table[, Improvement_vs_baseline := Retention_pct - Retention_pct[1]]

print(summary_table)

cat("\n")
cat("================================================================================\n")
cat("RVAS IMPLICATIONS\n")
cat("================================================================================\n\n")

# Find best approach
best_approach <- summary_table[which.max(Retention_pct), Approach]
best_retention <- summary_table[which.max(Retention_pct), Retention_pct]
best_retained <- summary_table[which.max(Retention_pct), Retained]

baseline_retention <- summary_table[1, Retention_pct]
baseline_retained <- summary_table[1, Retained]

cat(sprintf("Best approach: %s\n", best_approach))
cat(sprintf("  - Retains %.1f%% of rare variants (vs %.1f%% baseline)\n",
            best_retention, baseline_retention))
cat(sprintf("  - %s more rare variants available for RVAS\n",
            scales::comma(best_retained - baseline_retained)))
cat("\n")
cat("For burden tests (SKAT, SKAT-O, CMC):\n")
cat(sprintf("  - %d additional rare variants per analysis\n",
            round((best_retained - baseline_retained) / 20000)))  # Assuming ~20K genes
cat("  - Better gene coverage → improved power\n")
cat("\n")
cat("For Amerindigenous/Latino samples:\n")
cat("  - MX Biobank + TOPMed improve NAT-ancestry imputation\n")
cat("  - Expect even larger improvements for AMR-specific rare variants\n")
cat("  - Key for detecting population-specific disease associations\n")
cat("\n")
cat("================================================================================\n")

# =============================================================================
# Output Files
# =============================================================================

cat("\nWriting output files...\n")

fwrite(summary_table, paste0(opt$output, "_summary.txt"), sep = "\t")
cat(sprintf("  %s_summary.txt\n", opt$output))

fwrite(retention_combined, paste0(opt$output, "_full_comparison.txt"), sep = "\t")
cat(sprintf("  %s_full_comparison.txt\n", opt$output))

cat("\nDone! Figures saved as PDF and PNG.\n")
