#!/usr/bin/env Rscript
################################################################################
# compare_approaches.R
#
# Compare concordance metrics across all imputation approaches.
# Generates publication-ready figures and summary tables.
#
# Usage:
#   Rscript compare_approaches.R \
#       --results-dir benchmarking/results/ \
#       --output comparison_report
#
################################################################################

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(cowplot)
    library(optparse)
    library(scales)
})

# =============================================================================
# Parse Arguments
# =============================================================================

option_list <- list(
    make_option(c("-r", "--results-dir"), type = "character", default = "benchmarking/results",
                help = "Directory containing approach results [default: %default]"),
    make_option(c("-o", "--output"), type = "character", default = "comparison_report",
                help = "Output prefix [default: %default]"),
    make_option(c("--approaches"), type = "character",
                default = "approach_a_michigan,approach_a_topmed,approach_a_allofus,approach_b_topmed,approach_c_topmed,approach_d_topmed,our_pipeline_topmed,our_pipeline_allofus",
                help = "Comma-separated list of approaches to compare"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Verbose output")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

approaches <- strsplit(opt$approaches, ",")[[1]]

log_msg <- function(msg) {
    if (opt$verbose) cat(sprintf("[%s] %s\n", Sys.time(), msg))
}

# =============================================================================
# Load Results
# =============================================================================

log_msg("Loading results from all approaches...")

all_results <- list()
all_maf_results <- list()
all_ancestry_results <- list()

for (approach in approaches) {
    result_dir <- file.path(opt$`results-dir`, approach)

    # Check for concordance results
    overall_file <- file.path(result_dir, "concordance_results_overall.txt")
    maf_file <- file.path(result_dir, "concordance_results_by_maf.txt")
    ancestry_file <- file.path(result_dir, "concordance_results_by_ancestry.txt")

    if (file.exists(overall_file)) {
        log_msg(sprintf("  Loading %s...", approach))

        overall <- fread(overall_file)
        overall[, Approach := approach]
        all_results[[approach]] <- overall

        if (file.exists(maf_file)) {
            maf <- fread(maf_file)
            maf[, Approach := approach]
            all_maf_results[[approach]] <- maf
        }

        if (file.exists(ancestry_file)) {
            anc <- fread(ancestry_file)
            anc[, Approach := approach]
            all_ancestry_results[[approach]] <- anc
        }
    } else {
        log_msg(sprintf("  WARNING: Results not found for %s", approach))
    }
}

# Combine results
results <- rbindlist(all_results, fill = TRUE)
maf_results <- rbindlist(all_maf_results, fill = TRUE)
ancestry_results <- rbindlist(all_ancestry_results, fill = TRUE)

# Parse approach names for better labels
results[, `:=`(
    QC_Timing = fifelse(grepl("approach_[ab]", Approach), "Before",
                        fifelse(grepl("approach_[cd]", Approach), "After", "Hybrid")),
    Server = fifelse(grepl("michigan", Approach), "Michigan/1KG",
                     fifelse(grepl("topmed", Approach), "TOPMed",
                             fifelse(grepl("allofus", Approach), "All of Us", "Unknown"))),
    Pipeline = fifelse(grepl("our_pipeline", Approach), "Our Pipeline", "Traditional")
)]

maf_results[, `:=`(
    QC_Timing = fifelse(grepl("approach_[ab]", Approach), "Before",
                        fifelse(grepl("approach_[cd]", Approach), "After", "Hybrid")),
    Server = fifelse(grepl("michigan", Approach), "Michigan/1KG",
                     fifelse(grepl("topmed", Approach), "TOPMed",
                             fifelse(grepl("allofus", Approach), "All of Us", "Unknown"))),
    Pipeline = fifelse(grepl("our_pipeline", Approach), "Our Pipeline", "Traditional")
)]

# =============================================================================
# Generate Comparison Figures
# =============================================================================

log_msg("Generating comparison figures...")

# Color palette
approach_colors <- c(
    "approach_a_michigan" = "#E41A1C",
    "approach_a_topmed" = "#377EB8",
    "approach_a_allofus" = "#4DAF4A",
    "approach_b_topmed" = "#984EA3",
    "approach_c_topmed" = "#FF7F00",
    "approach_d_topmed" = "#FFFF33",
    "our_pipeline_topmed" = "#A65628",
    "our_pipeline_allofus" = "#F781BF"
)

server_colors <- c(
    "Michigan/1KG" = "#E41A1C",
    "TOPMed" = "#377EB8",
    "All of Us" = "#4DAF4A"
)

# Figure 1: Overall concordance by approach
if (nrow(results) > 0) {
    p1 <- ggplot(results, aes(x = reorder(Approach, Mean_Concordance), y = Mean_Concordance, fill = Server)) +
        geom_bar(stat = "identity") +
        geom_errorbar(aes(ymin = Mean_Concordance - SD_Concordance,
                          ymax = Mean_Concordance + SD_Concordance),
                      width = 0.2) +
        coord_flip() +
        scale_fill_manual(values = server_colors) +
        labs(title = "Genotype Concordance by Approach",
             x = "", y = "Mean Concordance") +
        theme_bw() +
        theme(legend.position = "bottom")

    ggsave(paste0(opt$output, "_fig1_overall_concordance.pdf"), p1, width = 10, height = 6)
    log_msg("  Saved fig1_overall_concordance.pdf")
}

# Figure 2: Concordance by MAF bin (heatmap)
if (nrow(maf_results) > 0) {
    p2 <- ggplot(maf_results, aes(x = MAF_BIN, y = Approach, fill = Mean_Concordance)) +
        geom_tile() +
        geom_text(aes(label = sprintf("%.3f", Mean_Concordance)), size = 3) +
        scale_fill_gradient2(low = "red", mid = "yellow", high = "green",
                             midpoint = 0.9, limits = c(0.7, 1)) +
        labs(title = "Concordance by MAF Bin and Approach",
             x = "MAF Bin", y = "", fill = "Concordance") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(paste0(opt$output, "_fig2_maf_heatmap.pdf"), p2, width = 12, height = 8)
    log_msg("  Saved fig2_maf_heatmap.pdf")
}

# Figure 3: Dosage R² by MAF bin
if (nrow(maf_results) > 0) {
    p3 <- ggplot(maf_results, aes(x = MAF_BIN, y = Mean_Dosage_R2, color = Approach, group = Approach)) +
        geom_line(size = 1) +
        geom_point(size = 2) +
        facet_wrap(~Server) +
        labs(title = "Dosage R² by MAF Bin",
             x = "MAF Bin", y = "Mean Dosage R²") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom")

    ggsave(paste0(opt$output, "_fig3_dosage_r2.pdf"), p3, width = 12, height = 6)
    log_msg("  Saved fig3_dosage_r2.pdf")
}

# Figure 4: Concordance by ancestry (if available)
if (nrow(ancestry_results) > 0) {
    # Focus on our pipeline vs best traditional approach
    compare_approaches <- c("approach_a_topmed", "our_pipeline_topmed")
    ancestry_compare <- ancestry_results[Approach %in% compare_approaches]

    if (nrow(ancestry_compare) > 0) {
        p4 <- ggplot(ancestry_compare, aes(x = Ancestry, y = Mean_Concordance, fill = Approach)) +
            geom_bar(stat = "identity", position = "dodge") +
            geom_errorbar(aes(ymin = Mean_Concordance - SD_Concordance,
                              ymax = Mean_Concordance + SD_Concordance),
                          position = position_dodge(0.9), width = 0.2) +
            facet_wrap(~MAF_BIN) +
            labs(title = "Concordance by Ancestry: Our Pipeline vs Traditional",
                 x = "Ancestry Group", y = "Mean Concordance") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "bottom")

        ggsave(paste0(opt$output, "_fig4_ancestry_comparison.pdf"), p4, width = 14, height = 10)
        log_msg("  Saved fig4_ancestry_comparison.pdf")
    }
}

# Figure 5: Summary panel
if (nrow(results) > 0 && nrow(maf_results) > 0) {
    # Our pipeline vs all others for rare variants
    rare_maf <- maf_results[MAF_BIN == "0-0.5%"]

    if (nrow(rare_maf) > 0) {
        p5a <- ggplot(rare_maf, aes(x = reorder(Approach, Mean_Concordance), y = Mean_Concordance, fill = Pipeline)) +
            geom_bar(stat = "identity") +
            coord_flip() +
            scale_fill_manual(values = c("Our Pipeline" = "#2166AC", "Traditional" = "#B2182B")) +
            labs(title = "Rare Variant Concordance (MAF < 0.5%)",
                 x = "", y = "Concordance") +
            theme_bw()

        p5b <- ggplot(rare_maf, aes(x = reorder(Approach, Mean_Dosage_R2), y = Mean_Dosage_R2, fill = Pipeline)) +
            geom_bar(stat = "identity") +
            coord_flip() +
            scale_fill_manual(values = c("Our Pipeline" = "#2166AC", "Traditional" = "#B2182B")) +
            labs(title = "Rare Variant Dosage R² (MAF < 0.5%)",
                 x = "", y = "Dosage R²") +
            theme_bw()

        p5_combined <- plot_grid(p5a, p5b, ncol = 2, labels = c("A", "B"))
        ggsave(paste0(opt$output, "_fig5_rare_variants.pdf"), p5_combined, width = 14, height = 6)
        log_msg("  Saved fig5_rare_variants.pdf")
    }
}

# =============================================================================
# Generate Summary Tables
# =============================================================================

log_msg("Generating summary tables...")

# Table 1: Overall comparison
summary_table <- results[, .(
    Approach, Server, QC_Timing, Pipeline,
    N_variants, N_samples,
    Concordance = sprintf("%.4f (%.4f)", Mean_Concordance, SD_Concordance),
    NRC = sprintf("%.4f (%.4f)", Mean_NRC, SD_NRC),
    Dosage_R2 = sprintf("%.4f (%.4f)", Mean_Dosage_R2, SD_Dosage_R2)
)]

fwrite(summary_table, paste0(opt$output, "_table1_overall.txt"), sep = "\t")
log_msg("  Saved table1_overall.txt")

# Table 2: By MAF bin
if (nrow(maf_results) > 0) {
    maf_wide <- dcast(maf_results,
                      MAF_BIN ~ Approach,
                      value.var = "Mean_Concordance")
    fwrite(maf_wide, paste0(opt$output, "_table2_maf.txt"), sep = "\t")
    log_msg("  Saved table2_maf.txt")
}

# Table 3: By ancestry (if available)
if (nrow(ancestry_results) > 0) {
    ancestry_wide <- dcast(ancestry_results,
                           Ancestry + MAF_BIN ~ Approach,
                           value.var = "Mean_Concordance")
    fwrite(ancestry_wide, paste0(opt$output, "_table3_ancestry.txt"), sep = "\t")
    log_msg("  Saved table3_ancestry.txt")
}

# =============================================================================
# Statistical Tests
# =============================================================================

log_msg("Performing statistical comparisons...")

# Compare our pipeline to traditional approaches
if ("our_pipeline_topmed" %in% results$Approach && "approach_a_topmed" %in% results$Approach) {

    our_result <- results[Approach == "our_pipeline_topmed"]
    trad_result <- results[Approach == "approach_a_topmed"]

    improvement <- data.table(
        Metric = c("Concordance", "NRC", "Dosage_R2"),
        Our_Pipeline = c(our_result$Mean_Concordance, our_result$Mean_NRC, our_result$Mean_Dosage_R2),
        Traditional = c(trad_result$Mean_Concordance, trad_result$Mean_NRC, trad_result$Mean_Dosage_R2)
    )
    improvement[, Difference := Our_Pipeline - Traditional]
    improvement[, Pct_Improvement := (Difference / Traditional) * 100]

    fwrite(improvement, paste0(opt$output, "_improvement.txt"), sep = "\t")
    log_msg("  Saved improvement.txt")
}

# =============================================================================
# Print Summary
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("APPROACH COMPARISON SUMMARY\n")
cat("================================================================================\n")
cat("\n")
cat("Approaches compared:", length(unique(results$Approach)), "\n")
cat("\n")
cat("Overall Results:\n")
print(summary_table)
cat("\n")

if (exists("improvement")) {
    cat("Our Pipeline vs Traditional (Approach A + TOPMed):\n")
    print(improvement)
    cat("\n")
}

cat("================================================================================\n")
cat(sprintf("Figures saved to: %s_fig*.pdf\n", opt$output))
cat(sprintf("Tables saved to: %s_table*.txt\n", opt$output))
cat("================================================================================\n")

log_msg("Done!")
