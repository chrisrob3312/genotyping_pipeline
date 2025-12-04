#!/usr/bin/env Rscript
################################################################################
# summarize_benchmark_timings.R
#
# Generates timing summary tables and figures from benchmark matrix results.
#
# Usage:
#   Rscript summarize_benchmark_timings.R \
#       --timings benchmark_results/timing/all_timings.tsv \
#       --counts benchmark_results/final_counts.tsv \
#       --output benchmark_timing_report
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
    make_option(c("-t", "--timings"), type = "character", default = NULL,
                help = "Timing data file (TSV)"),
    make_option(c("-c", "--counts"), type = "character", default = NULL,
                help = "Final variant/sample counts file"),
    make_option(c("-o", "--output"), type = "character", default = "timing_report",
                help = "Output prefix [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# =============================================================================
# Theme
# =============================================================================

theme_benchmark <- function() {
    theme_minimal() +
    theme(
        text = element_text(size = 11),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold")
    )
}

# Color schemes
approach_colors <- c(
    "A" = "#E69F00", "B" = "#F0E442",
    "C" = "#56B4E9", "D" = "#009E73",
    "E" = "#0072B2", "F" = "#D55E00"
)

server_colors <- c(
    "topmed" = "#1f77b4",
    "allofus" = "#2ca02c",
    "michigan_hrc" = "#ff7f0e",
    "michigan_1kg" = "#d62728"
)

# =============================================================================
# Load Data
# =============================================================================

cat("================================================================================\n")
cat("BENCHMARK TIMING SUMMARY\n")
cat("================================================================================\n\n")

if (!is.null(opt$timings) && file.exists(opt$timings)) {
    timings <- fread(opt$timings)
    cat(sprintf("Loaded %d timing records\n", nrow(timings)))
} else {
    stop("Timing file not found")
}

if (!is.null(opt$counts) && file.exists(opt$counts)) {
    counts <- fread(opt$counts, header = FALSE,
                    col.names = c("approach", "server", "n_variants", "n_samples"))
    cat(sprintf("Loaded %d final count records\n", nrow(counts)))
} else {
    counts <- NULL
    cat("No counts file provided\n")
}

# =============================================================================
# Calculate Total Times
# =============================================================================

cat("\n=== Total Pipeline Times ===\n\n")

total_times <- timings[, .(
    total_seconds = sum(duration_seconds)
), by = .(approach, server)]

total_times[, total_human := sprintf("%.1f min", total_seconds / 60)]

# Wide format for comparison
total_wide <- dcast(total_times, approach ~ server, value.var = "total_seconds")

cat("Total runtime (seconds) by approach × server:\n")
print(total_wide)

# =============================================================================
# Figure 1: Total Runtime Heatmap
# =============================================================================

cat("\n=== Generating Figures ===\n")

fig1 <- ggplot(total_times, aes(x = server, y = approach, fill = total_seconds / 60)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.0f", total_seconds / 60)), color = "white", size = 4) +
    scale_fill_viridis_c(option = "plasma", name = "Minutes") +
    labs(
        title = "Total Pipeline Runtime (minutes)",
        subtitle = "Full benchmark matrix: 6 approaches × 4 servers",
        x = "Imputation Server",
        y = "Approach"
    ) +
    theme_benchmark() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(paste0(opt$output, "_fig1_heatmap.pdf"), fig1, width = 10, height = 6)
ggsave(paste0(opt$output, "_fig1_heatmap.png"), fig1, width = 10, height = 6, dpi = 300)
cat("  Fig 1: Runtime heatmap saved\n")

# =============================================================================
# Figure 2: Runtime by Step
# =============================================================================

# Aggregate by step
step_times <- timings[, .(
    mean_seconds = mean(duration_seconds),
    sd_seconds = sd(duration_seconds)
), by = .(approach, step)]

fig2 <- ggplot(timings, aes(x = step, y = duration_seconds / 60, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = approach_colors) +
    labs(
        title = "Runtime by Pipeline Step",
        x = "Step",
        y = "Duration (minutes)",
        fill = "Approach"
    ) +
    theme_benchmark()

ggsave(paste0(opt$output, "_fig2_by_step.pdf"), fig2, width = 12, height = 6)
ggsave(paste0(opt$output, "_fig2_by_step.png"), fig2, width = 12, height = 6, dpi = 300)
cat("  Fig 2: By-step comparison saved\n")

# =============================================================================
# Figure 3: Server Comparison
# =============================================================================

fig3 <- ggplot(total_times, aes(x = approach, y = total_seconds / 60, fill = server)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = server_colors) +
    labs(
        title = "Total Runtime by Approach and Server",
        x = "Approach",
        y = "Total Time (minutes)",
        fill = "Server"
    ) +
    theme_benchmark()

ggsave(paste0(opt$output, "_fig3_server_comparison.pdf"), fig3, width = 10, height = 6)
ggsave(paste0(opt$output, "_fig3_server_comparison.png"), fig3, width = 10, height = 6, dpi = 300)
cat("  Fig 3: Server comparison saved\n")

# =============================================================================
# Figure 4: Traditional vs Our Pipeline
# =============================================================================

total_times[, pipeline_type := ifelse(approach %in% c("A", "B", "C", "D"),
                                       "Traditional (QC-before)",
                                       "Our Pipeline (QC-after)")]

fig4 <- ggplot(total_times, aes(x = pipeline_type, y = total_seconds / 60, fill = server)) +
    geom_boxplot() +
    scale_fill_manual(values = server_colors) +
    labs(
        title = "Traditional vs Our Pipeline: Runtime Comparison",
        x = "",
        y = "Total Time (minutes)",
        fill = "Server"
    ) +
    theme_benchmark() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(paste0(opt$output, "_fig4_trad_vs_ours.pdf"), fig4, width = 8, height = 6)
ggsave(paste0(opt$output, "_fig4_trad_vs_ours.png"), fig4, width = 8, height = 6, dpi = 300)
cat("  Fig 4: Traditional vs Our Pipeline saved\n")

# =============================================================================
# Figure 5: Variant Yield (if counts available)
# =============================================================================

if (!is.null(counts)) {
    fig5 <- ggplot(counts, aes(x = approach, y = n_variants / 1e6, fill = server)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
        scale_fill_manual(values = server_colors) +
        labs(
            title = "Final Variant Count by Approach and Server",
            subtitle = "Higher is better (more variants retained)",
            x = "Approach",
            y = "Variants (millions)",
            fill = "Server"
        ) +
        theme_benchmark()

    ggsave(paste0(opt$output, "_fig5_variant_yield.pdf"), fig5, width = 10, height = 6)
    ggsave(paste0(opt$output, "_fig5_variant_yield.png"), fig5, width = 10, height = 6, dpi = 300)
    cat("  Fig 5: Variant yield saved\n")

    # Efficiency metric: variants per minute
    efficiency <- merge(total_times, counts, by = c("approach", "server"))
    efficiency[, variants_per_minute := n_variants / (total_seconds / 60)]

    fig6 <- ggplot(efficiency, aes(x = approach, y = variants_per_minute / 1e6, fill = server)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
        scale_fill_manual(values = server_colors) +
        labs(
            title = "Efficiency: Variants Retained per Minute of Compute",
            subtitle = "Higher is better",
            x = "Approach",
            y = "Million Variants / Minute",
            fill = "Server"
        ) +
        theme_benchmark()

    ggsave(paste0(opt$output, "_fig6_efficiency.pdf"), fig6, width = 10, height = 6)
    ggsave(paste0(opt$output, "_fig6_efficiency.png"), fig6, width = 10, height = 6, dpi = 300)
    cat("  Fig 6: Efficiency metric saved\n")
}

# =============================================================================
# Summary Tables
# =============================================================================

cat("\n=== Summary Statistics ===\n\n")

# Best performers
cat("Fastest approach by server:\n")
fastest <- total_times[, .SD[which.min(total_seconds)], by = server]
print(fastest[, .(server, approach, total_human)])

cat("\nHighest variant yield by server:\n")
if (!is.null(counts)) {
    most_variants <- counts[, .SD[which.max(n_variants)], by = server]
    print(most_variants)
}

# Write summary table
fwrite(total_times, paste0(opt$output, "_summary.tsv"), sep = "\t")
cat(sprintf("\nSummary table saved: %s_summary.tsv\n", opt$output))

cat("\n================================================================================\n")
cat("TIMING REPORT COMPLETE\n")
cat("================================================================================\n")
cat(sprintf("\nFigures saved with prefix: %s\n", opt$output))
cat("  - fig1: Runtime heatmap\n")
cat("  - fig2: By-step comparison\n")
cat("  - fig3: Server comparison\n")
cat("  - fig4: Traditional vs Our Pipeline\n")
if (!is.null(counts)) {
    cat("  - fig5: Variant yield\n")
    cat("  - fig6: Efficiency metric\n")
}
cat("\n")
