#!/usr/bin/env Rscript
# MagicalRsq-X filtering script
# Applies machine-learning based imputation quality filtering

library(MagicalRsq)  # or MagicalRsqX if available
library(data.table)
library(optparse)

# Parse command line arguments
option_list <- list(
  make_option("--info", type="character", help="Input .info file from imputation"),
  make_option("--vcf", type="character", help="Input VCF file"),
  make_option("--out", type="character", help="Output prefix"),
  make_option("--threshold", type="numeric", default=0.3, help="MagicalRsq threshold"),
  make_option("--use_maf_stratified", type="logical", default=FALSE, 
              help="Use MAF-stratified thresholds"),
  make_option("--common_thresh", type="numeric", default=0.3),
  make_option("--lowfreq_thresh", type="numeric", default=0.5),
  make_option("--rare_thresh", type="numeric", default=0.8),
  make_option("--model", type="character", help="Pre-trained MagicalRsq-X model path")
)

opt <- parse_args(OptionParser(option_list=option_list))

cat("=================================================\n")
cat("MagicalRsq-X Filtering\n")
cat("=================================================\n")
cat("Input INFO:", opt$info, "\n")
cat("Input VCF:", opt$vcf, "\n")
cat("Threshold:", opt$threshold, "\n")

if (opt$use_maf_stratified) {
  cat("Using MAF-stratified thresholds:\n")
  cat("  Common (MAF>0.05):", opt$common_thresh, "\n")
  cat("  Low-freq (0.01<MAF≤0.05):", opt$lowfreq_thresh, "\n")
  cat("  Rare (MAF≤0.01):", opt$rare_thresh, "\n")
}

# Read imputation info file
info <- fread(opt$info)

# Calculate MagicalRsq-X if model provided
if (!is.null(opt$model)) {
  cat("\nCalculating MagicalRsq-X using pre-trained model...\n")
  
  # Apply MagicalRsq-X model (function depends on package version)
  # This is a placeholder - adjust based on actual MagicalRsqX package API
  info$MagicalRsq <- predict_MagicalRsq(info, model=opt$model)
  
} else {
  # Fall back to standard R² from imputation server
  cat("\nWARNING: No MagicalRsq-X model provided\n")
  cat("Using standard R² from imputation INFO file\n")
  
  if ("R2" %in% names(info)) {
    info$MagicalRsq <- info$R2
  } else if ("Rsq" %in% names(info)) {
    info$MagicalRsq <- info$Rsq
  } else {
    stop("Cannot find R2 or Rsq column in INFO file")
  }
}

# Apply filtering
if (opt$use_maf_stratified) {
  # MAF-stratified filtering
  info[, pass := ifelse(MAF > 0.05, MagicalRsq >= opt$common_thresh,
                 ifelse(MAF > 0.01, MagicalRsq >= opt$lowfreq_thresh,
                                     MagicalRsq >= opt$rare_thresh))]
} else {
  # Single threshold
  info[, pass := MagicalRsq >= opt$threshold]
}

# Statistics
n_total <- nrow(info)
n_pass <- sum(info$pass)
n_fail <- n_total - n_pass
pct_pass <- round(100 * n_pass / n_total, 2)

cat("\n=================================================\n")
cat("Filtering Results:\n")
cat("=================================================\n")
cat("Total variants:", n_total, "\n")
cat("Passing filter:", n_pass, "(", pct_pass, "%)\n")
cat("Failing filter:", n_fail, "\n")

# Save passing variants list
pass_vars <- info[pass == TRUE, .(SNP)]
fwrite(pass_vars, paste0(opt$out, "_pass_variants.txt"), col.names=FALSE)

# Save statistics
stats <- data.table(
  Metric = c("Total_variants", "Passing_filter", "Failing_filter", "Percent_passing"),
  Value = c(n_total, n_pass, n_fail, pct_pass)
)
fwrite(stats, paste0(opt$out, "_magicalrsq_stats.txt"), sep="\t")

cat("\nOutput files created:\n")
cat("  ", paste0(opt$out, "_pass_variants.txt"), "\n")
cat("  ", paste0(opt$out, "_magicalrsq_stats.txt"), "\n")
cat("=================================================\n")#Placeholder for R script for MagicalRSq
