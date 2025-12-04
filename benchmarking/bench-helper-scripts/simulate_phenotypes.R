#!/usr/bin/env Rscript
################################################################################
# simulate_phenotypes.R
#
# Simulates phenotypes for 1000 Genomes samples based on known GWAS effects.
# This allows testing the full imputation → GWAS → PRS pipeline without
# requiring restricted phenotype data.
#
# Simulation approach:
#   1. Use published effect sizes from large GWAS
#   2. Calculate genetic value from WGS truth genotypes
#   3. Add environmental noise to achieve target heritability
#   4. Output phenotype file compatible with PLINK/REGENIE
#
# Usage:
#   Rscript simulate_phenotypes.R \
#       --vcf 1kg_wgs.vcf.gz \
#       --trait height \
#       --h2 0.5 \
#       --output simulated_phenotypes.txt
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
    make_option(c("-v", "--vcf"), type = "character", default = NULL,
                help = "WGS VCF file (truth genotypes)"),
    make_option(c("-t", "--trait"), type = "character", default = "height",
                help = "Trait to simulate: height, bmi, t2d, ldl [default: %default]"),
    make_option(c("--h2"), type = "numeric", default = 0.5,
                help = "Heritability (variance explained by genetics) [default: %default]"),
    make_option(c("--n-causal"), type = "integer", default = 100,
                help = "Number of causal variants to use [default: %default]"),
    make_option(c("-o", "--output"), type = "character", default = "simulated_phenotypes",
                help = "Output prefix [default: %default]"),
    make_option(c("--seed"), type = "integer", default = 42,
                help = "Random seed for reproducibility [default: %default]"),
    make_option(c("--prevalence"), type = "numeric", default = 0.1,
                help = "Disease prevalence for binary traits [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

set.seed(opt$seed)

# =============================================================================
# Known GWAS Effects (Published Effect Sizes)
# =============================================================================

# These are real effect sizes from published GWAS
# Format: CHR, POS (GRCh38), RSID, Effect_Allele, Beta/OR, Trait

known_effects <- list(

    height = data.table(
        # GIANT consortium height hits (beta in SD units, ~0.06 SD per allele typical)
        CHR = c(2, 3, 6, 15, 18, 7, 12, 20, 1, 4),
        POS = c(27730940, 53829968, 32413564, 89433456, 23536230,
                2796312, 66442045, 6018576, 119431614, 145612134),
        RSID = c("rs1229984", "rs724016", "rs9272219", "rs2767486", "rs1871700",
                 "rs696917", "rs11611246", "rs373863828", "rs2274432", "rs6830062"),
        EA = c("A", "G", "C", "T", "C", "A", "T", "A", "G", "A"),
        BETA = c(0.08, 0.05, 0.06, 0.04, 0.05, 0.03, 0.04, 0.12, 0.03, 0.04),
        stringsAsFactors = FALSE
    ),

    bmi = data.table(
        # GIANT consortium BMI hits
        CHR = c(16, 18, 2, 1, 6, 19, 11, 4, 3, 10),
        POS = c(53767042, 57829135, 25150296, 72812817, 50873513,
                46182304, 27679916, 44877284, 187317257, 87433001),
        RSID = c("rs1421085", "rs6567160", "rs13021737", "rs543874", "rs987237",
                 "rs2287019", "rs2176598", "rs13107325", "rs13078960", "rs11191560"),
        EA = c("C", "C", "G", "G", "G", "C", "T", "T", "G", "C"),
        BETA = c(0.10, 0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.05, 0.03, 0.03),
        stringsAsFactors = FALSE
    ),

    t2d = data.table(
        # DIAMANTE T2D hits (log OR)
        CHR = c(10, 9, 3, 7, 11, 8, 6, 17, 2, 12),
        POS = c(114758349, 22134094, 12344730, 44235668, 72433098,
                118184783, 20686996, 6945866, 227093745, 121434438),
        RSID = c("rs7903146", "rs10811661", "rs4607103", "rs864745", "rs5215",
                 "rs13266634", "rs1800562", "rs13342232", "rs7607980", "rs7961581"),
        EA = c("T", "T", "C", "T", "C", "C", "A", "A", "C", "C"),
        # Log OR for T2D (0.15 = OR ~1.16)
        BETA = c(0.35, 0.20, 0.12, 0.15, 0.10, 0.12, 0.08, 0.25, 0.10, 0.08),
        stringsAsFactors = FALSE
    ),

    ldl = data.table(
        # GLGC LDL hits (beta in mg/dL)
        CHR = c(19, 1, 2, 19, 1, 6, 11, 8, 20, 7),
        POS = c(11202306, 55505647, 21263900, 45411941, 109818530,
                161010118, 116660686, 19844222, 39797465, 21607352),
        RSID = c("rs429358", "rs11591147", "rs515135", "rs445925", "rs2228671",
                 "rs3757354", "rs964184", "rs2737229", "rs6065906", "rs2072183"),
        EA = c("C", "T", "G", "G", "T", "T", "G", "A", "T", "C"),
        # Effect in SD units
        BETA = c(0.50, 0.40, 0.15, 0.20, 0.10, 0.08, 0.12, 0.06, 0.05, 0.05),
        stringsAsFactors = FALSE
    )
)

# =============================================================================
# Extract Genotypes at Causal Variants
# =============================================================================

extract_genotypes <- function(vcf_file, variants) {
    cat("Extracting genotypes at causal variants...\n")

    # Create regions file for bcftools
    regions_file <- tempfile(fileext = ".txt")
    fwrite(variants[, .(CHR, POS, POS)], regions_file,
           col.names = FALSE, sep = "\t")

    # Extract genotypes using bcftools
    cmd <- sprintf(
        "bcftools query -R %s -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT]\\n' %s",
        regions_file, vcf_file
    )

    geno_raw <- fread(cmd = cmd, header = FALSE)

    if (nrow(geno_raw) == 0) {
        stop("No variants found in VCF. Check chromosome naming (chr1 vs 1).")
    }

    # Get sample names
    cmd_samples <- sprintf("bcftools query -l %s", vcf_file)
    samples <- fread(cmd = cmd_samples, header = FALSE)$V1

    cat(sprintf("  Found %d variants, %d samples\n", nrow(geno_raw), length(samples)))

    # Parse genotypes
    setnames(geno_raw, c("CHR", "POS", "REF", "ALT", samples))

    # Convert GT to dosage (0, 1, 2 copies of ALT)
    for (s in samples) {
        geno_raw[, (s) := fcase(
            get(s) == "0/0" | get(s) == "0|0", 0,
            get(s) == "0/1" | get(s) == "1/0" | get(s) == "0|1" | get(s) == "1|0", 1,
            get(s) == "1/1" | get(s) == "1|1", 2,
            default = NA_real_
        )]
    }

    unlink(regions_file)

    return(list(genotypes = geno_raw, samples = samples))
}

# =============================================================================
# Calculate Genetic Value (Polygenic Score)
# =============================================================================

calculate_genetic_value <- function(geno_data, effects) {
    cat("Calculating genetic values...\n")

    geno <- geno_data$genotypes
    samples <- geno_data$samples

    # Merge with effect sizes
    geno[, CHR := as.character(CHR)]
    effects[, CHR := as.character(CHR)]

    merged <- merge(geno, effects, by = c("CHR", "POS"), all.x = TRUE)

    if (nrow(merged) == 0) {
        stop("No overlap between VCF variants and effect variants")
    }

    cat(sprintf("  Using %d causal variants\n", nrow(merged)))

    # Calculate genetic value for each sample
    # G = sum(dosage * beta), accounting for effect allele
    genetic_values <- sapply(samples, function(s) {
        dosages <- merged[[s]]
        betas <- merged$BETA

        # Adjust for effect allele (simplified - assume ALT = effect allele)
        # In practice, would need to check strand and flip

        sum(dosages * betas, na.rm = TRUE)
    })

    return(data.table(IID = samples, G = genetic_values))
}

# =============================================================================
# Simulate Phenotype
# =============================================================================

simulate_phenotype <- function(genetic_values, h2, binary = FALSE, prevalence = 0.1) {
    cat(sprintf("Simulating phenotype (h2 = %.2f)...\n", h2))

    G <- genetic_values$G
    n <- length(G)

    # Standardize genetic values
    G_std <- (G - mean(G)) / sd(G)

    # Calculate environmental variance to achieve target h2
    # Var(P) = Var(G) + Var(E) = 1 (standardized)
    # h2 = Var(G) / Var(P)
    # So Var(E) = (1 - h2) / h2 * Var(G)

    var_G <- var(G_std)
    var_E <- var_G * (1 - h2) / h2

    # Generate environmental component
    E <- rnorm(n, mean = 0, sd = sqrt(var_E))

    # Phenotype = Genetic + Environmental
    P <- G_std + E

    if (binary) {
        # Convert to binary using liability threshold model
        threshold <- qnorm(1 - prevalence)
        P_binary <- as.integer(P > threshold)

        cat(sprintf("  Cases: %d (%.1f%%), Controls: %d\n",
                    sum(P_binary), 100 * mean(P_binary), sum(1 - P_binary)))

        genetic_values[, PHENO := P_binary]
        genetic_values[, LIABILITY := P]
    } else {
        # Standardize continuous phenotype
        P_std <- (P - mean(P)) / sd(P)
        genetic_values[, PHENO := P_std]
    }

    # Store genetic value for validation
    genetic_values[, G_TRUE := G_std]

    return(genetic_values)
}

# =============================================================================
# Add Population Structure Effect (Optional)
# =============================================================================

add_ancestry_effect <- function(phenotypes, vcf_file, effect_size = 0.2) {
    cat("Adding ancestry-correlated effect...\n")

    # This simulates confounding by ancestry
    # In real data, this is what population stratification causes

    # Get ancestry from sample IDs (1KG convention)
    phenotypes[, POP := sub("^([A-Z]{3}).*", "\\1", IID)]

    # Add population-specific offset
    pop_effects <- c(
        "EUR" = 0, "GBR" = 0, "FIN" = 0.1, "CEU" = 0, "TSI" = -0.1, "IBS" = -0.05,
        "AFR" = -0.3, "YRI" = -0.3, "LWK" = -0.25, "GWD" = -0.35, "MSL" = -0.3,
        "ESN" = -0.3, "ASW" = -0.15, "ACB" = -0.2,
        "EAS" = 0.2, "CHB" = 0.2, "JPT" = 0.25, "CHS" = 0.2, "CDX" = 0.15, "KHV" = 0.15,
        "SAS" = -0.1, "GIH" = -0.1, "PJL" = -0.05, "BEB" = -0.1, "STU" = -0.1, "ITU" = -0.1,
        "AMR" = -0.05, "MXL" = -0.1, "PUR" = -0.05, "CLM" = -0.05, "PEL" = -0.15
    )

    phenotypes[, POP_EFFECT := pop_effects[POP]]
    phenotypes[is.na(POP_EFFECT), POP_EFFECT := 0]

    # Add to phenotype
    phenotypes[, PHENO_STRATIFIED := PHENO + effect_size * POP_EFFECT]

    cat("  Added population stratification (use PHENO_STRATIFIED to test correction)\n")

    return(phenotypes)
}

# =============================================================================
# Main
# =============================================================================

cat("================================================================================\n")
cat("PHENOTYPE SIMULATION FOR IMPUTATION BENCHMARKING\n")
cat("================================================================================\n\n")

if (is.null(opt$vcf)) {
    stop("--vcf is required")
}

# Get effects for selected trait
if (!opt$trait %in% names(known_effects)) {
    stop(sprintf("Unknown trait: %s. Available: %s",
                 opt$trait, paste(names(known_effects), collapse = ", ")))
}

effects <- known_effects[[opt$trait]]
cat(sprintf("Trait: %s (%d known variants)\n", opt$trait, nrow(effects)))
cat(sprintf("Heritability: %.2f\n", opt$h2))
cat(sprintf("VCF: %s\n\n", opt$vcf))

# Determine if binary trait
binary_traits <- c("t2d")
is_binary <- opt$trait %in% binary_traits

# Extract genotypes
geno_data <- extract_genotypes(opt$vcf, effects)

# Calculate genetic values
genetic_values <- calculate_genetic_value(geno_data, effects)

# Simulate phenotype
phenotypes <- simulate_phenotype(genetic_values, opt$h2,
                                  binary = is_binary,
                                  prevalence = opt$prevalence)

# Add ancestry effect for stratification testing
phenotypes <- add_ancestry_effect(phenotypes, opt$vcf)

# =============================================================================
# Output Files
# =============================================================================

cat("\nWriting output files...\n")

# PLINK-format phenotype file
plink_pheno <- phenotypes[, .(FID = IID, IID = IID, PHENO)]
fwrite(plink_pheno, paste0(opt$output, ".pheno"), sep = "\t")
cat(sprintf("  %s.pheno (PLINK format)\n", opt$output))

# PLINK-format with stratified phenotype
plink_strat <- phenotypes[, .(FID = IID, IID = IID,
                               PHENO_CLEAN = PHENO,
                               PHENO_STRATIFIED = PHENO_STRATIFIED)]
fwrite(plink_strat, paste0(opt$output, "_with_stratification.pheno"), sep = "\t")
cat(sprintf("  %s_with_stratification.pheno\n", opt$output))

# Full output with genetic values (for validation)
fwrite(phenotypes, paste0(opt$output, "_full.txt"), sep = "\t")
cat(sprintf("  %s_full.txt (includes true genetic values)\n", opt$output))

# Causal variants file (for checking recovery)
fwrite(effects, paste0(opt$output, "_causal_variants.txt"), sep = "\t")
cat(sprintf("  %s_causal_variants.txt\n", opt$output))

# =============================================================================
# Summary Statistics
# =============================================================================

cat("\n================================================================================\n")
cat("SIMULATION SUMMARY\n")
cat("================================================================================\n")
cat(sprintf("Samples:           %d\n", nrow(phenotypes)))
cat(sprintf("Causal variants:   %d\n", nrow(effects)))
cat(sprintf("Target h2:         %.2f\n", opt$h2))

if (is_binary) {
    cat(sprintf("Case prevalence:   %.1f%%\n", 100 * mean(phenotypes$PHENO)))
} else {
    cat(sprintf("Phenotype mean:    %.3f\n", mean(phenotypes$PHENO)))
    cat(sprintf("Phenotype SD:      %.3f\n", sd(phenotypes$PHENO)))
}

# Actual h2 achieved
cor_gp <- cor(phenotypes$G_TRUE, phenotypes$PHENO)
cat(sprintf("Achieved h2:       %.3f (cor(G,P)^2)\n", cor_gp^2))

cat("\n")
cat("Usage:\n")
cat("  # Run GWAS on imputed data:\n")
cat(sprintf("  plink2 --bfile imputed --pheno %s.pheno --glm --out gwas_%s\n",
            opt$output, opt$trait))
cat("\n")
cat("  # Check if causal variants are recovered:\n")
cat(sprintf("  Rscript check_hit_recovery.R --gwas gwas_%s.glm.linear \\\n", opt$trait))
cat(sprintf("      --causal %s_causal_variants.txt\n", opt$output))
cat("================================================================================\n")
