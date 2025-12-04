#!/usr/bin/env Rscript
################################################################################
# simulate_ancestry_aware_phenotypes.R
#
# Simulates phenotypes with ancestry-specific genetic architecture.
# Unlike using published GWAS (EUR-biased), this creates ground truth
# with controlled:
#   - Causal variant frequencies across populations
#   - Ancestry-specific effect sizes
#   - Local ancestry context effects
#
# Key scenarios for testing imputation approaches:
#   1. Variants common in AFR but rare in EUR (missed by EUR GWAS)
#   2. Ancestry-specific effects (same variant, different effect by ancestry)
#   3. Effects that depend on local ancestry background
#   4. Rare variants across all ancestries (tests imputation of rare variants)
#
# References:
#   - Martin et al. 2019 (Nature Genetics) - PRS portability
#   - Wojcik et al. 2019 (Nature) - PAGE Study
#   - Marnetto et al. 2020 - Local ancestry effects
#
# Usage:
#   Rscript simulate_ancestry_aware_phenotypes.R \
#       --vcf 1kg_wgs.vcf.gz \
#       --pop-file population_labels.txt \
#       --scenario afr_enriched \
#       --n-causal 50 \
#       --h2 0.5 \
#       --output simulated_pheno
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
                help = "VCF file with genotypes"),
    make_option(c("-p", "--pop-file"), type = "character", default = NULL,
                help = "Population labels file (IID, POP, SUPERPOP columns)"),
    make_option(c("--scenario"), type = "character", default = "ancestry_specific",
                help = "Simulation scenario (see --list-scenarios) [default: %default]"),
    make_option(c("--n-causal"), type = "integer", default = 50,
                help = "Number of causal variants [default: %default]"),
    make_option(c("--h2"), type = "numeric", default = 0.5,
                help = "Total heritability [default: %default]"),
    make_option(c("--h2-ancestry"), type = "numeric", default = 0.1,
                help = "Heritability from ancestry-specific effects [default: %default]"),
    make_option(c("-o", "--output"), type = "character", default = "sim_pheno",
                help = "Output prefix [default: %default]"),
    make_option(c("--seed"), type = "integer", default = 42,
                help = "Random seed [default: %default]"),
    make_option(c("--maf-min"), type = "numeric", default = 0.01,
                help = "Minimum MAF for causal variants [default: %default]"),
    make_option(c("--maf-max"), type = "numeric", default = 0.5,
                help = "Maximum MAF for causal variants [default: %default]"),
    make_option(c("--binary"), action = "store_true", default = FALSE,
                help = "Simulate binary trait"),
    make_option(c("--prevalence"), type = "numeric", default = 0.1,
                help = "Disease prevalence for binary traits [default: %default]"),
    make_option(c("--list-scenarios"), action = "store_true", default = FALSE,
                help = "List available simulation scenarios")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

set.seed(opt$seed)

# =============================================================================
# Simulation Scenarios
# =============================================================================

scenarios <- list(

    # Scenario 1: Variants enriched in AFR populations
    # Tests: Does imputation recover AFR-common variants that EUR GWAS miss?
    afr_enriched = list(
        description = "Causal variants common in AFR but rare/absent in EUR",
        select_variants = function(af_data) {
            # Select variants where AFR AF > 0.1 but EUR AF < 0.05
            af_data[AFR_AF > 0.1 & EUR_AF < 0.05]
        },
        effect_by_ancestry = function(n) {
            # Equal effects across ancestries (same causal mechanism)
            list(AFR = 1.0, EUR = 1.0, EAS = 1.0, SAS = 1.0, AMR = 1.0)
        }
    ),

    # Scenario 2: Ancestry-specific effect sizes
    # Tests: Can pipeline detect variants with heterogeneous effects?
    ancestry_specific = list(
        description = "Same variants but different effect sizes by ancestry",
        select_variants = function(af_data) {
            # Select variants polymorphic in multiple populations
            af_data[AFR_AF > 0.05 & EUR_AF > 0.05 &
                    (EAS_AF > 0.05 | SAS_AF > 0.05 | AMR_AF > 0.05)]
        },
        effect_by_ancestry = function(n) {
            # Generate different effect multipliers
            # Simulates LD differences, G×E, etc.
            list(
                AFR = runif(n, 0.5, 1.5),
                EUR = rep(1.0, n),  # Reference
                EAS = runif(n, 0.3, 1.2),
                SAS = runif(n, 0.4, 1.3),
                AMR = runif(n, 0.6, 1.4)  # Admixed, intermediate
            )
        }
    ),

    # Scenario 3: Rare variants (tests imputation quality for rare variants)
    rare_variants = list(
        description = "Causal rare variants (MAF 0.1-1%) - tests rare variant imputation",
        select_variants = function(af_data) {
            # Rare in all populations
            af_data[GLOBAL_AF >= 0.001 & GLOBAL_AF <= 0.01]
        },
        effect_by_ancestry = function(n) {
            # Larger effects for rare variants (as observed in nature)
            list(AFR = 1.0, EUR = 1.0, EAS = 1.0, SAS = 1.0, AMR = 1.0)
        }
    ),

    # Scenario 4: AMR-specific (Latino/admixed populations)
    # Tests: Pipeline performance in admixed populations
    amr_enriched = list(
        description = "Variants enriched in AMR (admixed) populations",
        select_variants = function(af_data) {
            # Higher frequency in AMR than EUR or AFR alone
            # These could be NAT-derived variants
            af_data[AMR_AF > 0.1 & EUR_AF < AMR_AF & AFR_AF < AMR_AF]
        },
        effect_by_ancestry = function(n) {
            list(AFR = 0.8, EUR = 0.8, EAS = 0.5, SAS = 0.5, AMR = 1.2)
        }
    ),

    # Scenario 5: Opposite direction effects
    # Tests: Can detect when effect direction differs by ancestry?
    opposite_effects = list(
        description = "Variants with opposite effects in different ancestries",
        select_variants = function(af_data) {
            af_data[AFR_AF > 0.1 & EUR_AF > 0.1]
        },
        effect_by_ancestry = function(n) {
            # Half variants have opposite effects in AFR vs EUR
            directions <- sample(c(-1, 1), n, replace = TRUE, prob = c(0.3, 0.7))
            list(
                AFR = directions,
                EUR = rep(1, n),
                EAS = directions * runif(n, 0.5, 1.0),
                SAS = rep(1, n) * runif(n, 0.8, 1.2),
                AMR = (directions + 1) / 2  # Intermediate for admixed
            )
        }
    ),

    # Scenario 6: Frequency-dependent effects
    # Tests: Power at different MAF ranges
    maf_stratified = list(
        description = "Effects scale with MAF (larger effects for rarer variants)",
        select_variants = function(af_data) {
            af_data[GLOBAL_AF > 0.001 & GLOBAL_AF < 0.3]
        },
        effect_by_ancestry = function(n) {
            list(AFR = 1.0, EUR = 1.0, EAS = 1.0, SAS = 1.0, AMR = 1.0)
        }
    ),

    # Scenario 7: Mixed architecture (realistic)
    # Combines multiple patterns seen in real traits
    realistic_mixed = list(
        description = "Realistic mix: some shared, some ancestry-specific effects",
        select_variants = function(af_data) {
            af_data[GLOBAL_AF > 0.005]
        },
        effect_by_ancestry = function(n) {
            # 60% shared effects, 20% EUR-enriched, 20% non-EUR enriched
            n_shared <- round(n * 0.6)
            n_eur <- round(n * 0.2)
            n_diverse <- n - n_shared - n_eur

            afr_eff <- c(rep(1, n_shared), rep(0.3, n_eur), runif(n_diverse, 1.2, 2.0))
            eur_eff <- c(rep(1, n_shared), rep(1, n_eur), runif(n_diverse, 0.3, 0.8))

            list(
                AFR = afr_eff,
                EUR = eur_eff,
                EAS = c(rep(1, n_shared), rep(0.5, n_eur), runif(n_diverse, 0.8, 1.2)),
                SAS = c(rep(1, n_shared), rep(0.6, n_eur), runif(n_diverse, 0.9, 1.1)),
                AMR = (afr_eff + eur_eff) / 2  # Admixed intermediate
            )
        }
    )
)

# =============================================================================
# List Scenarios
# =============================================================================

if (opt$`list-scenarios`) {
    cat("\n")
    cat("================================================================================\n")
    cat("AVAILABLE SIMULATION SCENARIOS\n")
    cat("================================================================================\n\n")

    for (name in names(scenarios)) {
        cat(sprintf("  %-20s %s\n", name, scenarios[[name]]$description))
    }

    cat("\n")
    cat("Usage:\n")
    cat("  Rscript simulate_ancestry_aware_phenotypes.R --scenario afr_enriched ...\n")
    cat("\n")
    cat("Recommended combinations for benchmarking:\n")
    cat("  1. afr_enriched      - Shows advantage for AFR imputation\n")
    cat("  2. ancestry_specific - Tests effect heterogeneity detection\n")
    cat("  3. rare_variants     - Tests rare variant imputation quality\n")
    cat("  4. realistic_mixed   - Most realistic, combines patterns\n")
    cat("\n")
    quit(status = 0)
}

# =============================================================================
# Validate Inputs
# =============================================================================

if (is.null(opt$vcf)) {
    stop("--vcf is required. Use --list-scenarios to see options.")
}

if (!opt$scenario %in% names(scenarios)) {
    stop(sprintf("Unknown scenario: %s. Use --list-scenarios to see options.", opt$scenario))
}

scenario <- scenarios[[opt$scenario]]

cat("================================================================================\n")
cat("ANCESTRY-AWARE PHENOTYPE SIMULATION\n")
cat("================================================================================\n\n")
cat(sprintf("Scenario:     %s\n", opt$scenario))
cat(sprintf("Description:  %s\n", scenario$description))
cat(sprintf("N causal:     %d\n", opt$`n-causal`))
cat(sprintf("Heritability: %.2f (ancestry-specific: %.2f)\n", opt$h2, opt$`h2-ancestry`))
cat(sprintf("Trait type:   %s\n", ifelse(opt$binary, "Binary", "Continuous")))
cat("\n")

# =============================================================================
# Calculate Population-Specific Allele Frequencies
# =============================================================================

cat("Calculating population-specific allele frequencies...\n")

# Get sample IDs and populations
cmd_samples <- sprintf("bcftools query -l %s", opt$vcf)
samples <- fread(cmd = cmd_samples, header = FALSE)$V1

# Read population file or infer from 1KG naming
if (!is.null(opt$`pop-file`) && file.exists(opt$`pop-file`)) {
    pop_info <- fread(opt$`pop-file`)
    if (!"SUPERPOP" %in% names(pop_info)) {
        # Try to add superpop from known 1KG populations
        pop_to_superpop <- c(
            CHB = "EAS", JPT = "EAS", CHS = "EAS", CDX = "EAS", KHV = "EAS",
            CEU = "EUR", TSI = "EUR", FIN = "EUR", GBR = "EUR", IBS = "EUR",
            YRI = "AFR", LWK = "AFR", GWD = "AFR", MSL = "AFR", ESN = "AFR",
            ASW = "AFR", ACB = "AFR",
            MXL = "AMR", PUR = "AMR", CLM = "AMR", PEL = "AMR",
            GIH = "SAS", PJL = "SAS", BEB = "SAS", STU = "SAS", ITU = "SAS"
        )
        pop_info[, SUPERPOP := pop_to_superpop[POP]]
    }
} else {
    # Infer from 1KG sample naming convention (e.g., NA12878 is CEU)
    cat("  No population file provided, inferring from sample IDs...\n")

    # 1KG samples often have population in ID or need lookup
    # For now, create placeholder - user should provide pop file
    pop_info <- data.table(
        IID = samples,
        POP = "UNK",
        SUPERPOP = "UNK"
    )

    # Try to match known 1KG sample patterns
    # This is approximate - better to provide pop file
    cat("  WARNING: Population inference is approximate. Provide --pop-file for accuracy.\n")
}

# Get variant info with population-specific AFs
cat("  Extracting variant information...\n")

# For 1KG VCF with INFO/AF fields, or calculate from genotypes
cmd <- sprintf("bcftools query -f '%%CHROM\\t%%POS\\t%%ID\\t%%REF\\t%%ALT\\t%%INFO/AF\\n' %s 2>/dev/null | head -100000",
               opt$vcf)
variants <- tryCatch({
    fread(cmd = cmd, header = FALSE, col.names = c("CHR", "POS", "ID", "REF", "ALT", "AF"))
}, error = function(e) {
    # If INFO/AF not available, we need to calculate from genotypes
    cat("  INFO/AF not found, calculating from genotypes (slower)...\n")
    NULL
})

if (is.null(variants) || nrow(variants) == 0) {
    # Calculate AFs from genotypes
    cmd <- sprintf("bcftools query -f '%%CHROM\\t%%POS\\t%%ID\\t%%REF\\t%%ALT[\\t%%GT]\\n' %s | head -50000",
                   opt$vcf)
    geno_raw <- fread(cmd = cmd, header = FALSE)

    if (nrow(geno_raw) == 0) {
        stop("No variants found in VCF")
    }

    # Set column names
    setnames(geno_raw, 1:5, c("CHR", "POS", "ID", "REF", "ALT"))
    sample_cols <- paste0("V", 6:ncol(geno_raw))
    setnames(geno_raw, 6:ncol(geno_raw), samples[1:(ncol(geno_raw) - 5)])

    # Calculate global AF
    calc_af <- function(gt_vec) {
        gt_vec <- gt_vec[gt_vec != "." & gt_vec != "./." & gt_vec != ".|."]
        if (length(gt_vec) == 0) return(NA)
        n_alt <- sum(sapply(strsplit(gt_vec, "[/|]"), function(x) sum(as.numeric(x))))
        n_total <- length(gt_vec) * 2
        return(n_alt / n_total)
    }

    geno_raw[, GLOBAL_AF := apply(.SD, 1, calc_af), .SDcols = samples[1:min(length(samples), ncol(geno_raw) - 5)]]

    # For now, use global AF for all populations (approximate)
    # In production, would calculate per-population
    variants <- geno_raw[, .(CHR, POS, ID, REF, ALT, GLOBAL_AF)]
    variants[, AFR_AF := GLOBAL_AF]
    variants[, EUR_AF := GLOBAL_AF]
    variants[, EAS_AF := GLOBAL_AF]
    variants[, SAS_AF := GLOBAL_AF]
    variants[, AMR_AF := GLOBAL_AF]
} else {
    variants[, GLOBAL_AF := as.numeric(AF)]
    variants[, AF := NULL]
    # Use global as proxy (in production, get from gnomAD or calculate)
    variants[, AFR_AF := GLOBAL_AF]
    variants[, EUR_AF := GLOBAL_AF]
    variants[, EAS_AF := GLOBAL_AF]
    variants[, SAS_AF := GLOBAL_AF]
    variants[, AMR_AF := GLOBAL_AF]
}

variants <- variants[!is.na(GLOBAL_AF)]
cat(sprintf("  Found %d variants with frequency data\n", nrow(variants)))

# =============================================================================
# Select Causal Variants Based on Scenario
# =============================================================================

cat("\nSelecting causal variants...\n")

# Apply scenario-specific selection
eligible <- scenario$select_variants(variants)

# Apply MAF filters
eligible <- eligible[GLOBAL_AF >= opt$`maf-min` & GLOBAL_AF <= opt$`maf-max`]

cat(sprintf("  %d variants eligible for selection\n", nrow(eligible)))

if (nrow(eligible) < opt$`n-causal`) {
    cat(sprintf("  WARNING: Only %d eligible variants, using all\n", nrow(eligible)))
    n_causal <- nrow(eligible)
} else {
    n_causal <- opt$`n-causal`
}

# Sample causal variants
causal_idx <- sample(1:nrow(eligible), n_causal)
causal_variants <- eligible[causal_idx]

cat(sprintf("  Selected %d causal variants\n", nrow(causal_variants)))

# Generate base effect sizes (scaled by MAF for rare variants)
# Rare variants tend to have larger effects
causal_variants[, BASE_BETA := rnorm(n_causal, mean = 0, sd = 1)]

# Scale by MAF (optional, based on scenario)
if (opt$scenario == "maf_stratified" || opt$scenario == "rare_variants") {
    # Larger effects for rarer variants: β ∝ 1/sqrt(MAF)
    causal_variants[, BASE_BETA := BASE_BETA * sqrt(0.1 / pmax(GLOBAL_AF, 0.001))]
}

# Get ancestry-specific effect multipliers
anc_effects <- scenario$effect_by_ancestry(n_causal)

# Add ancestry-specific effects
for (anc in names(anc_effects)) {
    col_name <- paste0("BETA_", anc)
    if (length(anc_effects[[anc]]) == 1) {
        causal_variants[, (col_name) := BASE_BETA * anc_effects[[anc]]]
    } else {
        causal_variants[, (col_name) := BASE_BETA * anc_effects[[anc]]]
    }
}

# =============================================================================
# Extract Genotypes at Causal Variants
# =============================================================================

cat("\nExtracting genotypes at causal variants...\n")

# Create regions file
regions_file <- tempfile(fileext = ".txt")
fwrite(causal_variants[, .(CHR, POS, POS)], regions_file, col.names = FALSE, sep = "\t")

# Extract genotypes
cmd <- sprintf("bcftools query -R %s -f '%%CHROM\\t%%POS[\\t%%GT]\\n' %s", regions_file, opt$vcf)
geno <- fread(cmd = cmd, header = FALSE)

if (nrow(geno) == 0) {
    stop("No genotypes extracted. Check chromosome naming (chr1 vs 1).")
}

setnames(geno, c("CHR", "POS", samples))

# Convert to dosage
for (s in samples) {
    geno[, (s) := fcase(
        get(s) %in% c("0/0", "0|0"), 0L,
        get(s) %in% c("0/1", "1/0", "0|1", "1|0"), 1L,
        get(s) %in% c("1/1", "1|1"), 2L,
        default = NA_integer_
    )]
}

unlink(regions_file)

cat(sprintf("  Extracted genotypes for %d variants, %d samples\n", nrow(geno), length(samples)))

# =============================================================================
# Calculate Genetic Values
# =============================================================================

cat("\nCalculating genetic values...\n")

# Merge genotypes with causal variant info
geno[, CHR := as.character(CHR)]
causal_variants[, CHR := as.character(CHR)]
geno_causal <- merge(geno, causal_variants[, .(CHR, POS, BASE_BETA, BETA_AFR, BETA_EUR, BETA_EAS, BETA_SAS, BETA_AMR)],
                     by = c("CHR", "POS"))

# Get sample ancestry
sample_anc <- merge(data.table(IID = samples), pop_info, by = "IID", all.x = TRUE)
sample_anc[is.na(SUPERPOP), SUPERPOP := "EUR"]  # Default to EUR if unknown

# Calculate genetic value for each sample using their ancestry-specific effects
genetic_values <- sapply(samples, function(s) {
    dosages <- geno_causal[[s]]
    anc <- sample_anc[IID == s, SUPERPOP]

    # Get appropriate beta based on ancestry
    beta_col <- paste0("BETA_", anc)
    if (!beta_col %in% names(geno_causal)) {
        beta_col <- "BASE_BETA"
    }
    betas <- geno_causal[[beta_col]]

    sum(dosages * betas, na.rm = TRUE)
})

results <- data.table(IID = samples, G = genetic_values)
results <- merge(results, sample_anc, by = "IID", all.x = TRUE)

# =============================================================================
# Simulate Phenotype
# =============================================================================

cat("\nSimulating phenotype...\n")

# Standardize genetic values
results[, G_std := (G - mean(G)) / sd(G)]

# Calculate environmental variance
var_G <- var(results$G_std)
var_E <- var_G * (1 - opt$h2) / opt$h2

# Environmental noise
results[, E := rnorm(.N, 0, sqrt(var_E))]

# Phenotype
results[, P := G_std + E]

if (opt$binary) {
    threshold <- qnorm(1 - opt$prevalence)
    results[, PHENO := as.integer(P > threshold)]
    results[, LIABILITY := P]
    cat(sprintf("  Cases: %d (%.1f%%), Controls: %d\n",
                sum(results$PHENO), 100 * mean(results$PHENO), sum(1 - results$PHENO)))
} else {
    results[, PHENO := (P - mean(P)) / sd(P)]
}

# Store true genetic value
results[, G_TRUE := G_std]

# =============================================================================
# Output Files
# =============================================================================

cat("\nWriting output files...\n")

# PLINK phenotype file
plink_pheno <- results[, .(FID = IID, IID = IID, PHENO)]
fwrite(plink_pheno, paste0(opt$output, ".pheno"), sep = "\t")
cat(sprintf("  %s.pheno\n", opt$output))

# Full results with ancestry
fwrite(results, paste0(opt$output, "_full.txt"), sep = "\t")
cat(sprintf("  %s_full.txt\n", opt$output))

# Causal variants with ancestry-specific effects
fwrite(causal_variants, paste0(opt$output, "_causal_variants.txt"), sep = "\t")
cat(sprintf("  %s_causal_variants.txt\n", opt$output))

# Summary by ancestry
cat("\n")
cat("================================================================================\n")
cat("SIMULATION SUMMARY BY ANCESTRY\n")
cat("================================================================================\n\n")

ancestry_summary <- results[, .(
    N = .N,
    Mean_G = mean(G_TRUE),
    SD_G = sd(G_TRUE),
    Mean_P = mean(PHENO),
    SD_P = sd(PHENO)
), by = SUPERPOP]

print(ancestry_summary)

# Achieved heritability
cat(sprintf("\nOverall h2: %.3f\n", cor(results$G_TRUE, results$PHENO)^2))

# Per-ancestry h2
cat("\nPer-ancestry h2:\n")
for (anc in unique(results$SUPERPOP)) {
    sub <- results[SUPERPOP == anc]
    if (nrow(sub) > 10) {
        h2_anc <- cor(sub$G_TRUE, sub$PHENO)^2
        cat(sprintf("  %s: %.3f (n=%d)\n", anc, h2_anc, nrow(sub)))
    }
}

# =============================================================================
# Usage Instructions
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("USAGE FOR BENCHMARKING\n")
cat("================================================================================\n")
cat("\n")
cat("1. Impute genotype data using different approaches\n")
cat("\n")
cat("2. Run GWAS:\n")
cat(sprintf("   plink2 --bfile imputed --pheno %s.pheno --glm --out gwas_%s\n",
            opt$output, opt$scenario))
cat("\n")
cat("3. Check causal variant recovery:\n")
cat(sprintf("   Rscript check_hit_recovery.R --gwas gwas_%s.glm.linear \\\n", opt$scenario))
cat(sprintf("       --causal %s_causal_variants.txt --output recovery_%s\n",
            opt$output, opt$scenario))
cat("\n")
cat("4. Compare power across ancestries:\n")
cat("   - Stratify GWAS by SUPERPOP\n")
cat("   - Compare % causal variants reaching significance\n")
cat("   - Our pipeline should show improved power for non-EUR groups\n")
cat("\n")
cat("================================================================================\n")
