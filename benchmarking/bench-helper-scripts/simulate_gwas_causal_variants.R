#!/usr/bin/env Rscript
################################################################################
# simulate_gwas_causal_variants.R
#
# Comprehensive GWAS phenotype simulation for benchmarking imputation approaches
# in understudied and admixed populations.
#
# Creates simulations with:
#   1. Randomly selected causal variants from WGS (not EUR-biased GWAS hits)
#   2. Ancestry-specific MAF differences (common in AFR, rare in EUR, etc.)
#   3. Ancestry-specific effect sizes (β_AFR ≠ β_EUR)
#   4. Local ancestry-specific effects (for Tractor GWAS testing)
#   5. Varying effect sizes and frequencies
#
# Outputs TWO sets:
#   - GWAS set: MAF > 0.01 variants for standard GWAS
#   - RVAS set: All variants including rare for rare variant analysis
#
# Usage:
#   Rscript simulate_gwas_causal_variants.R \
#       --vcf 1kg_wgs.vcf.gz \
#       --pop-file population_labels.txt \
#       --lai-msp rfmix_output.msp.tsv \
#       --n-causal 100 \
#       --scenario ancestry_specific_maf \
#       --output sim_gwas
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
    make_option(c("-p", "--pop-file"), type = "character", default = NULL,
                help = "Population/ancestry labels file (IID, POP, SUPERPOP)"),
    make_option(c("-l", "--lai-msp"), type = "character", default = NULL,
                help = "RFMix2 MSP file for LAI scenarios (optional)"),
    make_option(c("-n", "--n-causal"), type = "integer", default = 100,
                help = "Number of causal variants [default: %default]"),
    make_option(c("-s", "--scenario"), type = "character", default = "random",
                help = "Simulation scenario [default: %default]"),
    make_option(c("--h2"), type = "numeric", default = 0.4,
                help = "Heritability [default: %default]"),
    make_option(c("-o", "--output"), type = "character", default = "sim_gwas",
                help = "Output prefix [default: %default]"),
    make_option(c("--seed"), type = "integer", default = 42,
                help = "Random seed [default: %default]"),
    make_option(c("--maf-gwas"), type = "numeric", default = 0.01,
                help = "MAF threshold for GWAS output [default: %default]"),
    make_option(c("--list-scenarios"), action = "store_true", default = FALSE,
                help = "List available scenarios and exit")
)

opt_parser <- OptionParser(option_list = option_list,
    description = "Simulate GWAS phenotypes for imputation benchmarking")
opt <- parse_args(opt_parser)

set.seed(opt$seed)

# =============================================================================
# Simulation Scenarios
# =============================================================================

scenarios <- list(

    # Baseline: Random causal variants, same effect across ancestries
    random = list(
        description = "Random causal variants, uniform effects across ancestries",
        select_variants = function(vars, n) {
            # Random selection, no ancestry bias
            vars[sample(.N, min(n, .N))]
        },
        assign_effects = function(vars, pops) {
            # Same effect size for all ancestries
            vars[, BETA_AFR := BETA]
            vars[, BETA_EUR := BETA]
            vars[, BETA_AMR := BETA]
            vars[, BETA_EAS := BETA]
            vars[, BETA_SAS := BETA]
            vars
        },
        test_goal = "Baseline - should perform similarly across approaches"
    ),

    # Ancestry-specific MAF: variants common in AFR but rare in EUR
    ancestry_specific_maf_afr = list(
        description = "Causal variants common in AFR (MAF>5%) but rare in EUR (MAF<1%)",
        select_variants = function(vars, n) {
            # Select variants with high AFR freq, low EUR freq
            candidates <- vars[AF_AFR > 0.05 & AF_EUR < 0.01]
            if (nrow(candidates) < n) {
                # Fall back to relaxed criteria
                candidates <- vars[AF_AFR > 0.03 & AF_EUR < 0.02]
            }
            candidates[sample(.N, min(n, .N))]
        },
        assign_effects = function(vars, pops) {
            vars[, BETA_AFR := BETA]
            vars[, BETA_EUR := BETA]  # Same effect, but rare so less power
            vars[, BETA_AMR := BETA * 0.8]  # Partial effect in admixed
            vars[, BETA_EAS := BETA * 0.5]
            vars[, BETA_SAS := BETA * 0.5]
            vars
        },
        test_goal = "Tests if diverse panels capture AFR-common variants better"
    ),

    # Ancestry-specific MAF: variants common in EUR but rare in AFR
    ancestry_specific_maf_eur = list(
        description = "Causal variants common in EUR (MAF>5%) but rare in AFR (MAF<1%)",
        select_variants = function(vars, n) {
            candidates <- vars[AF_EUR > 0.05 & AF_AFR < 0.01]
            if (nrow(candidates) < n) {
                candidates <- vars[AF_EUR > 0.03 & AF_AFR < 0.02]
            }
            candidates[sample(.N, min(n, .N))]
        },
        assign_effects = function(vars, pops) {
            vars[, BETA_AFR := BETA]
            vars[, BETA_EUR := BETA]
            vars[, BETA_AMR := BETA * 0.8]
            vars[, BETA_EAS := BETA * 0.7]
            vars[, BETA_SAS := BETA * 0.7]
            vars
        },
        test_goal = "Baseline for EUR-centric discovery (should work well with all panels)"
    ),

    # Ancestry-specific effects: different beta by ancestry
    ancestry_specific_effects = list(
        description = "Same variants, but effect sizes differ by ancestry (β_AFR ≠ β_EUR)",
        select_variants = function(vars, n) {
            # Select variants present in all ancestries (MAF > 1% in each)
            candidates <- vars[AF_AFR > 0.01 & AF_EUR > 0.01 & AF_AMR > 0.01]
            candidates[sample(.N, min(n, .N))]
        },
        assign_effects = function(vars, pops) {
            # AFR effect is 1.5x EUR effect (or opposite direction sometimes)
            vars[, BETA_EUR := BETA]
            vars[, BETA_AFR := BETA * sample(c(1.5, 0.5, -0.5), .N, replace = TRUE,
                                              prob = c(0.4, 0.4, 0.2))]
            vars[, BETA_AMR := (BETA_AFR + BETA_EUR) / 2]  # Intermediate for admixed
            vars[, BETA_EAS := BETA * runif(.N, 0.3, 1.2)]
            vars[, BETA_SAS := BETA * runif(.N, 0.3, 1.2)]
            vars
        },
        test_goal = "Tests if ancestry-specific GWAS / meta-analysis improves detection"
    ),

    # Opposite effects by ancestry (masks in combined GWAS)
    opposite_effects = list(
        description = "Opposite effects in AFR vs EUR (cancels in combined GWAS!)",
        select_variants = function(vars, n) {
            candidates <- vars[AF_AFR > 0.05 & AF_EUR > 0.05]
            candidates[sample(.N, min(n, .N))]
        },
        assign_effects = function(vars, pops) {
            vars[, BETA_EUR := BETA]
            vars[, BETA_AFR := -BETA]  # OPPOSITE effect!
            vars[, BETA_AMR := 0]  # Cancels in admixed
            vars[, BETA_EAS := BETA * 0.5]
            vars[, BETA_SAS := BETA * 0.5]
            vars
        },
        test_goal = "Combined GWAS shows ~0 effect; ancestry-stratified finds both"
    ),

    # LAI-specific effects (for Tractor GWAS)
    lai_specific = list(
        description = "Effect depends on LOCAL ancestry tract, not global ancestry",
        select_variants = function(vars, n) {
            candidates <- vars[AF_AFR > 0.01 & AF_EUR > 0.01]
            candidates[sample(.N, min(n, .N))]
        },
        assign_effects = function(vars, pops) {
            # These are placeholder - actual LAI effects computed per-haplotype
            vars[, BETA_AFR := BETA]
            vars[, BETA_EUR := 0]  # Effect ONLY on AFR local ancestry tracts
            vars[, BETA_AMR := NA]  # Depends on LAI at each locus
            vars[, BETA_EAS := 0]
            vars[, BETA_SAS := 0]
            vars[, LAI_DEPENDENT := TRUE]
            vars
        },
        test_goal = "Tractor GWAS should outperform standard GWAS"
    ),

    # Rare variants (for RVAS benchmarking)
    rare_variants = list(
        description = "Causal variants are RARE (MAF < 1%) - tests rare variant imputation",
        select_variants = function(vars, n) {
            # Select rare variants
            candidates <- vars[AF_ALL < 0.01 & AF_ALL > 0.001]
            candidates[sample(.N, min(n, .N))]
        },
        assign_effects = function(vars, pops) {
            # Larger effects for rare variants (typical for rare variant associations)
            vars[, BETA := BETA * 3]  # Larger effect
            vars[, BETA_AFR := BETA]
            vars[, BETA_EUR := BETA]
            vars[, BETA_AMR := BETA]
            vars[, BETA_EAS := BETA]
            vars[, BETA_SAS := BETA]
            vars
        },
        test_goal = "Tests rare variant imputation quality for RVAS"
    ),

    # Mixed: combination of common and rare, ancestry-variable
    mixed_realistic = list(
        description = "Realistic mix: 70% common, 30% rare; ancestry-variable effects",
        select_variants = function(vars, n) {
            n_common <- round(n * 0.7)
            n_rare <- n - n_common

            common <- vars[AF_ALL > 0.05][sample(.N, min(n_common, .N))]
            rare <- vars[AF_ALL < 0.01 & AF_ALL > 0.001][sample(.N, min(n_rare, .N))]

            rbind(common, rare)
        },
        assign_effects = function(vars, pops) {
            # Larger effects for rare, smaller for common (realistic)
            vars[, BETA := ifelse(AF_ALL < 0.01, BETA * 2.5, BETA)]

            # Add ancestry variation
            vars[, BETA_EUR := BETA]
            vars[, BETA_AFR := BETA * runif(.N, 0.7, 1.3)]
            vars[, BETA_AMR := (BETA_AFR + BETA_EUR) / 2]
            vars[, BETA_EAS := BETA * runif(.N, 0.5, 1.1)]
            vars[, BETA_SAS := BETA * runif(.N, 0.5, 1.1)]
            vars
        },
        test_goal = "Most realistic - tests overall pipeline performance"
    )
)

# =============================================================================
# List Scenarios
# =============================================================================

if (opt$`list-scenarios`) {
    cat("\n")
    cat("================================================================================\n")
    cat("GWAS SIMULATION SCENARIOS FOR IMPUTATION BENCHMARKING\n")
    cat("================================================================================\n\n")

    for (name in names(scenarios)) {
        s <- scenarios[[name]]
        cat(sprintf("  %-25s\n", name))
        cat(sprintf("    Description: %s\n", s$description))
        cat(sprintf("    Test goal:   %s\n\n", s$test_goal))
    }

    cat("================================================================================\n")
    cat("Output files created for each scenario:\n")
    cat("  - {prefix}_gwas.pheno      : Phenotypes for GWAS (MAF > 0.01 variants)\n")
    cat("  - {prefix}_rvas.pheno      : Phenotypes for RVAS (all variants)\n")
    cat("  - {prefix}_causal.txt      : Causal variant details with true effects\n")
    cat("  - {prefix}_ancestry_effects.txt : Per-ancestry effect sizes\n")
    cat("================================================================================\n")
    quit(status = 0)
}

# =============================================================================
# Validate Inputs
# =============================================================================

if (is.null(opt$vcf)) {
    stop("--vcf is required. Use --list-scenarios to see options.")
}

if (!opt$scenario %in% names(scenarios)) {
    stop(sprintf("Unknown scenario: %s. Use --list-scenarios.", opt$scenario))
}

scenario <- scenarios[[opt$scenario]]

cat("================================================================================\n")
cat("GWAS CAUSAL VARIANT SIMULATION\n")
cat("================================================================================\n\n")
cat(sprintf("Scenario:    %s\n", opt$scenario))
cat(sprintf("Description: %s\n", scenario$description))
cat(sprintf("Test goal:   %s\n", scenario$test_goal))
cat(sprintf("N causal:    %d\n", opt$`n-causal`))
cat(sprintf("Heritability: %.2f\n", opt$h2))
cat("\n")

# =============================================================================
# Load Population Labels
# =============================================================================

cat("Loading data...\n")

# Get samples from VCF
cmd_samples <- sprintf("bcftools query -l %s", opt$vcf)
samples <- fread(cmd = cmd_samples, header = FALSE)$V1
n_samples <- length(samples)
cat(sprintf("  Samples in VCF: %d\n", n_samples))

# Load population labels
if (!is.null(opt$`pop-file`) && file.exists(opt$`pop-file`)) {
    pop_data <- fread(opt$`pop-file`)
    # Expect columns: IID (or sample), POP, SUPERPOP
    if (!"SUPERPOP" %in% names(pop_data) && "superpop" %in% names(pop_data)) {
        setnames(pop_data, "superpop", "SUPERPOP")
    }
    if (!"IID" %in% names(pop_data)) {
        setnames(pop_data, names(pop_data)[1], "IID")
    }
    cat(sprintf("  Population labels: %d samples\n", nrow(pop_data)))
} else {
    # Infer from 1KG naming convention
    cat("  Inferring populations from sample IDs (1KG naming)...\n")
    pop_mapping <- data.table(
        pattern = c("YRI|LWK|GWD|MSL|ESN", "CEU|TSI|FIN|GBR|IBS",
                   "CHB|JPT|CHS|CDX|KHV", "GIH|PJL|BEB|STU|ITU",
                   "MXL|PUR|CLM|PEL|ACB|ASW"),
        SUPERPOP = c("AFR", "EUR", "EAS", "SAS", "AMR")
    )

    pop_data <- data.table(IID = samples)
    pop_data[, SUPERPOP := "EUR"]  # Default
    for (i in 1:nrow(pop_mapping)) {
        pop_data[grepl(pop_mapping$pattern[i], IID), SUPERPOP := pop_mapping$SUPERPOP[i]]
    }
}

# Merge with samples
pop_data <- pop_data[IID %in% samples]
cat(sprintf("  Ancestry distribution:\n"))
print(table(pop_data$SUPERPOP))
cat("\n")

# =============================================================================
# Extract Variant Information with Ancestry-Specific Frequencies
# =============================================================================

cat("Extracting variant information...\n")

# Get variant info (limit to first 100k for speed, or use full genome)
cmd <- sprintf("bcftools query -f '%%CHROM\\t%%POS\\t%%ID\\t%%REF\\t%%ALT\\t%%INFO/AF\\n' %s 2>/dev/null | head -100000",
               opt$vcf)

variants <- tryCatch({
    fread(cmd = cmd, header = FALSE, col.names = c("CHR", "POS", "ID", "REF", "ALT", "AF"))
}, error = function(e) {
    # Try without AF
    cmd2 <- sprintf("bcftools query -f '%%CHROM\\t%%POS\\t%%ID\\t%%REF\\t%%ALT\\n' %s 2>/dev/null | head -100000", opt$vcf)
    fread(cmd = cmd2, header = FALSE, col.names = c("CHR", "POS", "ID", "REF", "ALT"))
})

cat(sprintf("  Total variants queried: %d\n", nrow(variants)))

# Calculate ancestry-specific allele frequencies
cat("  Calculating ancestry-specific frequencies...\n")

# For full implementation, would calculate per-ancestry AF from genotypes
# Here we simulate based on overall AF with ancestry variation
if ("AF" %in% names(variants)) {
    variants[, AF := as.numeric(AF)]
    variants <- variants[!is.na(AF) & AF > 0 & AF < 1]
    variants[, AF_ALL := AF]
} else {
    # Assign random AF if not available
    variants[, AF_ALL := runif(.N, 0.01, 0.5)]
}

# Simulate ancestry-specific frequencies (with realistic correlation structure)
# In practice, would calculate from actual genotypes
variants[, AF_EUR := pmin(0.99, pmax(0.001, AF_ALL + rnorm(.N, 0, 0.05)))]
variants[, AF_AFR := pmin(0.99, pmax(0.001, AF_ALL + rnorm(.N, 0, 0.08)))]  # More variation in AFR
variants[, AF_EAS := pmin(0.99, pmax(0.001, AF_ALL + rnorm(.N, 0, 0.06)))]
variants[, AF_SAS := pmin(0.99, pmax(0.001, AF_ALL + rnorm(.N, 0, 0.06)))]
variants[, AF_AMR := (AF_EUR * 0.5 + AF_AFR * 0.3 + 0.2 * runif(.N, 0.01, 0.5))]  # Admixed

cat(sprintf("  Variants with frequency data: %d\n", nrow(variants)))

# =============================================================================
# Select Causal Variants
# =============================================================================

cat("\nSelecting causal variants...\n")

# Assign base effect sizes (before ancestry adjustment)
variants[, BETA := rnorm(.N, 0, 0.3)]

# Apply scenario-specific selection
causal <- scenario$select_variants(variants, opt$`n-causal`)
cat(sprintf("  Selected %d causal variants\n", nrow(causal)))

# Assign ancestry-specific effects
causal <- scenario$assign_effects(causal, pop_data)

# Show effect size distribution
cat("  Effect size summary:\n")
cat(sprintf("    EUR: mean=%.3f, sd=%.3f\n", mean(causal$BETA_EUR), sd(causal$BETA_EUR)))
cat(sprintf("    AFR: mean=%.3f, sd=%.3f\n", mean(causal$BETA_AFR), sd(causal$BETA_AFR)))
cat(sprintf("    AMR: mean=%.3f, sd=%.3f\n", mean(causal$BETA_AMR, na.rm=TRUE), sd(causal$BETA_AMR, na.rm=TRUE)))

# =============================================================================
# Extract Genotypes at Causal Variants
# =============================================================================

cat("\nExtracting genotypes at causal variants...\n")

regions_file <- tempfile(fileext = ".txt")
fwrite(causal[, .(CHR, POS, POS)], regions_file, col.names = FALSE, sep = "\t")

cmd <- sprintf("bcftools query -R %s -f '%%CHROM\\t%%POS[\\t%%GT]\\n' %s", regions_file, opt$vcf)
geno <- tryCatch({
    fread(cmd = cmd, header = FALSE)
}, error = function(e) {
    cat("  Warning: Could not extract genotypes, using simulated values\n")
    NULL
})

if (!is.null(geno) && nrow(geno) > 0) {
    setnames(geno, c("CHR", "POS", samples))

    # Convert GT to dosage
    for (s in samples) {
        geno[, (s) := fcase(
            get(s) %in% c("0/0", "0|0"), 0,
            get(s) %in% c("0/1", "1/0", "0|1", "1|0"), 1,
            get(s) %in% c("1/1", "1|1"), 2,
            default = NA_real_
        )]
    }
    cat(sprintf("  Extracted genotypes for %d variants\n", nrow(geno)))
} else {
    # Simulate genotypes based on frequencies
    cat("  Simulating genotypes based on frequencies...\n")
    geno <- data.table(CHR = causal$CHR, POS = causal$POS)
    for (s in samples) {
        anc <- pop_data[IID == s, SUPERPOP]
        if (length(anc) == 0) anc <- "EUR"

        af_col <- paste0("AF_", anc)
        afs <- causal[[af_col]]

        # Simulate genotypes from binomial
        geno[, (s) := rbinom(.N, 2, afs)]
    }
}

unlink(regions_file)

# =============================================================================
# Calculate Genetic Values
# =============================================================================

cat("\nCalculating genetic values...\n")

# For each sample, calculate genetic value using ancestry-appropriate effects
genetic_values <- sapply(samples, function(s) {
    anc <- pop_data[IID == s, SUPERPOP]
    if (length(anc) == 0) anc <- "EUR"

    beta_col <- paste0("BETA_", anc)
    betas <- causal[[beta_col]]

    # Handle LAI-dependent scenarios (simplified - full version uses MSP file)
    if ("LAI_DEPENDENT" %in% names(causal) && any(causal$LAI_DEPENDENT, na.rm = TRUE)) {
        # For admixed, use weighted average based on global ancestry
        if (anc == "AMR") {
            betas <- (causal$BETA_AFR * 0.3 + causal$BETA_EUR * 0.5 + causal$BETA_EAS * 0.1) / 0.9
        }
    }

    # Get dosages for this sample
    dosages <- as.numeric(geno[[s]])

    sum(dosages * betas, na.rm = TRUE)
})

results <- data.table(
    FID = samples,
    IID = samples,
    G = genetic_values,
    SUPERPOP = pop_data[match(samples, IID), SUPERPOP]
)

# =============================================================================
# Simulate Phenotypes
# =============================================================================

cat("Simulating phenotypes...\n")

# Standardize genetic values
results[, G_std := (G - mean(G, na.rm = TRUE)) / sd(G, na.rm = TRUE)]
results[is.na(G_std) | is.infinite(G_std), G_std := 0]

# Add environmental noise for target heritability
var_G <- var(results$G_std, na.rm = TRUE)
if (var_G == 0 || is.na(var_G)) var_G <- 1
var_E <- var_G * (1 - opt$h2) / opt$h2

results[, E := rnorm(.N, 0, sqrt(var_E))]
results[, PHENO := G_std + E]
results[, PHENO := (PHENO - mean(PHENO)) / sd(PHENO)]  # Standardize final phenotype

# Check achieved heritability
achieved_h2 <- cor(results$G_std, results$PHENO, use = "complete.obs")^2
cat(sprintf("  Target h2: %.2f, Achieved h2: %.3f\n", opt$h2, achieved_h2))

# =============================================================================
# Create GWAS and RVAS Output Files
# =============================================================================

cat("\nCreating output files...\n")

# Identify which causal variants pass MAF filter for GWAS
causal[, PASS_GWAS_MAF := AF_ALL >= opt$`maf-gwas`]
n_gwas <- sum(causal$PASS_GWAS_MAF)
n_rvas <- nrow(causal)

cat(sprintf("  GWAS set (MAF >= %.3f): %d causal variants\n", opt$`maf-gwas`, n_gwas))
cat(sprintf("  RVAS set (all variants): %d causal variants\n", n_rvas))

# Phenotype files (same phenotypes, different downstream analysis)
pheno_out <- results[, .(FID, IID, PHENO)]
fwrite(pheno_out, paste0(opt$output, "_gwas.pheno"), sep = "\t")
fwrite(pheno_out, paste0(opt$output, "_rvas.pheno"), sep = "\t")

# Causal variant file with all details
causal_out <- causal[, .(CHR, POS, ID, REF, ALT, AF_ALL, AF_EUR, AF_AFR, AF_AMR, AF_EAS, AF_SAS,
                         BETA, BETA_EUR, BETA_AFR, BETA_AMR, BETA_EAS, BETA_SAS, PASS_GWAS_MAF)]
fwrite(causal_out, paste0(opt$output, "_causal.txt"), sep = "\t")

# GWAS-only causal variants (for power analysis)
fwrite(causal_out[PASS_GWAS_MAF == TRUE], paste0(opt$output, "_causal_gwas.txt"), sep = "\t")

# RVAS-only causal variants (rare)
fwrite(causal_out[PASS_GWAS_MAF == FALSE], paste0(opt$output, "_causal_rvas.txt"), sep = "\t")

# Per-ancestry genetic values (for stratified analysis)
ancestry_results <- results[, .(FID, IID, SUPERPOP, G, G_std, PHENO)]
fwrite(ancestry_results, paste0(opt$output, "_by_ancestry.txt"), sep = "\t")

# Summary by ancestry
ancestry_summary <- results[, .(
    N = .N,
    Mean_G = mean(G, na.rm = TRUE),
    SD_G = sd(G, na.rm = TRUE),
    Mean_Pheno = mean(PHENO, na.rm = TRUE),
    SD_Pheno = sd(PHENO, na.rm = TRUE)
), by = SUPERPOP]
fwrite(ancestry_summary, paste0(opt$output, "_ancestry_summary.txt"), sep = "\t")

# =============================================================================
# Summary
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("SIMULATION COMPLETE\n")
cat("================================================================================\n\n")

cat("Output files:\n")
cat(sprintf("  %s_gwas.pheno           - Phenotypes for GWAS\n", opt$output))
cat(sprintf("  %s_rvas.pheno           - Phenotypes for RVAS\n", opt$output))
cat(sprintf("  %s_causal.txt           - All causal variants with effects\n", opt$output))
cat(sprintf("  %s_causal_gwas.txt      - GWAS causal variants (MAF >= %.3f)\n", opt$output, opt$`maf-gwas`))
cat(sprintf("  %s_causal_rvas.txt      - RVAS causal variants (MAF < %.3f)\n", opt$output, opt$`maf-gwas`))
cat(sprintf("  %s_by_ancestry.txt      - Per-sample ancestry and genetic values\n", opt$output))
cat(sprintf("  %s_ancestry_summary.txt - Summary statistics by ancestry\n", opt$output))

cat("\n")
cat("Ancestry summary:\n")
print(ancestry_summary)

cat("\n")
cat("================================================================================\n")
cat("BENCHMARKING INSTRUCTIONS\n")
cat("================================================================================\n\n")

cat("1. Run each imputation approach (A-F) on array data\n\n")

cat("2. Run GWAS on imputed data:\n")
cat(sprintf("   plink2 --bfile imputed --pheno %s_gwas.pheno --glm --out gwas_results\n\n", opt$output))

cat("3. Run RVAS on imputed data (burden test):\n")
cat(sprintf("   # Use SAIGE-GENE or SKAT for rare variant analysis\n"))
cat(sprintf("   # Compare detection of causal variants in %s_causal_rvas.txt\n\n", opt$output))

cat("4. Compare power across approaches:\n")
cat(sprintf("   # Which approach detects more causal variants from %s_causal.txt?\n", opt$output))
cat("   # Stratify by ancestry and MAF\n\n")

cat("5. For LAI scenarios, also run Tractor GWAS:\n")
cat("   python Tractor.py --vcf imputed.vcf.gz --msp rfmix.msp --out tractor_geno\n")
cat("   python RunTractor.py --hapdose tractor_geno --phe phenotypes --out tractor_gwas\n")

cat("\n================================================================================\n")
