#!/usr/bin/env Rscript
################################################################################
# simulate_lai_phenotypes.R
#
# Simulates phenotypes where effects depend on LOCAL ANCESTRY context.
# This tests whether LAI-aware GWAS (Tractor) improves over standard GWAS.
#
# Key concept: The same variant can have different effects depending on
# the local ancestry background it's on. This happens due to:
#   - Different LD structure (causal variant vs tag SNP correlation)
#   - G×G interactions with ancestry-specific variants
#   - Local regulatory differences
#
# Scenarios tested:
#   1. Effects only on AFR background (variant only causal in AFR haplotypes)
#   2. Effects only on EUR background
#   3. Opposite effects by local ancestry (masking in standard GWAS)
#   4. LD-driven differences (tag SNP on array, causal nearby)
#
# For Tractor GWAS comparison:
#   - Standard GWAS may miss or attenuate effects
#   - Tractor should detect ancestry-specific effects
#   - Shows value of LAI module in pipeline
#
# Usage:
#   Rscript simulate_lai_phenotypes.R \
#       --vcf 1kg_wgs.vcf.gz \
#       --lai-msp rfmix_output.msp.tsv \
#       --scenario lai_afr_only \
#       --output lai_sim
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
                help = "VCF file with phased genotypes"),
    make_option(c("-l", "--lai-msp"), type = "character", default = NULL,
                help = "RFMix2 MSP file with local ancestry calls"),
    make_option(c("--scenario"), type = "character", default = "lai_afr_only",
                help = "Simulation scenario [default: %default]"),
    make_option(c("--n-causal"), type = "integer", default = 30,
                help = "Number of causal variants [default: %default]"),
    make_option(c("--h2"), type = "numeric", default = 0.4,
                help = "Heritability [default: %default]"),
    make_option(c("-o", "--output"), type = "character", default = "lai_sim",
                help = "Output prefix [default: %default]"),
    make_option(c("--seed"), type = "integer", default = 42,
                help = "Random seed [default: %default]"),
    make_option(c("--ancestry-codes"), type = "character", default = "0=AFR,1=EUR,2=NAT",
                help = "Ancestry code mapping for MSP file [default: %default]"),
    make_option(c("--list-scenarios"), action = "store_true", default = FALSE,
                help = "List available scenarios")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

set.seed(opt$seed)

# =============================================================================
# LAI-Specific Simulation Scenarios
# =============================================================================

lai_scenarios <- list(

    # Effects only manifest on AFR local ancestry background
    lai_afr_only = list(
        description = "Causal effect only on AFR local ancestry tracts",
        effect_by_lai = function(lai) {
            # lai is 0=AFR, 1=EUR, 2=NAT (or as specified)
            ifelse(lai == 0, 1.0, 0.0)  # Effect only if AFR
        },
        expected_tractor_advantage = "High - standard GWAS diluted by EUR/NAT carriers"
    ),

    # Effects only on EUR background
    lai_eur_only = list(
        description = "Causal effect only on EUR local ancestry tracts",
        effect_by_lai = function(lai) {
            ifelse(lai == 1, 1.0, 0.0)
        },
        expected_tractor_advantage = "Moderate - EUR well-represented but admixed samples dilute"
    ),

    # Effects only on NAT (Native American) background
    lai_nat_only = list(
        description = "Causal effect only on NAT local ancestry tracts",
        effect_by_lai = function(lai) {
            ifelse(lai == 2, 1.0, 0.0)
        },
        expected_tractor_advantage = "Very high - NAT tracts often missed by EUR-centric panels"
    ),

    # Opposite effects by local ancestry (cancel in standard GWAS!)
    lai_opposite = list(
        description = "Opposite effects on AFR vs EUR backgrounds (masks in std GWAS)",
        effect_by_lai = function(lai) {
            case_when(
                lai == 0 ~ 1.0,   # AFR: positive effect
                lai == 1 ~ -1.0, # EUR: negative effect
                TRUE ~ 0.0       # Other: no effect
            )
        },
        expected_tractor_advantage = "Very high - standard GWAS shows ~0 effect"
    ),

    # Gradient effects by local ancestry
    lai_gradient = list(
        description = "Effect size varies continuously with local ancestry dosage",
        effect_by_lai = function(lai) {
            # Stronger effect on AFR background, weaker on EUR
            case_when(
                lai == 0 ~ 1.5,   # AFR
                lai == 1 ~ 0.5,   # EUR
                lai == 2 ~ 1.0,   # NAT (intermediate)
                TRUE ~ 0.0
            )
        },
        expected_tractor_advantage = "Moderate - effect visible but attenuated in std GWAS"
    ),

    # LD-driven scenario (simulates tag SNP vs causal variant)
    lai_ld_difference = list(
        description = "Variant tags causal in AFR (high LD) but not EUR (low LD)",
        effect_by_lai = function(lai) {
            # Tag SNP correlates with causal only in AFR
            # Simulates shorter LD in AFR vs EUR
            ifelse(lai == 0, 1.0, runif(1, 0.1, 0.3))  # EUR: weak due to low LD
        },
        expected_tractor_advantage = "High - explains why EUR GWAS hits don't replicate in AFR"
    )
)

# =============================================================================
# List Scenarios
# =============================================================================

if (opt$`list-scenarios`) {
    cat("\n")
    cat("================================================================================\n")
    cat("LOCAL ANCESTRY-AWARE SIMULATION SCENARIOS (For Tractor GWAS Testing)\n")
    cat("================================================================================\n\n")

    for (name in names(lai_scenarios)) {
        cat(sprintf("  %-20s %s\n", name, lai_scenarios[[name]]$description))
        cat(sprintf("  %20s Tractor advantage: %s\n", "", lai_scenarios[[name]]$expected_tractor_advantage))
        cat("\n")
    }

    cat("Key concept:\n")
    cat("  Standard GWAS: E[Y | G] = β × G\n")
    cat("  Tractor GWAS:  E[Y | G, LAI] = β_AFR × G_AFR + β_EUR × G_EUR + β_NAT × G_NAT\n")
    cat("\n")
    cat("When effects differ by local ancestry, Tractor detects what standard GWAS misses.\n")
    cat("\n")
    quit(status = 0)
}

# =============================================================================
# Validate Inputs
# =============================================================================

cat("================================================================================\n")
cat("LOCAL ANCESTRY-INFORMED PHENOTYPE SIMULATION\n")
cat("================================================================================\n\n")

if (is.null(opt$vcf)) {
    stop("--vcf is required")
}

if (is.null(opt$`lai-msp`)) {
    cat("WARNING: No LAI MSP file provided. Will simulate LAI from population labels.\n")
    cat("         For real benchmarking, run RFMix2 first and provide --lai-msp\n\n")
}

if (!opt$scenario %in% names(lai_scenarios)) {
    stop(sprintf("Unknown scenario: %s. Use --list-scenarios.", opt$scenario))
}

scenario <- lai_scenarios[[opt$scenario]]

cat(sprintf("Scenario:      %s\n", opt$scenario))
cat(sprintf("Description:   %s\n", scenario$description))
cat(sprintf("Tractor gain:  %s\n", scenario$expected_tractor_advantage))
cat("\n")

# Parse ancestry codes
ancestry_map <- setNames(
    sapply(strsplit(strsplit(opt$`ancestry-codes`, ",")[[1]], "="), `[`, 2),
    sapply(strsplit(strsplit(opt$`ancestry-codes`, ",")[[1]], "="), `[`, 1)
)

# =============================================================================
# Read or Simulate Local Ancestry
# =============================================================================

cat("Loading local ancestry information...\n")

# Get samples from VCF
cmd_samples <- sprintf("bcftools query -l %s", opt$vcf)
samples <- fread(cmd = cmd_samples, header = FALSE)$V1
n_samples <- length(samples)

if (!is.null(opt$`lai-msp`) && file.exists(opt$`lai-msp`)) {
    # Read RFMix2 MSP file
    cat(sprintf("  Reading LAI from %s\n", basename(opt$`lai-msp`)))
    msp <- fread(opt$`lai-msp`)

    # MSP format: #chm, spos, epos, sgpos, egpos, n_snps, sample1.0, sample1.1, ...
    msp_samples <- unique(gsub("\\.[01]$", "", names(msp)[7:ncol(msp)]))
    cat(sprintf("  Found %d samples, %d LAI tracts\n", length(msp_samples), nrow(msp)))

} else {
    # Simulate LAI based on global ancestry
    # For real use, should run RFMix2 first!
    cat("  Simulating LAI from global ancestry (for testing only)\n")

    # Infer global ancestry from sample IDs or provide
    # For 1KG: sample naming convention
    get_simulated_lai <- function(sample_id) {
        # Simplified: assign based on 1KG population codes
        # In reality, admixed individuals have mixed LAI across genome
        pop_patterns <- list(
            AFR = c("YRI", "LWK", "GWD", "MSL", "ESN"),
            EUR = c("CEU", "TSI", "FIN", "GBR", "IBS"),
            EAS = c("CHB", "JPT", "CHS", "CDX", "KHV"),
            SAS = c("GIH", "PJL", "BEB", "STU", "ITU"),
            AMR = c("MXL", "PUR", "CLM", "PEL")  # Admixed
        )

        # Check which population the sample belongs to
        for (anc in names(pop_patterns)) {
            for (pop in pop_patterns[[anc]]) {
                if (grepl(pop, sample_id)) {
                    if (anc == "AMR") {
                        # Admixed: simulate mixed LAI
                        return(sample(c(0, 1, 2), 1, prob = c(0.3, 0.5, 0.2)))
                    } else if (anc == "AFR") {
                        return(0)  # AFR
                    } else {
                        return(1)  # EUR/EAS/SAS treated as non-AFR
                    }
                }
            }
        }
        return(1)  # Default EUR
    }

    # Create fake MSP structure
    msp <- data.table(
        `#chm` = 22,
        spos = 1,
        epos = 50000000,
        sgpos = 0,
        egpos = 50,
        n_snps = 1000
    )

    # Add sample columns
    for (s in samples) {
        lai <- get_simulated_lai(s)
        msp[, paste0(s, ".0") := lai]
        msp[, paste0(s, ".1") := lai]  # Same for both haplotypes in this sim
    }

    msp_samples <- samples
}

# =============================================================================
# Select Causal Variants and Extract Genotypes
# =============================================================================

cat("\nSelecting causal variants...\n")

# Get variant list from VCF
cmd <- sprintf("bcftools query -f '%%CHROM\\t%%POS\\t%%ID\\t%%REF\\t%%ALT\\t%%INFO/AF\\n' %s | head -50000",
               opt$vcf)
variants <- tryCatch({
    fread(cmd = cmd, header = FALSE, col.names = c("CHR", "POS", "ID", "REF", "ALT", "AF"))
}, error = function(e) {
    # Simpler query without AF
    cmd2 <- sprintf("bcftools query -f '%%CHROM\\t%%POS\\t%%ID\\t%%REF\\t%%ALT\\n' %s | head -50000", opt$vcf)
    fread(cmd = cmd2, header = FALSE, col.names = c("CHR", "POS", "ID", "REF", "ALT"))
})

if ("AF" %in% names(variants)) {
    variants[, AF := as.numeric(AF)]
    variants <- variants[!is.na(AF) & AF > 0.01 & AF < 0.5]
}

cat(sprintf("  %d eligible variants\n", nrow(variants)))

# Sample causal variants
n_causal <- min(opt$`n-causal`, nrow(variants))
causal_idx <- sample(1:nrow(variants), n_causal)
causal_variants <- variants[causal_idx]

# Extract phased genotypes
cat("  Extracting phased genotypes...\n")
regions_file <- tempfile(fileext = ".txt")
fwrite(causal_variants[, .(CHR, POS, POS)], regions_file, col.names = FALSE, sep = "\t")

cmd <- sprintf("bcftools query -R %s -f '%%CHROM\\t%%POS[\\t%%GT]\\n' %s", regions_file, opt$vcf)
geno <- fread(cmd = cmd, header = FALSE)
setnames(geno, c("CHR", "POS", samples))

unlink(regions_file)

cat(sprintf("  Extracted %d variants\n", nrow(geno)))

# =============================================================================
# Calculate LAI-Aware Genetic Values
# =============================================================================

cat("\nCalculating LAI-aware genetic values...\n")

# For each variant, we need:
# 1. Genotype on each haplotype (from phased VCF)
# 2. Local ancestry on each haplotype (from MSP)
# 3. Apply effect function based on LAI

# Generate base effect sizes
causal_variants[, BASE_BETA := rnorm(n_causal, mean = 0, sd = 1)]

# Function to get LAI at position
get_lai_at_pos <- function(msp_data, chrom, pos, sample, hap) {
    # Find tract containing this position
    tract <- msp_data[`#chm` == chrom & spos <= pos & epos >= pos]

    if (nrow(tract) == 0) {
        return(NA)
    }

    col_name <- paste0(sample, ".", hap)
    if (!col_name %in% names(tract)) {
        return(NA)
    }

    return(tract[[col_name]][1])
}

# Calculate genetic value for each sample
genetic_values <- sapply(samples, function(s) {
    g_value <- 0

    for (i in 1:nrow(geno)) {
        var_chr <- geno$CHR[i]
        var_pos <- geno$POS[i]

        # Get phased genotype
        gt <- geno[[s]][i]
        if (is.na(gt) || gt %in% c(".", "./.", ".|.")) next

        # Parse phased genotype
        alleles <- as.numeric(strsplit(gt, "[|/]")[[1]])
        if (length(alleles) != 2) next

        # Get base effect
        beta <- causal_variants$BASE_BETA[i]

        # For each haplotype, get LAI and apply effect function
        for (hap in 0:1) {
            if (alleles[hap + 1] == 1) {  # Carries alt allele on this haplotype
                # Get LAI for this haplotype
                lai <- get_lai_at_pos(msp, var_chr, var_pos, s, hap)
                if (is.na(lai)) lai <- 1  # Default to EUR if unknown

                # Apply scenario-specific effect
                effect_mult <- scenario$effect_by_lai(lai)
                g_value <- g_value + beta * effect_mult
            }
        }
    }

    return(g_value)
})

results <- data.table(IID = samples, G = genetic_values)

# =============================================================================
# Simulate Phenotype
# =============================================================================

cat("Simulating phenotype...\n")

# Standardize
results[, G_std := (G - mean(G, na.rm = TRUE)) / sd(G, na.rm = TRUE)]
results[G_std == Inf | G_std == -Inf | is.na(G_std), G_std := 0]

# Environmental variance
var_G <- var(results$G_std, na.rm = TRUE)
if (var_G == 0) var_G <- 1
var_E <- var_G * (1 - opt$h2) / opt$h2

results[, E := rnorm(.N, 0, sqrt(var_E))]
results[, P := G_std + E]
results[, PHENO := (P - mean(P)) / sd(P)]
results[, G_TRUE := G_std]

# =============================================================================
# Create Tractor Input Files
# =============================================================================

cat("\nCreating Tractor GWAS input files...\n")

# Tractor needs:
# 1. Phenotype file
# 2. Local ancestry dosage file (from MSP)
# 3. Ancestry-specific genotype files

# Phenotype file
plink_pheno <- results[, .(FID = IID, IID = IID, PHENO)]
fwrite(plink_pheno, paste0(opt$output, ".pheno"), sep = "\t")
cat(sprintf("  %s.pheno\n", opt$output))

# Full results
fwrite(results, paste0(opt$output, "_full.txt"), sep = "\t")

# Causal variants with effects
fwrite(causal_variants, paste0(opt$output, "_causal_variants.txt"), sep = "\t")
cat(sprintf("  %s_causal_variants.txt\n", opt$output))

# LAI dosage summary for Tractor
# Tractor expects: sample, lancAFR, lancEUR, lancNAT (global ancestry proportions)
# And per-SNP local ancestry dosages
if (exists("msp") && nrow(msp) > 0) {
    # Calculate global ancestry proportions from LAI
    lai_dosage <- data.table(IID = samples)

    for (anc_code in names(ancestry_map)) {
        anc_name <- ancestry_map[anc_code]
        anc_code_num <- as.numeric(anc_code)

        # Count fraction of genome on this ancestry
        proportions <- sapply(samples, function(s) {
            hap0_col <- paste0(s, ".0")
            hap1_col <- paste0(s, ".1")

            if (!hap0_col %in% names(msp)) return(NA)

            hap0 <- msp[[hap0_col]]
            hap1 <- msp[[hap1_col]]

            mean(c(hap0 == anc_code_num, hap1 == anc_code_num), na.rm = TRUE)
        })

        lai_dosage[, paste0("prop_", anc_name) := proportions]
    }

    fwrite(lai_dosage, paste0(opt$output, "_lai_dosage.txt"), sep = "\t")
    cat(sprintf("  %s_lai_dosage.txt (for Tractor)\n", opt$output))
}

# =============================================================================
# Summary
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("SIMULATION SUMMARY\n")
cat("================================================================================\n\n")

cat(sprintf("Samples:        %d\n", nrow(results)))
cat(sprintf("Causal variants: %d\n", n_causal))
cat(sprintf("Target h2:      %.2f\n", opt$h2))
cat(sprintf("Achieved h2:    %.3f\n", cor(results$G_TRUE, results$PHENO, use = "complete.obs")^2))

cat("\n")
cat("================================================================================\n")
cat("BENCHMARKING: STANDARD GWAS vs TRACTOR GWAS\n")
cat("================================================================================\n")
cat("\n")
cat("1. Run STANDARD GWAS (ignores local ancestry):\n")
cat(sprintf("   plink2 --bfile imputed --pheno %s.pheno --glm --out gwas_standard\n\n", opt$output))
cat("2. Run TRACTOR GWAS (uses local ancestry):\n")
cat("   # First, run Tractor to create ancestry-specific genotypes:\n")
cat("   python Tractor.py \\\n")
cat("       --vcf imputed.vcf.gz \\\n")
cat("       --msp rfmix_output.msp.tsv \\\n")
cat("       --out tractor_geno\n")
cat("\n")
cat("   # Then run Tractor GWAS:\n")
cat("   python RunTractor.py \\\n")
cat("       --hapdose tractor_geno \\\n")
cat(sprintf("       --phe %s.pheno \\\n", opt$output))
cat("       --method linear \\\n")
cat("       --out gwas_tractor\n")
cat("\n")
cat("3. Compare results:\n")
cat("   - Standard GWAS: Check if causal variants reach significance\n")
cat("   - Tractor GWAS: Check ancestry-specific effects detected\n")
cat(sprintf("   - Expected: Tractor advantage is %s for this scenario\n",
            scenario$expected_tractor_advantage))
cat("\n")
cat("================================================================================\n")
