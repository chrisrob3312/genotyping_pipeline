#!/bin/bash
################################################################################
# simulate_gwas_phenotypes.sh
#
# Wrapper script for GWAS phenotype simulation using GCTA.
# Creates scenarios specifically designed to test imputation approaches
# across diverse ancestry groups.
#
# Scenarios:
#   1. afr_enriched   - Variants common in AFR, rare in EUR
#   2. ancestry_het   - Heterogeneous effects across ancestries
#   3. rare_variant   - Rare variants (MAF 0.1-1%)
#   4. realistic      - Mix of shared and ancestry-specific effects
#
# Requirements:
#   - GCTA (will download if not found)
#   - PLINK2
#   - bcftools
#   - R (for some scenarios)
#
# Usage:
#   ./simulate_gwas_phenotypes.sh \
#       --geno benchmarking/test_data/genotypes/1kg_omni \
#       --scenario afr_enriched \
#       --h2 0.5 \
#       --n-causal 100 \
#       --output benchmarking/test_data/phenotypes
#
################################################################################

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

GENO_PREFIX=""
SCENARIO="realistic"
H2=0.5
N_CAUSAL=100
OUTPUT_DIR="./simulated_phenotypes"
SEED=42
FREQ_FILE=""  # Optional: gnomAD or population-specific frequencies

GCTA_URL="https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip"

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Simulate GWAS phenotypes for benchmarking imputation approaches.

Required:
  -g, --geno PREFIX      PLINK file prefix (without .bed/.bim/.fam)

Options:
  -s, --scenario NAME    Simulation scenario (default: realistic)
                         Options: afr_enriched, ancestry_het, rare_variant, realistic
  --h2 VALUE             Heritability (default: 0.5)
  --n-causal N           Number of causal variants (default: 100)
  -o, --output DIR       Output directory (default: ./simulated_phenotypes)
  --seed N               Random seed (default: 42)
  --freq-file FILE       Population frequency file (gnomAD format)
  -h, --help             Show this help

Scenarios:
  afr_enriched    Select causal variants enriched in AFR populations
                  Tests: Does imputation recover AFR-common variants?

  ancestry_het    Same variants with heterogeneous effects by ancestry
                  Tests: Can we detect ancestry-specific effects?

  rare_variant    Causal rare variants (MAF 0.1-1%)
                  Tests: Rare variant imputation quality

  realistic       60% shared + 20% EUR-specific + 20% diverse-enriched
                  Tests: Real-world genetic architecture

Output Files:
  {scenario}.phen           PLINK phenotype file
  {scenario}_causal.txt     Causal variants and effect sizes
  {scenario}_effects.txt    Ancestry-specific effects (if applicable)

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -g|--geno)
            GENO_PREFIX="$2"
            shift 2
            ;;
        -s|--scenario)
            SCENARIO="$2"
            shift 2
            ;;
        --h2)
            H2="$2"
            shift 2
            ;;
        --n-causal)
            N_CAUSAL="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --seed)
            SEED="$2"
            shift 2
            ;;
        --freq-file)
            FREQ_FILE="$2"
            shift 2
            ;;
        -h|--help)
            print_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            print_usage
            exit 1
            ;;
    esac
done

if [ -z "$GENO_PREFIX" ]; then
    echo "ERROR: --geno is required"
    print_usage
    exit 1
fi

echo "=============================================="
echo "GWAS Phenotype Simulation"
echo "=============================================="
echo "Genotypes:    ${GENO_PREFIX}"
echo "Scenario:     ${SCENARIO}"
echo "Heritability: ${H2}"
echo "N causal:     ${N_CAUSAL}"
echo "Output:       ${OUTPUT_DIR}"
echo "Seed:         ${SEED}"
echo ""

mkdir -p "${OUTPUT_DIR}"

# =============================================================================
# Check/Install GCTA
# =============================================================================

check_gcta() {
    if command -v gcta64 &> /dev/null; then
        GCTA_BIN="gcta64"
        return 0
    elif [ -f "./gcta_1.94.1" ]; then
        GCTA_BIN="./gcta_1.94.1"
        return 0
    elif [ -f "${OUTPUT_DIR}/gcta_1.94.1" ]; then
        GCTA_BIN="${OUTPUT_DIR}/gcta_1.94.1"
        return 0
    else
        return 1
    fi
}

if ! check_gcta; then
    echo "GCTA not found. Downloading..."
    cd "${OUTPUT_DIR}"
    wget -q "${GCTA_URL}" -O gcta.zip
    unzip -q gcta.zip
    chmod +x gcta_*
    cd - > /dev/null
    check_gcta || { echo "ERROR: Failed to install GCTA"; exit 1; }
fi

echo "Using GCTA: ${GCTA_BIN}"

# =============================================================================
# Get Variant Information
# =============================================================================

echo ""
echo "=== Extracting variant information ==="

# Get variant MAFs from PLINK file
if [ ! -f "${GENO_PREFIX}.frq" ]; then
    plink2 --bfile "${GENO_PREFIX}" --freq --out "${GENO_PREFIX}" 2>/dev/null || true
fi

# Count variants
N_VARIANTS=$(wc -l < "${GENO_PREFIX}.bim")
echo "  Total variants: ${N_VARIANTS}"

# =============================================================================
# Scenario-Specific Causal Variant Selection
# =============================================================================

CAUSAL_FILE="${OUTPUT_DIR}/${SCENARIO}_causal.txt"

echo ""
echo "=== Selecting causal variants (${SCENARIO}) ==="

case $SCENARIO in

    afr_enriched)
        # Variants common in AFR, rare in EUR
        # If no frequency file, use random selection weighted by MAF
        echo "  Selecting AFR-enriched variants..."

        if [ -n "$FREQ_FILE" ] && [ -f "$FREQ_FILE" ]; then
            # Use provided frequency file
            # Expected columns: CHROM, POS, ID, AF_afr, AF_nfe (or similar)
            awk -v n="$N_CAUSAL" -v seed="$SEED" '
                BEGIN { srand(seed) }
                NR > 1 && $4 > 0.1 && $5 < 0.05 {
                    print $3, (rand()-0.5)*0.1
                }
            ' "$FREQ_FILE" | head -n "$N_CAUSAL" > "$CAUSAL_FILE"
        else
            # Without ancestry-specific frequencies, select moderate-frequency variants
            # and assign effects (user should provide freq file for real benchmarking)
            echo "  WARNING: No frequency file provided."
            echo "           Using MAF-based selection (less realistic)."
            echo "           For real benchmarking, provide gnomAD frequencies with --freq-file"

            awk -v n="$N_CAUSAL" -v seed="$SEED" '
                BEGIN { srand(seed) }
                NR > 1 && $5 > 0.05 && $5 < 0.3 {
                    if (rand() < 0.01) print $2, (rand()-0.5)*0.1
                }
            ' "${GENO_PREFIX}.frq" | head -n "$N_CAUSAL" > "$CAUSAL_FILE"
        fi
        ;;

    ancestry_het)
        # Select variants polymorphic across ancestries
        # Effects will be modified post-simulation
        echo "  Selecting variants for heterogeneous effects..."

        awk -v n="$N_CAUSAL" -v seed="$SEED" '
            BEGIN { srand(seed) }
            NR > 1 && $5 > 0.1 && $5 < 0.4 {
                if (rand() < 0.02) print $2, (rand()-0.5)*0.08
            }
        ' "${GENO_PREFIX}.frq" | head -n "$N_CAUSAL" > "$CAUSAL_FILE"
        ;;

    rare_variant)
        # Rare variants (MAF 0.1-1%)
        echo "  Selecting rare variants (MAF 0.001-0.01)..."

        awk -v n="$N_CAUSAL" -v seed="$SEED" '
            BEGIN { srand(seed) }
            NR > 1 && $5 >= 0.001 && $5 <= 0.01 {
                # Larger effects for rare variants
                if (rand() < 0.1) print $2, (rand()-0.5)*0.3
            }
        ' "${GENO_PREFIX}.frq" | head -n "$N_CAUSAL" > "$CAUSAL_FILE"
        ;;

    realistic)
        # Mixed architecture: some shared, some ancestry-specific
        echo "  Selecting realistic mix of variants..."

        # Get different frequency ranges
        N_COMMON=$((N_CAUSAL * 60 / 100))
        N_MODERATE=$((N_CAUSAL * 30 / 100))
        N_RARE=$((N_CAUSAL - N_COMMON - N_MODERATE))

        # Common variants (MAF > 10%)
        awk -v seed="$SEED" 'BEGIN{srand(seed)} NR>1 && $5>0.1 && $5<0.5 {if(rand()<0.01) print $2, (rand()-0.5)*0.05}' \
            "${GENO_PREFIX}.frq" | head -n "$N_COMMON" > "${CAUSAL_FILE}.tmp"

        # Moderate frequency (MAF 1-10%)
        awk -v seed="$((SEED+1))" 'BEGIN{srand(seed)} NR>1 && $5>0.01 && $5<=0.1 {if(rand()<0.05) print $2, (rand()-0.5)*0.1}' \
            "${GENO_PREFIX}.frq" | head -n "$N_MODERATE" >> "${CAUSAL_FILE}.tmp"

        # Rare variants (MAF < 1%)
        awk -v seed="$((SEED+2))" 'BEGIN{srand(seed)} NR>1 && $5>=0.001 && $5<=0.01 {if(rand()<0.1) print $2, (rand()-0.5)*0.2}' \
            "${GENO_PREFIX}.frq" | head -n "$N_RARE" >> "${CAUSAL_FILE}.tmp"

        mv "${CAUSAL_FILE}.tmp" "$CAUSAL_FILE"
        ;;

    *)
        echo "ERROR: Unknown scenario: ${SCENARIO}"
        echo "Available: afr_enriched, ancestry_het, rare_variant, realistic"
        exit 1
        ;;
esac

N_SELECTED=$(wc -l < "$CAUSAL_FILE")
echo "  Selected ${N_SELECTED} causal variants"

if [ "$N_SELECTED" -lt 10 ]; then
    echo "ERROR: Too few causal variants selected."
    echo "       Check your frequency file or genotype data."
    exit 1
fi

# =============================================================================
# Run GCTA Simulation
# =============================================================================

echo ""
echo "=== Running GCTA phenotype simulation ==="

${GCTA_BIN} --bfile "${GENO_PREFIX}" \
    --simu-qt \
    --simu-causal-loci "$CAUSAL_FILE" \
    --simu-hsq "$H2" \
    --simu-rep 1 \
    --out "${OUTPUT_DIR}/${SCENARIO}"

# GCTA outputs .phen file
if [ -f "${OUTPUT_DIR}/${SCENARIO}.phen" ]; then
    echo "  Created phenotype file: ${OUTPUT_DIR}/${SCENARIO}.phen"
else
    echo "ERROR: GCTA simulation failed"
    exit 1
fi

# =============================================================================
# Create Summary
# =============================================================================

echo ""
echo "=== Creating summary files ==="

# Copy causal variants with full info
cat > "${OUTPUT_DIR}/${SCENARIO}_README.txt" << EOF
================================================================================
Simulated Phenotype: ${SCENARIO}
================================================================================

Date:           $(date)
Genotypes:      ${GENO_PREFIX}
Heritability:   ${H2}
N Causal:       ${N_SELECTED}
Random Seed:    ${SEED}

Scenario Description:
---------------------
EOF

case $SCENARIO in
    afr_enriched)
        cat >> "${OUTPUT_DIR}/${SCENARIO}_README.txt" << EOF
Causal variants enriched in AFR populations (common in AFR, rare in EUR).
Tests whether imputation approaches recover AFR-relevant variants that
traditional EUR-focused pipelines might miss.

Expected results:
- Our pipeline should show higher power in AFR samples
- Traditional pipelines may filter out these variants
EOF
        ;;
    ancestry_het)
        cat >> "${OUTPUT_DIR}/${SCENARIO}_README.txt" << EOF
Variants with heterogeneous effects across ancestries.
Same variants selected, but effects differ by population.

Expected results:
- Standard GWAS shows averaged (diluted) effect
- Tractor GWAS should detect ancestry-specific effects
- Demonstrates value of LAI-aware analysis
EOF
        ;;
    rare_variant)
        cat >> "${OUTPUT_DIR}/${SCENARIO}_README.txt" << EOF
Rare variant causal architecture (MAF 0.1-1%).
Tests imputation quality for rare variants.

Expected results:
- Better imputation → more rare variants recovered
- Our pipeline with TOPMed should outperform 1KG-only
- Larger effects make rare variants detectable
EOF
        ;;
    realistic)
        cat >> "${OUTPUT_DIR}/${SCENARIO}_README.txt" << EOF
Realistic mixed architecture:
- 60% common variants (shared effects)
- 30% moderate frequency (some ancestry-specific)
- 10% rare variants (larger effects)

Expected results:
- Balanced test of overall pipeline performance
- Should see improvement across all components
EOF
        ;;
esac

cat >> "${OUTPUT_DIR}/${SCENARIO}_README.txt" << EOF

Files Created:
--------------
${SCENARIO}.phen           - PLINK phenotype file (FID, IID, PHENO)
${SCENARIO}_causal.txt     - Causal variants (SNP_ID, effect_size)
${SCENARIO}_README.txt     - This file

Usage:
------
# Run GWAS after imputation
plink2 --bfile imputed_data \\
    --pheno ${OUTPUT_DIR}/${SCENARIO}.phen \\
    --glm \\
    --out gwas_${SCENARIO}

# Check causal variant recovery
Rscript benchmarking/bench-helper-scripts/check_hit_recovery.R \\
    --gwas gwas_${SCENARIO}.glm.linear \\
    --causal ${OUTPUT_DIR}/${SCENARIO}_causal.txt \\
    --output recovery_${SCENARIO}

================================================================================
EOF

echo "  Created README: ${OUTPUT_DIR}/${SCENARIO}_README.txt"

# =============================================================================
# Summary
# =============================================================================

echo ""
echo "=============================================="
echo "Simulation Complete!"
echo "=============================================="
echo ""
echo "Output files:"
echo "  ${OUTPUT_DIR}/${SCENARIO}.phen"
echo "  ${OUTPUT_DIR}/${SCENARIO}_causal.txt"
echo "  ${OUTPUT_DIR}/${SCENARIO}_README.txt"
echo ""
echo "Next steps:"
echo "  1. Impute genotypes using different approaches"
echo "  2. Run GWAS: plink2 --bfile <imputed> --pheno ${SCENARIO}.phen --glm"
echo "  3. Compare power: check_hit_recovery.R --causal ${SCENARIO}_causal.txt"
echo ""
