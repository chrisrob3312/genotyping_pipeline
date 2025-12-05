#!/bin/bash
################################################################################
# approach_c_michigan_1kg.sh
#
# APPROACH C: Thorough QC BEFORE Imputation → Merge AFTER Imputation + Michigan 1KG
#
# Philosophy: Apply thorough QC on each platform BEFORE imputation, impute
# each platform separately, then merge AFTER imputation.
#
# KEY COMPARISON:
#   - C vs D: Thorough QC before (C) vs after (D) imputation
#   - Both merge AFTER imputation
#   - Shows effect of QC timing on imputation quality
#
# Workflow:
#   Platform 1 → THOROUGH QC → Ref Align → Impute ─┐
#   Platform 2 → THOROUGH QC → Ref Align → Impute ─┼─→ INTERSECT/MERGE → R² > 0.3 → QC after
#   Platform 3 → THOROUGH QC → Ref Align → Impute ─┘
#                (before impute)          (separately)  (AFTER impute)   (traditional)
#
# Usage:
#   ./approach_c_michigan_1kg.sh \
#       --inputs "platform1.bed,platform2.bed" \
#       --output /path/to/output \
#       --michigan-token YOUR_TOKEN
#
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/common_functions.sh"

# =============================================================================
# Configuration
# =============================================================================

INPUT_BEDS=""
OUTPUT_DIR="./results/approach_c_michigan_1kg"
MICHIGAN_TOKEN=""
MICHIGAN_PASSWORD=""
THREADS=4

# R² filter (traditional)
R2_THRESHOLD=0.3

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

APPROACH C: Thorough QC BEFORE Imputation → Merge AFTER Imputation + Michigan 1000G

KEY COMPARISON: Same as D but with thorough QC BEFORE imputation instead of after.
Tests whether QC timing affects imputation quality.

Workflow:
  1. Thorough QC on each platform (95% call rate, het, relatedness)
  2. Reference alignment (Rayner script) on each platform
  3. Impute each platform SEPARATELY via Michigan (1000G panel)
  4. MERGE across platforms AFTER imputation
  5. R² > 0.3 filter (traditional)
  6. Light post-imputation QC (call rate only)

Required:
  --inputs LIST              Comma-separated PLINK prefixes
  --michigan-token TOKEN     Michigan API token

Optional:
  -o, --output DIR           Output directory
  --michigan-password PASS   Password for decryption
  -t, --threads N            Threads (default: 4)
  --r2 FLOAT                 R² threshold (default: 0.3)
  -h, --help                 Show this help

C vs D comparison: Does QC before or after imputation matter?

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --inputs) INPUT_BEDS="$2"; shift 2 ;;
        -o|--output) OUTPUT_DIR="$2"; shift 2 ;;
        --michigan-token) MICHIGAN_TOKEN="$2"; shift 2 ;;
        --michigan-password) MICHIGAN_PASSWORD="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --r2) R2_THRESHOLD="$2"; shift 2 ;;
        -h|--help) print_usage; exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$INPUT_BEDS" ]] || [[ -z "$MICHIGAN_TOKEN" ]]; then
    echo "ERROR: --inputs and --michigan-token are required"
    exit 1
fi

# =============================================================================
# Setup
# =============================================================================

mkdir -p "${OUTPUT_DIR}"/{logs,per_platform,imputation,merge,qc_after,final}
setup_logging "${OUTPUT_DIR}"

log "=============================================="
log "APPROACH C: Thorough QC BEFORE → Merge AFTER"
log "=============================================="
log "Server: Michigan Imputation Server (1000G panel)"
log "KEY: Thorough QC before imputation (vs D which does QC after)"
log "Inputs: ${INPUT_BEDS}"
log "Output: ${OUTPUT_DIR}"
log "R² threshold: ${R2_THRESHOLD}"

TOTAL_START=$(date +%s)

IFS=',' read -ra PLATFORMS <<< "$INPUT_BEDS"
N_PLATFORMS=${#PLATFORMS[@]}
log "Number of platforms: ${N_PLATFORMS}"

# =============================================================================
# STEP 1: PER-PLATFORM THOROUGH QC + REFERENCE ALIGNMENT
# =============================================================================

log ""
log "=== STEP 1: Per-Platform THOROUGH QC and Reference Alignment ==="
log "  (KEY DIFFERENCE: Thorough QC BEFORE imputation)"
time_start "STEP1_PER_PLATFORM_PREP"

for i in "${!PLATFORMS[@]}"; do
    PLATFORM="${PLATFORMS[$i]}"
    PLATFORM_DIR="${OUTPUT_DIR}/per_platform/platform_${i}"
    mkdir -p "${PLATFORM_DIR}"

    log "  Processing platform ${i}: $(basename ${PLATFORM})"

    # THOROUGH QC (KEY DIFFERENCE from D)
    log "    Running THOROUGH QC before imputation..."
    run_thorough_qc "${PLATFORM}" "${PLATFORM_DIR}/thorough_qcd" ${THREADS}

    log "    After thorough QC: $(count_variants_samples ${PLATFORM_DIR}/thorough_qcd)"

    # Reference alignment (Rayner script)
    log "    Running reference alignment..."
    run_reference_alignment "${PLATFORM_DIR}/thorough_qcd" "${PLATFORM_DIR}/ref_aligned" ${THREADS}

    log "    After ref alignment: $(count_variants_samples ${PLATFORM_DIR}/ref_aligned)"
done

time_end "STEP1_PER_PLATFORM_PREP"

# =============================================================================
# STEP 2: PER-PLATFORM IMPUTATION (SEPARATELY)
# =============================================================================

log ""
log "=== STEP 2: Per-Platform Imputation via Michigan 1000G (SEPARATE) ==="
time_start "STEP2_IMPUTATION"

for i in "${!PLATFORMS[@]}"; do
    PLATFORM_DIR="${OUTPUT_DIR}/per_platform/platform_${i}"
    IMPUTE_DIR="${OUTPUT_DIR}/imputation/platform_${i}"
    mkdir -p "${IMPUTE_DIR}"

    log "  Imputing platform ${i}..."

    cd "${IMPUTE_DIR}"

    # Convert to VCF and liftover to hg38
    plink2 --bfile "${PLATFORM_DIR}/ref_aligned" \
        --export vcf-4.2 bgz \
        --out pre_liftover \
        --threads ${THREADS}

    bcftools index pre_liftover.vcf.gz
    run_liftover_hg38 pre_liftover.vcf.gz lifted_hg38.vcf.gz

    # Split by chromosome
    mkdir -p vcfs_for_imputation
    for chr in {1..22}; do
        bcftools view -r chr${chr} lifted_hg38.vcf.gz \
            -Oz -o vcfs_for_imputation/chr${chr}.vcf.gz
        bcftools index vcfs_for_imputation/chr${chr}.vcf.gz
    done

    # Submit to Michigan Imputation Server (1000G panel)
    submit_michigan_1kg vcfs_for_imputation michigan_results \
        "${MICHIGAN_TOKEN}" "${MICHIGAN_PASSWORD}" || {
        log "NOTE: Manual submission required for platform ${i}"
    }

    cd "${OUTPUT_DIR}"
done

time_end "STEP2_IMPUTATION"

# =============================================================================
# STEP 3: MERGE ACROSS PLATFORMS (AFTER IMPUTATION)
# =============================================================================

log ""
log "=== STEP 3: Merge Across Platforms (AFTER imputation) ==="
time_start "STEP3_MERGE"

cd "${OUTPUT_DIR}/merge"

# Check if imputation results exist
IMPUTED_PLATFORMS=()
for i in "${!PLATFORMS[@]}"; do
    IMPUTE_DIR="${OUTPUT_DIR}/imputation/platform_${i}/michigan_results"
    if [[ -d "${IMPUTE_DIR}" ]]; then
        IMPUTED_PLATFORMS+=("${IMPUTE_DIR}")
        log "  Found results for platform ${i}"
    fi
done

if [[ ${#IMPUTED_PLATFORMS[@]} -gt 0 ]]; then
    # Concat chromosomes for each platform first
    for i in "${!IMPUTED_PLATFORMS[@]}"; do
        IMPUTE_DIR="${IMPUTED_PLATFORMS[$i]}"
        ls "${IMPUTE_DIR}"/*.vcf.gz 2>/dev/null | sort -V > "vcf_list_${i}.txt"

        if [[ -s "vcf_list_${i}.txt" ]]; then
            bcftools concat -f "vcf_list_${i}.txt" -Oz -o "platform_${i}_imputed.vcf.gz"
            bcftools index "platform_${i}_imputed.vcf.gz"
        fi
    done

    # Find common variants across imputed platforms
    log "  Finding common variants..."
    for i in "${!IMPUTED_PLATFORMS[@]}"; do
        if [[ -f "platform_${i}_imputed.vcf.gz" ]]; then
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "platform_${i}_imputed.vcf.gz" > "variants_${i}.txt"
        fi
    done

    # Intersect variants
    if [[ ${#IMPUTED_PLATFORMS[@]} -eq 1 ]]; then
        cp variants_0.txt common_variants.txt
    else
        cp variants_0.txt common_tmp.txt
        for i in $(seq 1 $((${#IMPUTED_PLATFORMS[@]} - 1))); do
            if [[ -f "variants_${i}.txt" ]]; then
                comm -12 <(sort common_tmp.txt) <(sort "variants_${i}.txt") > common_tmp2.txt
                mv common_tmp2.txt common_tmp.txt
            fi
        done
        mv common_tmp.txt common_variants.txt
    fi

    N_COMMON=$(wc -l < common_variants.txt)
    log "  Common variants after imputation: ${N_COMMON}"

    # Create regions file for filtering
    awk '{print $1":"$2"-"$2}' common_variants.txt > common_regions.txt

    # Filter each platform to common variants and merge
    for i in "${!IMPUTED_PLATFORMS[@]}"; do
        if [[ -f "platform_${i}_imputed.vcf.gz" ]]; then
            bcftools view -R common_regions.txt "platform_${i}_imputed.vcf.gz" \
                -Oz -o "platform_${i}_common.vcf.gz"
            bcftools index "platform_${i}_common.vcf.gz"
        fi
    done

    # Merge samples from all platforms
    ls platform_*_common.vcf.gz > merge_list.txt
    bcftools merge -l merge_list.txt -Oz -o merged_imputed.vcf.gz
    bcftools index merged_imputed.vcf.gz

    log "  Merged: $(bcftools query -l merged_imputed.vcf.gz | wc -l) samples"
fi

cd "${OUTPUT_DIR}"
time_end "STEP3_MERGE"

# =============================================================================
# STEP 4: R² FILTER (TRADITIONAL)
# =============================================================================

log ""
log "=== STEP 4: R² Filter (Traditional) ==="
log "  Threshold: R² > ${R2_THRESHOLD}"
time_start "STEP4_R2_FILTER"

cd "${OUTPUT_DIR}/qc_after"

if [[ -f "../merge/merged_imputed.vcf.gz" ]]; then
    # Apply R² filter (traditional)
    filter_r2 "../merge/merged_imputed.vcf.gz" "imputed_r2filt.vcf.gz" ${R2_THRESHOLD}

    # Convert to PLINK
    plink2 --vcf imputed_r2filt.vcf.gz \
        --make-bed \
        --out imputed_r2filt \
        --threads ${THREADS}

    log "  After R² filter: $(count_variants_samples imputed_r2filt)"
fi

cd "${OUTPUT_DIR}"
time_end "STEP4_R2_FILTER"

# =============================================================================
# STEP 5: THOROUGH POST-IMPUTATION QC (QC both before AND after - standard practice)
# =============================================================================

log ""
log "=== STEP 5: Thorough Post-Imputation QC ==="
log "  (Standard practice: thorough QC both BEFORE and AFTER imputation)"
time_start "STEP5_QC_THOROUGH"

cd "${OUTPUT_DIR}/final"

if [[ -f "../qc_after/imputed_r2filt.bed" ]]; then
    # Thorough QC after imputation (standard practice even with pre-imputation QC)
    log "  Applying thorough QC (call rate, het, relatedness)..."

    run_thorough_qc "../qc_after/imputed_r2filt" "approach_c_michigan_1kg" ${THREADS}

    log "  Final: $(count_variants_samples approach_c_michigan_1kg)"

    # Create GWAS and RVAS output tracks
    log "  Creating output tracks..."

    # GWAS track: MAF > 1%
    plink2 --bfile approach_c_michigan_1kg \
        --maf 0.01 \
        --make-bed \
        --out approach_c_michigan_1kg_gwas \
        --threads ${THREADS}

    log "    GWAS track (MAF>1%): $(count_variants_samples approach_c_michigan_1kg_gwas)"

    # RVAS track: All variants (no MAF filter)
    cp approach_c_michigan_1kg.bed approach_c_michigan_1kg_rvas.bed
    cp approach_c_michigan_1kg.bim approach_c_michigan_1kg_rvas.bim
    cp approach_c_michigan_1kg.fam approach_c_michigan_1kg_rvas.fam
    log "    RVAS track (all variants): $(count_variants_samples approach_c_michigan_1kg_rvas)"
fi

cd "${OUTPUT_DIR}"
time_end "STEP5_QC_THOROUGH"

# =============================================================================
# SUMMARY
# =============================================================================

TOTAL_END=$(date +%s)
TOTAL_TIME=$((TOTAL_END - TOTAL_START))

log ""
log "=============================================="
log "APPROACH C (Michigan 1KG) Complete!"
log "=============================================="
log "Total time: $(seconds_to_human ${TOTAL_TIME})"
log ""
log "KEY WORKFLOW: Thorough QC BEFORE AND AFTER imputation, Merge AFTER imputation"
log ""
log "Steps performed:"
log "  1. Per-platform THOROUGH QC (call rate, het, relatedness)"
log "  2. Per-platform reference alignment (Rayner)"
log "  3. Per-platform imputation via Michigan 1000G (SEPARATELY)"
log "  4. MERGE across platforms (AFTER imputation)"
log "  5. R² > ${R2_THRESHOLD} filter (traditional)"
log "  6. THOROUGH post-imputation QC (call rate, het, relatedness)"
log ""
log "C vs D comparison: This does QC BEFORE AND AFTER, D does only AFTER"
log "  → Same merge timing (after impute), C adds pre-imputation QC"
log "  → Shows effect of pre-imputation QC on final quality"
log ""
log "Output tracks:"
log "  GWAS: ${OUTPUT_DIR}/final/approach_c_michigan_1kg_gwas"
log "  RVAS: ${OUTPUT_DIR}/final/approach_c_michigan_1kg_rvas"
