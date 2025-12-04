#!/bin/bash
################################################################################
# approach_d_topmed.sh
#
# APPROACH D: Intersect → QC-After (Verma 2014 / Charon 2021 style) + TOPMed
#
# Philosophy: Intersect for consistency, but delay QC until after imputation
# Based on: Verma 2014 workflow + Charon 2021 QC timing findings
#
# Key: Minimal pre-imputation QC, thorough POST-imputation QC
#
# Usage:
#   ./approach_d_topmed.sh \
#       --inputs "platform1.bed,platform2.bed" \
#       --output /path/to/output \
#       --topmed-token YOUR_TOKEN
#
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/common_functions.sh"

# =============================================================================
# Configuration
# =============================================================================

INPUT_BEDS=""
OUTPUT_DIR="./results/approach_d_topmed"
TOPMED_TOKEN=""
TOPMED_PASSWORD=""
THREADS=4

# Minimal pre-imputation QC (Charon 2021 recommendation)
PRE_GENO_THRESHOLD=0.10   # Lenient 90%
PRE_MIND_THRESHOLD=0.10   # Lenient 90%

# Thorough post-imputation QC
R2_THRESHOLD=0.3
POST_MAF_THRESHOLD=0.01
POST_HWE_PVALUE=1e-6
POST_GENO_THRESHOLD=0.02
POST_MIND_THRESHOLD=0.02

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

APPROACH D: Intersect → QC-After (Verma/Charon style) + TOPMed.

Based on:
  - Verma 2014: Multi-platform merging workflow
  - Charon 2021: Post-imputation filtering often superior

Required:
  --inputs LIST            Comma-separated PLINK prefixes
  --topmed-token TOKEN     TOPMed API token

Optional:
  -o, --output DIR         Output directory
  --topmed-password PASS   Password for decryption
  -t, --threads N          Threads (default: 4)
  --r2 FLOAT               R² threshold (default: 0.3)
  -h, --help               Show this help

Key difference from Approach C:
  - Minimal QC BEFORE imputation
  - Thorough QC AFTER imputation (per Charon 2021)

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --inputs) INPUT_BEDS="$2"; shift 2 ;;
        -o|--output) OUTPUT_DIR="$2"; shift 2 ;;
        --topmed-token) TOPMED_TOKEN="$2"; shift 2 ;;
        --topmed-password) TOPMED_PASSWORD="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --r2) R2_THRESHOLD="$2"; shift 2 ;;
        -h|--help) print_usage; exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$INPUT_BEDS" ]] || [[ -z "$TOPMED_TOKEN" ]]; then
    echo "ERROR: --inputs and --topmed-token are required"
    exit 1
fi

# =============================================================================
# Setup
# =============================================================================

mkdir -p "${OUTPUT_DIR}"/{logs,intersect,qc_before,liftover,imputation,qc_after,final}
setup_logging "${OUTPUT_DIR}"

log "=============================================="
log "APPROACH D: Intersect → QC-After + TOPMed"
log "=============================================="
log "Based on: Verma 2014 + Charon 2021"
log "Inputs: ${INPUT_BEDS}"
log "Output: ${OUTPUT_DIR}"

TOTAL_START=$(date +%s)

IFS=',' read -ra PLATFORMS <<< "$INPUT_BEDS"
N_PLATFORMS=${#PLATFORMS[@]}
log "Number of platforms: ${N_PLATFORMS}"

# =============================================================================
# STEP 1: INTERSECT VARIANTS (same as Approach C)
# =============================================================================

log ""
log "=== STEP 1: Intersect Variants Across Platforms ==="
time_start "STEP1_INTERSECT"

cd "${OUTPUT_DIR}/intersect"

for i in "${!PLATFORMS[@]}"; do
    cut -f2 "${PLATFORMS[$i]}.bim" > "variants_${i}.txt"
done

if [[ ${N_PLATFORMS} -eq 1 ]]; then
    cp variants_0.txt common_variants.txt
else
    cp variants_0.txt common_tmp.txt
    for i in $(seq 1 $((N_PLATFORMS - 1))); do
        comm -12 <(sort common_tmp.txt) <(sort "variants_${i}.txt") > common_tmp2.txt
        mv common_tmp2.txt common_tmp.txt
    done
    mv common_tmp.txt common_variants.txt
fi

N_COMMON=$(wc -l < common_variants.txt)
log "  Common variants: ${N_COMMON}"

for i in "${!PLATFORMS[@]}"; do
    plink2 --bfile "${PLATFORMS[$i]}" \
        --extract common_variants.txt \
        --make-bed \
        --out "platform_${i}" \
        --threads ${THREADS}
done

if [[ ${N_PLATFORMS} -gt 1 ]]; then
    for i in $(seq 1 $((N_PLATFORMS - 1))); do
        echo "platform_${i}"
    done > merge_list.txt

    plink2 --bfile platform_0 \
        --pmerge-list merge_list.txt bfile \
        --make-bed \
        --out merged \
        --threads ${THREADS}
else
    cp platform_0.* merged.*
fi

cd "${OUTPUT_DIR}"
time_end "STEP1_INTERSECT"

# =============================================================================
# STEP 2: MINIMAL PRE-IMPUTATION QC (KEY DIFFERENCE FROM APPROACH C)
# =============================================================================

log ""
log "=== STEP 2: Minimal Pre-Imputation QC (Charon 2021) ==="
time_start "STEP2_QC_MINIMAL"

cd "${OUTPUT_DIR}/qc_before"

# Only call rate filters - NO MAF, NO HWE (save for after imputation)
log "  Minimal filters: call rate only"

plink2 --bfile ../intersect/merged \
    --geno ${PRE_GENO_THRESHOLD} \
    --make-bed \
    --out step1_geno \
    --threads ${THREADS}

plink2 --bfile step1_geno \
    --mind ${PRE_MIND_THRESHOLD} \
    --make-bed \
    --out minimal_qcd \
    --threads ${THREADS}

cd "${OUTPUT_DIR}"
time_end "STEP2_QC_MINIMAL"

log "  After minimal QC: $(count_variants_samples ${OUTPUT_DIR}/qc_before/minimal_qcd)"

# =============================================================================
# STEP 3: LIFTOVER AND SUBMIT TO TOPMED
# =============================================================================

log ""
log "=== STEP 3: Liftover and Submit to TOPMed ==="
time_start "STEP3_IMPUTATION"

cd "${OUTPUT_DIR}/liftover"

plink2 --bfile ../qc_before/minimal_qcd \
    --export vcf-4.2 bgz \
    --out pre_liftover \
    --threads ${THREADS}

bcftools index pre_liftover.vcf.gz
run_liftover_hg38 pre_liftover.vcf.gz lifted_hg38.vcf.gz

cd "${OUTPUT_DIR}/imputation"

for chr in {1..22}; do
    bcftools view -r chr${chr} ../liftover/lifted_hg38.vcf.gz \
        -Oz -o chr${chr}.vcf.gz
    bcftools index chr${chr}.vcf.gz
done

mkdir -p vcfs_for_imputation
mv chr*.vcf.gz* vcfs_for_imputation/

submit_topmed vcfs_for_imputation topmed_results \
    "${TOPMED_TOKEN}" "${TOPMED_PASSWORD}" || {
    log "NOTE: Manual submission required"
}

cd "${OUTPUT_DIR}"
time_end "STEP3_IMPUTATION"

# =============================================================================
# STEP 4: R² FILTER
# =============================================================================

log ""
log "=== STEP 4: Post-Imputation R² Filtering ==="
time_start "STEP4_R2_FILTER"

cd "${OUTPUT_DIR}/qc_after"

if [[ -d "../imputation/topmed_results" ]]; then
    for vcf in ../imputation/topmed_results/*.vcf.gz; do
        base=$(basename "$vcf" .vcf.gz)
        filter_r2 "$vcf" "${base}_r2filt.vcf.gz" ${R2_THRESHOLD}
    done

    ls *_r2filt.vcf.gz 2>/dev/null | sort -V > vcf_list.txt
    if [[ -s vcf_list.txt ]]; then
        bcftools concat -f vcf_list.txt -Oz -o imputed_filtered.vcf.gz
        bcftools index imputed_filtered.vcf.gz

        plink2 --vcf imputed_filtered.vcf.gz \
            --make-bed \
            --out imputed_r2filt \
            --threads ${THREADS}
    fi
fi

cd "${OUTPUT_DIR}"
time_end "STEP4_R2_FILTER"

# =============================================================================
# STEP 5: THOROUGH POST-IMPUTATION QC (KEY FROM CHARON 2021)
# =============================================================================

log ""
log "=== STEP 5: Thorough Post-Imputation QC (Charon 2021) ==="
time_start "STEP5_QC_THOROUGH"

cd "${OUTPUT_DIR}/final"

if [[ -f "../qc_after/imputed_r2filt.bed" ]]; then
    log "  Applying MAF, HWE, call rate, and relatedness filters..."

    run_post_imputation_qc "../qc_after/imputed_r2filt" "approach_d_topmed" ${THREADS} \
        ${POST_MAF_THRESHOLD} ${POST_HWE_PVALUE} ${POST_GENO_THRESHOLD} ${POST_MIND_THRESHOLD}

    log "  Final: $(count_variants_samples approach_d_topmed)"
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
log "APPROACH D (TOPMed) Complete!"
log "=============================================="
log "Total time: $(seconds_to_human ${TOTAL_TIME})"
log ""
log "Key differences from Approach C:"
log "  - Minimal QC before imputation (call rate only)"
log "  - Thorough QC AFTER imputation (MAF, HWE, relatedness)"
log "  - Based on Charon 2021 finding that post-QC is superior"
log ""
log "Results: ${OUTPUT_DIR}/final/"
