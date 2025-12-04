#!/bin/bash
################################################################################
# approach_b_allofus.sh
#
# APPROACH B: Minimal QC-Before + All of Us (Southam 2011 Recommendation)
#
# Philosophy: Let imputation handle variant filtering; thorough QC AFTER
# Based on: Southam et al. 2011 finding that stringent pre-QC hurts imputation
#
# Usage:
#   ./approach_b_allofus.sh \
#       --input /path/to/plink_prefix \
#       --output /path/to/output
#
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/common_functions.sh"

# =============================================================================
# Configuration
# =============================================================================

INPUT_PLINK=""
OUTPUT_DIR="./results/approach_b_allofus"
THREADS=4

# Minimal QC thresholds (Southam 2011 style)
GENO_THRESHOLD=0.10
MIND_THRESHOLD=0.10

# Post-imputation (thorough)
R2_THRESHOLD=0.3
MAF_THRESHOLD=0.01
HWE_PVALUE=1e-6

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

APPROACH B: Minimal QC-before-imputation (Southam 2011) + All of Us.

Required:
  -i, --input PREFIX       Input PLINK file prefix

Optional:
  -o, --output DIR         Output directory
  -t, --threads N          Threads (default: 4)
  --r2 FLOAT               Post-imputation R² threshold (default: 0.3)
  -h, --help               Show this help

Note: Requires terralab CLI. Run 'terralab login' first.

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input) INPUT_PLINK="$2"; shift 2 ;;
        -o|--output) OUTPUT_DIR="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --r2) R2_THRESHOLD="$2"; shift 2 ;;
        -h|--help) print_usage; exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$INPUT_PLINK" ]]; then
    echo "ERROR: --input is required"
    print_usage
    exit 1
fi

# =============================================================================
# Setup
# =============================================================================

mkdir -p "${OUTPUT_DIR}"/{logs,qc_before,liftover,imputation,qc_after,final}
setup_logging "${OUTPUT_DIR}"

log "=============================================="
log "APPROACH B: Minimal QC + All of Us (Southam 2011)"
log "=============================================="
log "Input: ${INPUT_PLINK}"
log "Output: ${OUTPUT_DIR}"

TOTAL_START=$(date +%s)

# =============================================================================
# STEP 1: MINIMAL PRE-IMPUTATION QC
# =============================================================================

log ""
log "=== STEP 1: Minimal Pre-Imputation QC ==="
time_start "STEP1_QC_BEFORE"

run_minimal_qc "${INPUT_PLINK}" "${OUTPUT_DIR}/qc_before/qc" ${THREADS} \
    ${GENO_THRESHOLD} ${MIND_THRESHOLD}

time_end "STEP1_QC_BEFORE"
log "  After minimal QC: $(count_variants_samples ${OUTPUT_DIR}/qc_before/qc_qcd)"

# =============================================================================
# STEP 2: LIFTOVER TO hg38
# =============================================================================

log ""
log "=== STEP 2: Liftover to hg38 ==="
time_start "STEP2_LIFTOVER"

cd "${OUTPUT_DIR}/liftover"
plink2 --bfile ../qc_before/qc_qcd \
    --export vcf-4.2 bgz \
    --out pre_liftover \
    --threads ${THREADS}

bcftools index pre_liftover.vcf.gz
run_liftover_hg38 pre_liftover.vcf.gz lifted_hg38.vcf.gz

cd "${OUTPUT_DIR}"
time_end "STEP2_LIFTOVER"

# =============================================================================
# STEP 3: SUBMIT TO ALL OF US
# =============================================================================

log ""
log "=== STEP 3: Submit to All of Us Server ==="
time_start "STEP3_IMPUTATION"

cd "${OUTPUT_DIR}/imputation"

for chr in {1..22}; do
    bcftools view -r chr${chr} ../liftover/lifted_hg38.vcf.gz \
        -Oz -o chr${chr}.vcf.gz
    bcftools index chr${chr}.vcf.gz
done

mkdir -p vcfs_for_imputation
mv chr*.vcf.gz* vcfs_for_imputation/

submit_allofus vcfs_for_imputation allofus_results || {
    log "NOTE: Manual submission required via terralab"
}

cd "${OUTPUT_DIR}"
time_end "STEP3_IMPUTATION"

# =============================================================================
# STEP 4: POST-IMPUTATION R² FILTER
# =============================================================================

log ""
log "=== STEP 4: Post-Imputation R² Filtering ==="
time_start "STEP4_R2_FILTER"

cd "${OUTPUT_DIR}/qc_after"

if [[ -d "../imputation/allofus_results" ]]; then
    for vcf in ../imputation/allofus_results/*.vcf.gz; do
        base=$(basename "$vcf" .vcf.gz)
        filter_r2 "$vcf" "${base}_r2filt.vcf.gz" ${R2_THRESHOLD}
    done

    ls *_r2filt.vcf.gz | sort -V > vcf_list.txt
    bcftools concat -f vcf_list.txt -Oz -o imputed_filtered.vcf.gz
    bcftools index imputed_filtered.vcf.gz

    plink2 --vcf imputed_filtered.vcf.gz \
        --make-bed \
        --out imputed_r2filt \
        --threads ${THREADS}
fi

cd "${OUTPUT_DIR}"
time_end "STEP4_R2_FILTER"

# =============================================================================
# STEP 5: THOROUGH POST-IMPUTATION QC
# =============================================================================

log ""
log "=== STEP 5: Thorough Post-Imputation QC ==="
time_start "STEP5_QC_AFTER"

cd "${OUTPUT_DIR}/final"

if [[ -f "../qc_after/imputed_r2filt.bed" ]]; then
    run_post_imputation_qc "../qc_after/imputed_r2filt" "approach_b_allofus" ${THREADS} \
        ${MAF_THRESHOLD} ${HWE_PVALUE} 0.02 0.02

    log "  Final: $(count_variants_samples approach_b_allofus)"
fi

cd "${OUTPUT_DIR}"
time_end "STEP5_QC_AFTER"

# =============================================================================
# SUMMARY
# =============================================================================

TOTAL_END=$(date +%s)
TOTAL_TIME=$((TOTAL_END - TOTAL_START))

log ""
log "=============================================="
log "APPROACH B (All of Us) Complete!"
log "=============================================="
log "Total time: $(seconds_to_human ${TOTAL_TIME})"
log "Results: ${OUTPUT_DIR}/final/"
