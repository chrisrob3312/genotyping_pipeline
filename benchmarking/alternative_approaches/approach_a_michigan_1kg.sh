#!/bin/bash
################################################################################
# approach_a_michigan_1kg.sh
#
# APPROACH A: Stringent QC-Before + Michigan Imputation Server (1000G panel)
#
# Philosophy: Clean data thoroughly before imputation (traditional practice)
# Critique: May remove valid variants in admixed populations
#
# Usage:
#   ./approach_a_michigan_1kg.sh \
#       --input /path/to/plink_prefix \
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

INPUT_PLINK=""
OUTPUT_DIR="./results/approach_a_michigan_1kg"
MICHIGAN_TOKEN=""
MICHIGAN_PASSWORD=""
THREADS=4

# QC thresholds (stringent - traditional)
MAF_THRESHOLD=0.01
HWE_PVALUE=1e-6
GENO_THRESHOLD=0.02
MIND_THRESHOLD=0.02

# Post-imputation
R2_THRESHOLD=0.3

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

APPROACH A: Stringent QC-before-imputation with Michigan Server (1000G panel).

Required:
  -i, --input PREFIX         Input PLINK file prefix
  --michigan-token TOKEN     Michigan Imputation Server API token

Optional:
  -o, --output DIR           Output directory (default: ./results/approach_a_michigan_1kg)
  --michigan-password PASS   Michigan password for decryption
  -t, --threads N            Threads (default: 4)
  --maf FLOAT                MAF threshold (default: 0.01)
  --r2 FLOAT                 Post-imputation R² threshold (default: 0.3)
  -h, --help                 Show this help

Note: Uses imputationbot for automated submission.
      Install with: pip install imputationbot

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input) INPUT_PLINK="$2"; shift 2 ;;
        -o|--output) OUTPUT_DIR="$2"; shift 2 ;;
        --michigan-token) MICHIGAN_TOKEN="$2"; shift 2 ;;
        --michigan-password) MICHIGAN_PASSWORD="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --maf) MAF_THRESHOLD="$2"; shift 2 ;;
        --r2) R2_THRESHOLD="$2"; shift 2 ;;
        -h|--help) print_usage; exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$INPUT_PLINK" ]] || [[ -z "$MICHIGAN_TOKEN" ]]; then
    echo "ERROR: --input and --michigan-token are required"
    print_usage
    exit 1
fi

# =============================================================================
# Setup
# =============================================================================

mkdir -p "${OUTPUT_DIR}"/{logs,qc_before,liftover,imputation,qc_after,final}
setup_logging "${OUTPUT_DIR}"

log "=============================================="
log "APPROACH A: Stringent QC + Michigan (1000G)"
log "=============================================="
log "Input: ${INPUT_PLINK}"
log "Output: ${OUTPUT_DIR}"

TOTAL_START=$(date +%s)

# =============================================================================
# STEP 1: PRE-IMPUTATION QC (Stringent)
# =============================================================================

log ""
log "=== STEP 1: Pre-Imputation QC (Stringent) ==="
time_start "STEP1_QC_BEFORE"

run_stringent_qc "${INPUT_PLINK}" "${OUTPUT_DIR}/qc_before/qc" ${THREADS} \
    ${MAF_THRESHOLD} ${HWE_PVALUE} ${GENO_THRESHOLD} ${MIND_THRESHOLD}

time_end "STEP1_QC_BEFORE"
log "  After QC: $(count_variants_samples ${OUTPUT_DIR}/qc_before/qc_qcd)"

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
# STEP 3: SUBMIT TO MICHIGAN (1000G)
# =============================================================================

log ""
log "=== STEP 3: Submit to Michigan Server (1000G panel) ==="
time_start "STEP3_IMPUTATION"

cd "${OUTPUT_DIR}/imputation"

# Split by chromosome
for chr in {1..22}; do
    bcftools view -r chr${chr} ../liftover/lifted_hg38.vcf.gz \
        -Oz -o chr${chr}.vcf.gz
    bcftools index chr${chr}.vcf.gz
done

mkdir -p vcfs_for_imputation
mv chr*.vcf.gz* vcfs_for_imputation/

# Submit via imputationbot
submit_michigan_1kg vcfs_for_imputation michigan_results \
    "${MICHIGAN_TOKEN}" "${MICHIGAN_PASSWORD}" || {
    log "NOTE: Manual submission required"
    log "Upload: ${OUTPUT_DIR}/imputation/vcfs_for_imputation/"
    log "Server: https://imputationserver.sph.umich.edu/"
    log "Panel: 1000G Phase 3 v5"
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

if [[ -d "../imputation/michigan_results" ]]; then
    for vcf in ../imputation/michigan_results/*.vcf.gz; do
        base=$(basename "$vcf" .vcf.gz)
        filter_r2 "$vcf" "${base}_r2filt.vcf.gz" ${R2_THRESHOLD}
    done

    ls *_r2filt.vcf.gz | sort -V > vcf_list.txt
    bcftools concat -f vcf_list.txt -Oz -o imputed_filtered.vcf.gz
    bcftools index imputed_filtered.vcf.gz
else
    log "WARNING: Imputed results not found"
fi

cd "${OUTPUT_DIR}"
time_end "STEP4_R2_FILTER"

# =============================================================================
# STEP 5: BASIC POST-QC
# =============================================================================

log ""
log "=== STEP 5: Basic Post-Imputation QC ==="
time_start "STEP5_QC_AFTER"

cd "${OUTPUT_DIR}/final"

if [[ -f "../qc_after/imputed_filtered.vcf.gz" ]]; then
    plink2 --vcf ../qc_after/imputed_filtered.vcf.gz \
        --make-bed \
        --out approach_a_michigan_1kg_final \
        --threads ${THREADS}

    plink2 --bfile approach_a_michigan_1kg_final \
        --geno 0.05 \
        --mind 0.05 \
        --make-bed \
        --out approach_a_michigan_1kg_qcd \
        --threads ${THREADS}

    log "  Final: $(count_variants_samples approach_a_michigan_1kg_qcd)"
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
log "APPROACH A (Michigan 1000G) Complete!"
log "=============================================="
log "Total time: $(seconds_to_human ${TOTAL_TIME})"
log "Results: ${OUTPUT_DIR}/final/"
