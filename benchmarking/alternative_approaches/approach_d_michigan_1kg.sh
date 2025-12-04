#!/bin/bash
################################################################################
# approach_d_michigan_1kg.sh
#
# APPROACH D: Intersect → QC-After (Verma 2014 / Charon 2021 style) + Michigan 1KG
#
# Philosophy: Intersect for consistency, but delay QC until after imputation
# Based on: Verma 2014 workflow + Charon 2021 QC timing findings
#
# Usage:
#   ./approach_d_michigan_1kg.sh \
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
OUTPUT_DIR="./results/approach_d_michigan_1kg"
MICHIGAN_TOKEN=""
MICHIGAN_PASSWORD=""
THREADS=4

# Minimal pre-imputation QC
PRE_GENO_THRESHOLD=0.10
PRE_MIND_THRESHOLD=0.10

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

APPROACH D: Intersect → QC-After (Verma/Charon style) + Michigan 1000G.

Required:
  --inputs LIST              Comma-separated PLINK prefixes
  --michigan-token TOKEN     Michigan API token

Optional:
  -o, --output DIR           Output directory
  --michigan-password PASS   Password for decryption
  -t, --threads N            Threads (default: 4)
  --r2 FLOAT                 R² threshold (default: 0.3)
  -h, --help                 Show this help

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

mkdir -p "${OUTPUT_DIR}"/{logs,intersect,qc_before,liftover,imputation,qc_after,final}
setup_logging "${OUTPUT_DIR}"

log "=============================================="
log "APPROACH D: Intersect → QC-After + Michigan 1KG"
log "=============================================="
log "Inputs: ${INPUT_BEDS}"
log "Output: ${OUTPUT_DIR}"

TOTAL_START=$(date +%s)

IFS=',' read -ra PLATFORMS <<< "$INPUT_BEDS"
N_PLATFORMS=${#PLATFORMS[@]}

# =============================================================================
# STEP 1: INTERSECT VARIANTS
# =============================================================================

log ""
log "=== STEP 1: Intersect Variants ==="
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
# STEP 2: MINIMAL PRE-IMPUTATION QC
# =============================================================================

log ""
log "=== STEP 2: Minimal Pre-Imputation QC ==="
time_start "STEP2_QC_MINIMAL"

cd "${OUTPUT_DIR}/qc_before"

plink2 --bfile ../intersect/merged \
    --geno ${PRE_GENO_THRESHOLD} \
    --mind ${PRE_MIND_THRESHOLD} \
    --make-bed \
    --out minimal_qcd \
    --threads ${THREADS}

cd "${OUTPUT_DIR}"
time_end "STEP2_QC_MINIMAL"

# =============================================================================
# STEP 3: LIFTOVER AND SUBMIT
# =============================================================================

log ""
log "=== STEP 3: Liftover and Submit to Michigan 1KG ==="
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

submit_michigan_1kg vcfs_for_imputation michigan_results \
    "${MICHIGAN_TOKEN}" "${MICHIGAN_PASSWORD}" || {
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

if [[ -d "../imputation/michigan_results" ]]; then
    for vcf in ../imputation/michigan_results/*.vcf.gz; do
        base=$(basename "$vcf" .vcf.gz)
        filter_r2 "$vcf" "${base}_r2filt.vcf.gz" ${R2_THRESHOLD}
    done

    ls *_r2filt.vcf.gz 2>/dev/null | sort -V > vcf_list.txt
    if [[ -s vcf_list.txt ]]; then
        bcftools concat -f vcf_list.txt -Oz -o imputed_filtered.vcf.gz
        plink2 --vcf imputed_filtered.vcf.gz --make-bed --out imputed_r2filt --threads ${THREADS}
    fi
fi

cd "${OUTPUT_DIR}"
time_end "STEP4_R2_FILTER"

# =============================================================================
# STEP 5: THOROUGH POST-IMPUTATION QC
# =============================================================================

log ""
log "=== STEP 5: Thorough Post-Imputation QC ==="
time_start "STEP5_QC_THOROUGH"

cd "${OUTPUT_DIR}/final"

if [[ -f "../qc_after/imputed_r2filt.bed" ]]; then
    run_post_imputation_qc "../qc_after/imputed_r2filt" "approach_d_michigan_1kg" ${THREADS} \
        ${POST_MAF_THRESHOLD} ${POST_HWE_PVALUE} ${POST_GENO_THRESHOLD} ${POST_MIND_THRESHOLD}

    log "  Final: $(count_variants_samples approach_d_michigan_1kg)"
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
log "APPROACH D (Michigan 1KG) Complete!"
log "=============================================="
log "Total time: $(seconds_to_human ${TOTAL_TIME})"
log "Results: ${OUTPUT_DIR}/final/"
