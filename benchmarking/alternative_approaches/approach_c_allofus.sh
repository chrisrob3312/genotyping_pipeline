#!/bin/bash
################################################################################
# approach_c_allofus.sh
#
# APPROACH C: Intersect-First + QC-Before + All of Us
#
# Philosophy: Ensure all platforms have identical variants before processing
# Critique: Loses platform-specific variants that could be informative
#
# Usage:
#   ./approach_c_allofus.sh \
#       --inputs "platform1.bed,platform2.bed" \
#       --output /path/to/output
#
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/common_functions.sh"

# =============================================================================
# Configuration
# =============================================================================

INPUT_BEDS=""
OUTPUT_DIR="./results/approach_c_allofus"
THREADS=4

# QC thresholds
MAF_THRESHOLD=0.01
HWE_PVALUE=1e-6
GENO_THRESHOLD=0.02
MIND_THRESHOLD=0.02
R2_THRESHOLD=0.3

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

APPROACH C: Intersect-first + QC-before + All of Us.

Required:
  --inputs LIST            Comma-separated PLINK prefixes

Optional:
  -o, --output DIR         Output directory
  -t, --threads N          Threads (default: 4)
  --maf FLOAT              MAF threshold (default: 0.01)
  --r2 FLOAT               R² threshold (default: 0.3)
  -h, --help               Show this help

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --inputs) INPUT_BEDS="$2"; shift 2 ;;
        -o|--output) OUTPUT_DIR="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --maf) MAF_THRESHOLD="$2"; shift 2 ;;
        --r2) R2_THRESHOLD="$2"; shift 2 ;;
        -h|--help) print_usage; exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$INPUT_BEDS" ]]; then
    echo "ERROR: --inputs is required"
    exit 1
fi

# =============================================================================
# Setup
# =============================================================================

mkdir -p "${OUTPUT_DIR}"/{logs,intersect,qc_before,liftover,imputation,qc_after,final}
setup_logging "${OUTPUT_DIR}"

log "=============================================="
log "APPROACH C: Intersect-First + All of Us"
log "=============================================="
log "Inputs: ${INPUT_BEDS}"
log "Output: ${OUTPUT_DIR}"

TOTAL_START=$(date +%s)

IFS=',' read -ra PLATFORMS <<< "$INPUT_BEDS"
N_PLATFORMS=${#PLATFORMS[@]}
log "Number of platforms: ${N_PLATFORMS}"

# =============================================================================
# STEP 1: INTERSECT VARIANTS
# =============================================================================

log ""
log "=== STEP 1: Intersect Variants Across Platforms ==="
time_start "STEP1_INTERSECT"

cd "${OUTPUT_DIR}/intersect"

# Extract and intersect variant lists
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

# Extract and merge
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
# STEP 2: PRE-IMPUTATION QC
# =============================================================================

log ""
log "=== STEP 2: Pre-Imputation QC ==="
time_start "STEP2_QC"

run_stringent_qc "${OUTPUT_DIR}/intersect/merged" "${OUTPUT_DIR}/qc_before/qc" \
    ${THREADS} ${MAF_THRESHOLD} ${HWE_PVALUE} ${GENO_THRESHOLD} ${MIND_THRESHOLD}

time_end "STEP2_QC"
log "  After QC: $(count_variants_samples ${OUTPUT_DIR}/qc_before/qc_qcd)"

# =============================================================================
# STEP 3: LIFTOVER AND SUBMIT
# =============================================================================

log ""
log "=== STEP 3: Liftover and Submit to All of Us ==="
time_start "STEP3_IMPUTATION"

cd "${OUTPUT_DIR}/liftover"

plink2 --bfile ../qc_before/qc_qcd \
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

submit_allofus vcfs_for_imputation allofus_results || {
    log "NOTE: Manual submission required"
}

cd "${OUTPUT_DIR}"
time_end "STEP3_IMPUTATION"

# =============================================================================
# STEP 4-5: POST-IMPUTATION QC
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

log ""
log "=== STEP 5: Final QC ==="
time_start "STEP5_FINAL"

cd "${OUTPUT_DIR}/final"

if [[ -f "../qc_after/imputed_r2filt.bed" ]]; then
    plink2 --bfile ../qc_after/imputed_r2filt \
        --geno 0.05 \
        --mind 0.05 \
        --make-bed \
        --out approach_c_allofus \
        --threads ${THREADS}

    log "  Final: $(count_variants_samples approach_c_allofus)"
fi

cd "${OUTPUT_DIR}"
time_end "STEP5_FINAL"

# =============================================================================
# SUMMARY
# =============================================================================

TOTAL_END=$(date +%s)
TOTAL_TIME=$((TOTAL_END - TOTAL_START))

log ""
log "=============================================="
log "APPROACH C (All of Us) Complete!"
log "=============================================="
log "Total time: $(seconds_to_human ${TOTAL_TIME})"
log "Results: ${OUTPUT_DIR}/final/"
