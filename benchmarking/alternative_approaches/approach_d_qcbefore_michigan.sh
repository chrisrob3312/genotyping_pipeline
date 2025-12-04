#!/bin/bash
################################################################################
# approach_d_qcbefore_michigan.sh
#
# APPROACH D: QC-before + Intersect + Michigan (HRC)
#
# Combines traditional QC-before approach with Michigan server.
# Tests impact of using HRC panel vs TOPMed for diverse ancestries.
#
# Pipeline Steps:
#   1. Per-platform QC BEFORE merging
#   2. Intersect variants across platforms
#   3. Merge samples
#   4. Submit to Michigan (HRC panel)
#   5. Traditional R² filtering
#
# Usage:
#   ./approach_d_qcbefore_michigan.sh \
#       --inputs "platform1.bed,platform2.bed" \
#       --output /path/to/output \
#       --michigan-token YOUR_TOKEN
#
################################################################################

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

INPUT_BEDS=""
OUTPUT_DIR="./results/approach_d_qcbefore_michigan"
MICHIGAN_TOKEN=""
PANEL="hrc"
THREADS=4

# QC thresholds
MAF_THRESHOLD=0.01
HWE_PVALUE=1e-6
GENO_THRESHOLD=0.02
MIND_THRESHOLD=0.02
R2_THRESHOLD=0.3

TIMING_LOG=""

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

APPROACH D: Per-platform QC + Intersect + Michigan Imputation.

Required:
  --inputs LIST              Comma-separated PLINK prefixes
  --michigan-token TOKEN     Michigan Imputation Server API token

Optional:
  -o, --output DIR           Output directory
  --panel PANEL              Reference panel: hrc, 1000g-phase3 (default: hrc)
  -t, --threads N            Threads (default: 4)
  -h, --help                 Show this help

Key Difference:
  - QC applied per-platform BEFORE intersecting/merging
  - Uses Michigan/HRC instead of TOPMed
  - Tests whether HRC panel works for diverse cohorts

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --inputs)
            INPUT_BEDS="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --michigan-token)
            MICHIGAN_TOKEN="$2"
            shift 2
            ;;
        --panel)
            PANEL="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -h|--help)
            print_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

if [[ -z "$INPUT_BEDS" ]]; then
    echo "ERROR: Input PLINK files required"
    exit 1
fi

# =============================================================================
# Setup
# =============================================================================

mkdir -p "${OUTPUT_DIR}"/{logs,per_platform_qc,intersect_merge,imputation,final}
TIMING_LOG="${OUTPUT_DIR}/logs/timing.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${OUTPUT_DIR}/logs/pipeline.log"
}

time_start() { echo "$1 START $(date +%s)" >> "${TIMING_LOG}"; }
time_end() { echo "$1 END $(date +%s)" >> "${TIMING_LOG}"; }

log "=============================================="
log "APPROACH D: QC-Before + Intersect + Michigan"
log "=============================================="

TOTAL_START=$(date +%s)
IFS=',' read -ra PLATFORMS <<< "$INPUT_BEDS"

# =============================================================================
# STEP 1: PER-PLATFORM QC (Before Merging)
# =============================================================================

log ""
log "=== STEP 1: Per-Platform QC (Before Merging) ==="
time_start "STEP1_PERPLATFORM_QC"

cd "${OUTPUT_DIR}/per_platform_qc"

for i in "${!PLATFORMS[@]}"; do
    plink="${PLATFORMS[$i]}"
    prefix="platform_${i}"

    log "  Processing platform ${i}..."

    # Call rate filters
    plink2 --bfile "${plink}" \
        --geno ${GENO_THRESHOLD} \
        --mind ${MIND_THRESHOLD} \
        --make-bed \
        --out "${prefix}_callrate" \
        --threads ${THREADS}

    # MAF filter (traditional - before imputation)
    plink2 --bfile "${prefix}_callrate" \
        --maf ${MAF_THRESHOLD} \
        --make-bed \
        --out "${prefix}_maf" \
        --threads ${THREADS}

    # HWE filter
    plink2 --bfile "${prefix}_maf" \
        --hwe ${HWE_PVALUE} midp \
        --make-bed \
        --out "${prefix}_qcd" \
        --threads ${THREADS}

    n_vars=$(wc -l < "${prefix}_qcd.bim")
    n_samp=$(wc -l < "${prefix}_qcd.fam")
    log "    Platform ${i} after QC: ${n_vars} variants, ${n_samp} samples"
done

cd "${OUTPUT_DIR}"
time_end "STEP1_PERPLATFORM_QC"

# =============================================================================
# STEP 2: INTERSECT AND MERGE
# =============================================================================

log ""
log "=== STEP 2: Intersect and Merge Platforms ==="
time_start "STEP2_MERGE"

cd "${OUTPUT_DIR}/intersect_merge"

# Get variant lists
for i in "${!PLATFORMS[@]}"; do
    awk '{print $2}' "../per_platform_qc/platform_${i}_qcd.bim" > "vars_${i}.txt"
done

# Find intersection
cp vars_0.txt common_vars.txt
for i in $(seq 1 $((${#PLATFORMS[@]} - 1))); do
    comm -12 <(sort common_vars.txt) <(sort "vars_${i}.txt") > tmp.txt
    mv tmp.txt common_vars.txt
done

N_COMMON=$(wc -l < common_vars.txt)
log "  Common variants: ${N_COMMON}"

# Extract and merge
for i in "${!PLATFORMS[@]}"; do
    plink2 --bfile "../per_platform_qc/platform_${i}_qcd" \
        --extract common_vars.txt \
        --make-bed \
        --out "platform_${i}_common" \
        --threads ${THREADS}
done

if [[ ${#PLATFORMS[@]} -eq 1 ]]; then
    cp platform_0_common.* merged_final.
else
    for i in $(seq 1 $((${#PLATFORMS[@]} - 1))); do
        echo "platform_${i}_common"
    done > merge_list.txt

    plink2 --bfile platform_0_common \
        --pmerge-list merge_list.txt bfile \
        --make-bed \
        --out merged_final \
        --threads ${THREADS}
fi

cd "${OUTPUT_DIR}"
time_end "STEP2_MERGE"

# =============================================================================
# STEP 3: PREPARE FOR MICHIGAN
# =============================================================================

log ""
log "=== STEP 3: Prepare for Michigan Server ==="
time_start "STEP3_PREPARE"

cd "${OUTPUT_DIR}/imputation"

plink2 --bfile ../intersect_merge/merged_final \
    --export vcf-4.2 bgz \
    --out for_michigan \
    --threads ${THREADS}

bcftools index for_michigan.vcf.gz

# Split by chromosome
for chr in {1..22}; do
    bcftools view -r ${chr} for_michigan.vcf.gz -Oz -o chr${chr}.vcf.gz
    bcftools index chr${chr}.vcf.gz
done

log "  Files ready for Michigan submission: chr*.vcf.gz"
log "  Panel: ${PANEL}"

cd "${OUTPUT_DIR}"
time_end "STEP3_PREPARE"

# =============================================================================
# SUMMARY
# =============================================================================

TOTAL_END=$(date +%s)
TOTAL_TIME=$((TOTAL_END - TOTAL_START))

log ""
log "=============================================="
log "APPROACH D (QC-Before + Michigan) Complete!"
log "=============================================="
log ""
log "Total time: ${TOTAL_TIME} seconds"
log ""
log "Key comparison points:"
log "  - QC applied PER-PLATFORM before merging"
log "  - Intersected variants after QC"
log "  - Uses ${PANEL} panel (vs TOPMed)"
log ""
log "Expected differences from our pipeline:"
log "  - Lost rare variants (MAF filter before imputation)"
log "  - HRC panel has less diverse coverage than TOPMed"
log "  - Per-platform QC may remove different variants per platform"
log ""
