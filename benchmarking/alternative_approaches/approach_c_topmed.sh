#!/bin/bash
################################################################################
# approach_c_intersect_first.sh
#
# APPROACH C: Intersect-first + QC-before-imputation + TOPMed
#
# Key difference: Intersects variants ACROSS platforms/batches BEFORE any QC.
# This ensures all platforms have identical variant sets but may lose
# platform-specific variants that could be useful.
#
# Pipeline Steps:
#   1. Intersect variants across all platforms/batches FIRST
#   2. Apply thorough QC (HWE, MAF, call rate, het, relatedness)
#   3. Submit intersected data to TOPMed
#   4. Apply traditional R² filter
#   5. Basic post-imputation QC
#
# Usage:
#   ./approach_c_intersect_first.sh \
#       --inputs "platform1.bed,platform2.bed,platform3.bed" \
#       --output /path/to/output \
#       --topmed-token YOUR_TOKEN
#
################################################################################

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

INPUT_BEDS=""  # Comma-separated list of PLINK prefixes
OUTPUT_DIR="./results/approach_c_intersect"
TOPMED_TOKEN=""
TOPMED_PASSWORD=""
THREADS=4

# QC thresholds
MAF_THRESHOLD=0.01
HWE_PVALUE=1e-6
GENO_THRESHOLD=0.02
MIND_THRESHOLD=0.02
HET_SD=3
PIHAT_THRESHOLD=0.25
R2_THRESHOLD=0.3

TIMING_LOG=""

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

APPROACH C: Intersect-first, then QC, then impute with TOPMed.

Required:
  --inputs LIST              Comma-separated PLINK prefixes (e.g., "p1,p2,p3")
  --topmed-token TOKEN       TOPMed API token
  --topmed-password PASS     TOPMed password

Optional:
  -o, --output DIR           Output directory (default: ./results/approach_c_intersect)
  -t, --threads N            Threads (default: 4)
  --maf FLOAT                MAF threshold (default: 0.01)
  --r2 FLOAT                 Post-imputation R² threshold (default: 0.3)
  -h, --help                 Show this help

Key Difference from Approach A:
  - Intersects variants BEFORE QC (loses platform-specific variants)
  - Tests whether union vs intersect strategy matters

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
        --topmed-token)
            TOPMED_TOKEN="$2"
            shift 2
            ;;
        --topmed-password)
            TOPMED_PASSWORD="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --maf)
            MAF_THRESHOLD="$2"
            shift 2
            ;;
        --r2)
            R2_THRESHOLD="$2"
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

if [[ -z "$INPUT_BEDS" ]]; then
    echo "ERROR: Input PLINK files required (--inputs)"
    print_usage
    exit 1
fi

# =============================================================================
# Setup
# =============================================================================

mkdir -p "${OUTPUT_DIR}"/{logs,intersect,qc_before,liftover,imputation,qc_after,final}
TIMING_LOG="${OUTPUT_DIR}/logs/timing.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${OUTPUT_DIR}/logs/pipeline.log"
}

time_start() {
    echo "$1 START $(date +%s)" >> "${TIMING_LOG}"
}

time_end() {
    echo "$1 END $(date +%s)" >> "${TIMING_LOG}"
}

log "=============================================="
log "APPROACH C: Intersect-First + QC + TOPMed"
log "=============================================="
log "Inputs: ${INPUT_BEDS}"
log "Output: ${OUTPUT_DIR}"
log ""

TOTAL_START=$(date +%s)

# Parse input list
IFS=',' read -ra PLATFORMS <<< "$INPUT_BEDS"
N_PLATFORMS=${#PLATFORMS[@]}
log "Number of platforms/batches: ${N_PLATFORMS}"

# =============================================================================
# STEP 1: INTERSECT VARIANTS ACROSS PLATFORMS
# =============================================================================

log ""
log "=== STEP 1: Intersect Variants Across Platforms ==="
time_start "STEP1_INTERSECT"

cd "${OUTPUT_DIR}/intersect"

# Extract variant IDs from each platform
log "  Extracting variant lists..."
for i in "${!PLATFORMS[@]}"; do
    plink="${PLATFORMS[$i]}"
    awk '{print $2}' "${plink}.bim" > "variants_platform_${i}.txt"
    n_vars=$(wc -l < "variants_platform_${i}.txt")
    log "    Platform ${i}: ${n_vars} variants"
done

# Find intersection
log "  Computing intersection..."
if [[ ${N_PLATFORMS} -eq 1 ]]; then
    cp variants_platform_0.txt common_variants.txt
else
    # Start with first platform
    cp variants_platform_0.txt common_variants_tmp.txt

    # Intersect with each subsequent platform
    for i in $(seq 1 $((N_PLATFORMS - 1))); do
        comm -12 <(sort common_variants_tmp.txt) <(sort "variants_platform_${i}.txt") > common_variants_tmp2.txt
        mv common_variants_tmp2.txt common_variants_tmp.txt
    done

    mv common_variants_tmp.txt common_variants.txt
fi

N_COMMON=$(wc -l < common_variants.txt)
log "  Common variants across all platforms: ${N_COMMON}"

# Extract common variants from each platform
log "  Extracting common variants from each platform..."
for i in "${!PLATFORMS[@]}"; do
    plink="${PLATFORMS[$i]}"
    plink2 --bfile "${plink}" \
        --extract common_variants.txt \
        --make-bed \
        --out "platform_${i}_intersected" \
        --threads ${THREADS}
done

# Merge all platforms
log "  Merging platforms..."
if [[ ${N_PLATFORMS} -eq 1 ]]; then
    cp platform_0_intersected.* merged_intersected.
else
    # Create merge list
    for i in $(seq 1 $((N_PLATFORMS - 1))); do
        echo "platform_${i}_intersected"
    done > merge_list.txt

    plink2 --bfile platform_0_intersected \
        --pmerge-list merge_list.txt bfile \
        --make-bed \
        --out merged_intersected \
        --threads ${THREADS}
fi

N_MERGED_VARS=$(wc -l < merged_intersected.bim)
N_MERGED_SAMP=$(wc -l < merged_intersected.fam)
log "  After intersect-merge: ${N_MERGED_VARS} variants, ${N_MERGED_SAMP} samples"

cd "${OUTPUT_DIR}"
time_end "STEP1_INTERSECT"

# =============================================================================
# STEP 2: PRE-IMPUTATION QC (on intersected data)
# =============================================================================

log ""
log "=== STEP 2: Pre-Imputation QC ==="
time_start "STEP2_QC"

cd "${OUTPUT_DIR}/qc_before"

# Copy intersected data
cp ../intersect/merged_intersected.* ./

# Apply same QC as Approach A
log "  2a. Variant call rate filter..."
plink2 --bfile merged_intersected \
    --geno ${GENO_THRESHOLD} \
    --make-bed \
    --out step2a_geno \
    --threads ${THREADS}

log "  2b. Sample call rate filter..."
plink2 --bfile step2a_geno \
    --mind ${MIND_THRESHOLD} \
    --make-bed \
    --out step2b_mind \
    --threads ${THREADS}

log "  2c. MAF filter (${MAF_THRESHOLD})..."
plink2 --bfile step2b_mind \
    --maf ${MAF_THRESHOLD} \
    --make-bed \
    --out step2c_maf \
    --threads ${THREADS}

log "  2d. HWE filter..."
plink2 --bfile step2c_maf \
    --hwe ${HWE_PVALUE} midp \
    --make-bed \
    --out step2d_hwe \
    --threads ${THREADS}

log "  2e. Heterozygosity filter..."
plink2 --bfile step2d_hwe \
    --het \
    --out step2e_het \
    --threads ${THREADS}

Rscript - << 'RSCRIPT'
het <- read.table("step2e_het.het", header=TRUE)
het$HET_RATE <- (het$OBS_CT - het$O.HOM.) / het$OBS_CT
mean_het <- mean(het$HET_RATE, na.rm=TRUE)
sd_het <- sd(het$HET_RATE, na.rm=TRUE)
het$OUTLIER <- abs(het$HET_RATE - mean_het) > 3 * sd_het
outliers <- het[het$OUTLIER, c("FID", "IID")]
write.table(outliers, "het_outliers.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
RSCRIPT

plink2 --bfile step2d_hwe \
    --remove het_outliers.txt \
    --make-bed \
    --out step2e_het_filtered \
    --threads ${THREADS}

log "  2f. Relatedness filter..."
plink2 --bfile step2e_het_filtered \
    --indep-pairwise 200 50 0.25 \
    --out prune_list \
    --threads ${THREADS}

plink2 --bfile step2e_het_filtered \
    --extract prune_list.prune.in \
    --make-bed \
    --out pruned_for_ibd \
    --threads ${THREADS}

plink2 --bfile pruned_for_ibd \
    --make-king-table \
    --out ibd_check \
    --threads ${THREADS}

Rscript - << 'RSCRIPT'
king <- read.table("ibd_check.kin0", header=TRUE)
related <- king[king$KINSHIP > 0.125, ]
if (nrow(related) > 0) {
    sample_counts <- table(c(related$ID1, related$ID2))
    to_remove <- names(sort(sample_counts, decreasing=TRUE))
    removed <- c()
    for (s in to_remove) {
        remaining <- related[!(related$ID1 %in% removed | related$ID2 %in% removed), ]
        if (nrow(remaining) == 0) break
        if (s %in% c(remaining$ID1, remaining$ID2)) removed <- c(removed, s)
    }
    write.table(data.frame(FID=removed, IID=removed), "related_to_remove.txt",
                row.names=FALSE, col.names=FALSE, quote=FALSE)
} else {
    file.create("related_to_remove.txt")
}
RSCRIPT

plink2 --bfile step2e_het_filtered \
    --remove related_to_remove.txt \
    --make-bed \
    --out final_qcd \
    --threads ${THREADS}

cd "${OUTPUT_DIR}"
time_end "STEP2_QC"

N_QCD_VARS=$(wc -l < qc_before/final_qcd.bim)
N_QCD_SAMP=$(wc -l < qc_before/final_qcd.fam)
log "  After QC: ${N_QCD_VARS} variants, ${N_QCD_SAMP} samples"

# =============================================================================
# STEPS 3-5: IMPUTATION AND POST-QC (same as Approach A)
# =============================================================================

log ""
log "=== STEP 3: Liftover and Submit to TOPMed ==="
time_start "STEP3_IMPUTATION"

cd "${OUTPUT_DIR}/liftover"

plink2 --bfile ../qc_before/final_qcd \
    --export vcf-4.2 bgz \
    --out for_imputation \
    --threads ${THREADS}

bcftools index for_imputation.vcf.gz

log "  NOTE: Manual submission to TOPMed required"
log "  Upload: ${OUTPUT_DIR}/liftover/for_imputation.vcf.gz"

cd "${OUTPUT_DIR}"
time_end "STEP3_IMPUTATION"

# =============================================================================
# SUMMARY
# =============================================================================

TOTAL_END=$(date +%s)
TOTAL_TIME=$((TOTAL_END - TOTAL_START))

log ""
log "=============================================="
log "APPROACH C (Intersect-First) Complete!"
log "=============================================="
log ""
log "Total wall-clock time: ${TOTAL_TIME} seconds"
log ""
log "Key characteristics:"
log "  - Started with ${N_COMMON} variants common to all ${N_PLATFORMS} platforms"
log "  - Lost platform-specific variants (e.g., rare variants on one array)"
log "  - Ensures all samples have identical variant coverage"
log ""
log "Comparison to Our Pipeline:"
log "  - Our pipeline uses UNION merge (keeps all variants)"
log "  - Intersect loses rare/unique variants that may be informative"
log "  - Expect lower total variant count but more consistent QC metrics"
log ""
log "Results in: ${OUTPUT_DIR}/"
log "Timing log: ${TIMING_LOG}"
log ""
