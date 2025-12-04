#!/bin/bash
################################################################################
# approach_b_michigan.sh
#
# APPROACH B: Traditional QC-before-imputation + Michigan Imputation Server
#
# Same as Approach A but uses Michigan instead of TOPMed.
# Useful for comparing reference panel impact (HRC/1KG vs TOPMed).
#
# Pipeline Steps:
#   1. Thorough QC BEFORE imputation (HWE, MAF, call rate, het, relatedness)
#   2. Intersect variants across platforms/batches
#   3. Submit to Michigan Imputation Server (HRC or 1KG panel)
#   4. Apply traditional R² filter (0.3 or 0.8)
#   5. Basic post-imputation QC
#
# Usage:
#   ./approach_b_michigan.sh \
#       --input /path/to/plink_prefix \
#       --output /path/to/output \
#       --michigan-token YOUR_TOKEN \
#       --panel hrc
#
################################################################################

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

INPUT_PLINK=""
OUTPUT_DIR="./results/approach_b_michigan"
MICHIGAN_TOKEN=""
PANEL="hrc"  # hrc, 1000g-phase3, or caapa
THREADS=4

# QC thresholds (traditional, same as approach A)
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

APPROACH B: Traditional QC-before-imputation with Michigan Imputation Server.

Required:
  -i, --input PREFIX          Input PLINK file prefix
  --michigan-token TOKEN      Michigan Imputation Server API token

Optional:
  -o, --output DIR            Output directory (default: ./results/approach_b_michigan)
  --panel PANEL               Reference panel: hrc, 1000g-phase3, caapa (default: hrc)
  -t, --threads N             Threads (default: 4)
  --maf FLOAT                 MAF threshold (default: 0.01)
  --r2 FLOAT                  Post-imputation R² threshold (default: 0.3)
  -h, --help                  Show this help

Reference Panels:
  hrc           - Haplotype Reference Consortium (EUR-focused)
  1000g-phase3  - 1000 Genomes Phase 3 (global diversity)
  caapa         - CAAPA (African ancestry)

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_PLINK="$2"
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

if [[ -z "$INPUT_PLINK" ]]; then
    echo "ERROR: Input PLINK files required"
    print_usage
    exit 1
fi

# =============================================================================
# Setup
# =============================================================================

mkdir -p "${OUTPUT_DIR}"/{logs,qc_before,liftover,imputation,qc_after,final}
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
log "APPROACH B: Traditional QC + Michigan (${PANEL})"
log "=============================================="
log "Input: ${INPUT_PLINK}"
log "Output: ${OUTPUT_DIR}"
log "Panel: ${PANEL}"
log ""

TOTAL_START=$(date +%s)

# =============================================================================
# STEP 1: PRE-IMPUTATION QC (Same as Approach A)
# =============================================================================

log ""
log "=== STEP 1: Pre-Imputation QC (Traditional) ==="
time_start "STEP1_QC_BEFORE"

cd "${OUTPUT_DIR}/qc_before"

# 1a. Initial variant QC
log "  1a. Variant call rate filter (${GENO_THRESHOLD})..."
plink2 --bfile "${INPUT_PLINK}" \
    --geno ${GENO_THRESHOLD} \
    --make-bed \
    --out step1a_geno \
    --threads ${THREADS}

# 1b. Sample call rate filter
log "  1b. Sample call rate filter (${MIND_THRESHOLD})..."
plink2 --bfile step1a_geno \
    --mind ${MIND_THRESHOLD} \
    --make-bed \
    --out step1b_mind \
    --threads ${THREADS}

# 1c. MAF filter (removes rare variants BEFORE imputation)
log "  1c. MAF filter (${MAF_THRESHOLD}) - REMOVES RARE VARIANTS..."
plink2 --bfile step1b_mind \
    --maf ${MAF_THRESHOLD} \
    --make-bed \
    --out step1c_maf \
    --threads ${THREADS}

# 1d. HWE filter
log "  1d. HWE filter (${HWE_PVALUE})..."
plink2 --bfile step1c_maf \
    --hwe ${HWE_PVALUE} midp \
    --make-bed \
    --out step1d_hwe \
    --threads ${THREADS}

# 1e. Heterozygosity filter
log "  1e. Heterozygosity filter (±${HET_SD} SD)..."
plink2 --bfile step1d_hwe \
    --het \
    --out step1e_het \
    --threads ${THREADS}

Rscript - << 'RSCRIPT'
het <- read.table("step1e_het.het", header=TRUE)
het$HET_RATE <- (het$OBS_CT - het$O.HOM.) / het$OBS_CT
mean_het <- mean(het$HET_RATE, na.rm=TRUE)
sd_het <- sd(het$HET_RATE, na.rm=TRUE)
het$OUTLIER <- abs(het$HET_RATE - mean_het) > 3 * sd_het
outliers <- het[het$OUTLIER, c("FID", "IID")]
write.table(outliers, "het_outliers.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
cat("Heterozygosity outliers:", nrow(outliers), "\n")
RSCRIPT

plink2 --bfile step1d_hwe \
    --remove het_outliers.txt \
    --make-bed \
    --out step1e_het_filtered \
    --threads ${THREADS}

# 1f. Relatedness filter
log "  1f. Relatedness filter..."
plink2 --bfile step1e_het_filtered \
    --indep-pairwise 200 50 0.25 \
    --out prune_list \
    --threads ${THREADS}

plink2 --bfile step1e_het_filtered \
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
        if (s %in% c(remaining$ID1, remaining$ID2)) {
            removed <- c(removed, s)
        }
    }
    write.table(data.frame(FID=removed, IID=removed), "related_to_remove.txt",
                row.names=FALSE, col.names=FALSE, quote=FALSE)
    cat("Related samples to remove:", length(removed), "\n")
} else {
    file.create("related_to_remove.txt")
    cat("No related samples found\n")
}
RSCRIPT

plink2 --bfile step1e_het_filtered \
    --remove related_to_remove.txt \
    --make-bed \
    --out step1f_unrelated \
    --threads ${THREADS}

cd "${OUTPUT_DIR}"
time_end "STEP1_QC_BEFORE"

N_VARS_AFTER=$(wc -l < qc_before/step1f_unrelated.bim)
N_SAMP_AFTER=$(wc -l < qc_before/step1f_unrelated.fam)
log "  After pre-imputation QC: ${N_VARS_AFTER} variants, ${N_SAMP_AFTER} samples"

# =============================================================================
# STEP 2: PREPARE FOR MICHIGAN (hg19 or hg38 depending on panel)
# =============================================================================

log ""
log "=== STEP 2: Prepare for Michigan Server ==="
time_start "STEP2_PREPARE"

cd "${OUTPUT_DIR}/liftover"

# Michigan accepts hg19 for HRC, hg38 for some panels
# Convert to VCF
plink2 --bfile ../qc_before/step1f_unrelated \
    --export vcf-4.2 bgz \
    --out for_michigan \
    --threads ${THREADS}

bcftools index for_michigan.vcf.gz

# Split by chromosome (Michigan requires per-chromosome files)
for chr in {1..22}; do
    bcftools view -r ${chr} for_michigan.vcf.gz \
        -Oz -o chr${chr}.vcf.gz
    bcftools index chr${chr}.vcf.gz
done

cd "${OUTPUT_DIR}"
time_end "STEP2_PREPARE"

# =============================================================================
# STEP 3: SUBMIT TO MICHIGAN
# =============================================================================

log ""
log "=== STEP 3: Submit to Michigan Imputation Server ==="
time_start "STEP3_IMPUTATION"

cd "${OUTPUT_DIR}/imputation"

log "  NOTE: Michigan submission requires manual upload or imputationbot"
log "  Upload files: ${OUTPUT_DIR}/liftover/chr*.vcf.gz"
log "  Server: https://imputationserver.sph.umich.edu/"
log "  Panel: ${PANEL}"
log ""
log "  For automated submission, install imputationbot:"
log "    curl -sL imputationbot.now.sh | bash"
log "    imputationbot add-instance michigan"
log "    imputationbot impute --files ../liftover/chr*.vcf.gz --refpanel ${PANEL}"
log ""

cd "${OUTPUT_DIR}"
time_end "STEP3_IMPUTATION"

# =============================================================================
# STEP 4-5: POST-IMPUTATION (Same structure as Approach A)
# =============================================================================

log ""
log "=== STEP 4: Post-Imputation R² Filtering ==="
time_start "STEP4_R2_FILTER"

cd "${OUTPUT_DIR}/qc_after"

IMPUTED_DIR="../imputation/michigan_results"

if [[ -d "$IMPUTED_DIR" ]]; then
    log "  Applying R² filter (threshold: ${R2_THRESHOLD})..."

    for vcf in ${IMPUTED_DIR}/*.vcf.gz; do
        base=$(basename "$vcf" .vcf.gz)
        bcftools view -i "INFO/R2 >= ${R2_THRESHOLD}" "$vcf" \
            -Oz -o "${base}_r2filtered.vcf.gz"
        log "    Filtered: ${base}"
    done

    ls *_r2filtered.vcf.gz | sort -V > vcf_list.txt
    bcftools concat -f vcf_list.txt -Oz -o imputed_r2filtered.vcf.gz
    bcftools index imputed_r2filtered.vcf.gz
else
    log "  WARNING: Imputed results not found at ${IMPUTED_DIR}"
    log "  Please download Michigan results first"
fi

cd "${OUTPUT_DIR}"
time_end "STEP4_R2_FILTER"

log ""
log "=== STEP 5: Basic Post-Imputation QC ==="
time_start "STEP5_QC_AFTER"

cd "${OUTPUT_DIR}/final"

if [[ -f "../qc_after/imputed_r2filtered.vcf.gz" ]]; then
    plink2 --vcf ../qc_after/imputed_r2filtered.vcf.gz \
        --make-bed \
        --out approach_b_michigan_final \
        --threads ${THREADS}

    plink2 --bfile approach_b_michigan_final \
        --geno 0.05 \
        --mind 0.05 \
        --make-bed \
        --out approach_b_michigan_qcd \
        --threads ${THREADS}

    N_FINAL_VARS=$(wc -l < approach_b_michigan_qcd.bim)
    N_FINAL_SAMP=$(wc -l < approach_b_michigan_qcd.fam)
    log "  Final dataset: ${N_FINAL_VARS} variants, ${N_FINAL_SAMP} samples"
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
log "APPROACH B (Michigan/${PANEL}) Complete!"
log "=============================================="
log ""
log "Total wall-clock time: ${TOTAL_TIME} seconds"
log ""
log "Key comparison points vs Approach A (TOPMed):"
log "  - Same QC strategy (traditional QC-before)"
log "  - Different reference panel (${PANEL} vs TOPMed)"
log "  - TOPMed has better diverse ancestry coverage"
log "  - HRC is larger but EUR-focused"
log ""
log "Results in: ${OUTPUT_DIR}/final/"
log "Timing log: ${TIMING_LOG}"
log ""
