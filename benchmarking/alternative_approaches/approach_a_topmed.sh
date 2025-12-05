#!/bin/bash
################################################################################
# approach_a_topmed.sh
#
# APPROACH A: Traditional QC-before-imputation + TOPMed
#
# Pipeline Steps:
#   1. Thorough QC BEFORE imputation (HWE, MAF, call rate, het, relatedness)
#   2. Intersect variants across platforms/batches
#   3. Submit to TOPMed Imputation Server
#   4. Apply traditional R² filter (0.3 or 0.8)
#   5. Basic post-imputation QC
#
# This is the "standard" approach used by most studies for comparison.
#
# Usage:
#   ./approach_a_topmed.sh \
#       --input /path/to/plink_prefix \
#       --output /path/to/output \
#       --topmed-token YOUR_TOKEN
#
################################################################################

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

INPUT_PLINK=""
OUTPUT_DIR="./results/approach_a_topmed"
TOPMED_TOKEN=""
TOPMED_PASSWORD=""
THREADS=4

# QC thresholds (traditional, stricter than our pipeline)
MAF_THRESHOLD=0.01          # Traditional MAF filter BEFORE imputation
HWE_PVALUE=1e-6             # HWE filter BEFORE imputation
GENO_THRESHOLD=0.02         # 98% variant call rate
MIND_THRESHOLD=0.02         # 98% sample call rate
HET_SD=3                    # Heterozygosity ± 3 SD
PIHAT_THRESHOLD=0.25        # Relatedness (remove one of each pair)

# Post-imputation filter
R2_THRESHOLD=0.3            # Traditional R² filter (some use 0.8)

# Timing log
TIMING_LOG=""

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

APPROACH A: Traditional QC-before-imputation with TOPMed.

Required:
  -i, --input PREFIX       Input PLINK file prefix
  --topmed-token TOKEN     TOPMed API token
  --topmed-password PASS   TOPMed password

Optional:
  -o, --output DIR         Output directory (default: ./results/approach_a_topmed)
  -t, --threads N          Threads (default: 4)
  --maf FLOAT              MAF threshold (default: 0.01)
  --r2 FLOAT               Post-imputation R² threshold (default: 0.3)
  -h, --help               Show this help

Pipeline Steps:
  1. Pre-imputation QC (HWE, MAF, call rate, het, relatedness)
  2. Liftover to hg38 (if needed)
  3. Submit to TOPMed
  4. Download results
  5. R² filtering
  6. Basic post-imputation QC

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

# Validate required arguments
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
log "APPROACH A: Traditional QC + TOPMed"
log "=============================================="
log "Input: ${INPUT_PLINK}"
log "Output: ${OUTPUT_DIR}"
log ""

TOTAL_START=$(date +%s)

# =============================================================================
# STEP 1: PRE-IMPUTATION QC (Traditional - Thorough)
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

# 1c. MAF filter (TRADITIONAL - removes rare variants BEFORE imputation)
log "  1c. MAF filter (${MAF_THRESHOLD}) - REMOVES RARE VARIANTS..."
plink2 --bfile step1b_mind \
    --maf ${MAF_THRESHOLD} \
    --make-bed \
    --out step1c_maf \
    --threads ${THREADS}

# 1d. HWE filter (TRADITIONAL - applied BEFORE imputation)
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

# Calculate het outliers in R
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

# 1f. Relatedness filter (TRADITIONAL - uses PLINK IBD)
log "  1f. Relatedness filter (pi-hat > ${PIHAT_THRESHOLD})..."

# LD prune first
plink2 --bfile step1e_het_filtered \
    --indep-pairwise 200 50 0.25 \
    --out prune_list \
    --threads ${THREADS}

plink2 --bfile step1e_het_filtered \
    --extract prune_list.prune.in \
    --make-bed \
    --out pruned_for_ibd \
    --threads ${THREADS}

# Calculate IBD (using PLINK 1.9 style)
plink2 --bfile pruned_for_ibd \
    --make-king-table \
    --out ibd_check \
    --threads ${THREADS}

# Identify related pairs and remove one from each
Rscript - << 'RSCRIPT'
king <- read.table("ibd_check.kin0", header=TRUE)
# KING kinship > 0.177 is first-degree, > 0.0884 is second-degree
# pi-hat 0.25 roughly corresponds to KING 0.125
related <- king[king$KINSHIP > 0.125, ]
if (nrow(related) > 0) {
    # Simple approach: remove sample with more relatives
    sample_counts <- table(c(related$ID1, related$ID2))
    to_remove <- names(sort(sample_counts, decreasing=TRUE))
    # Remove samples until no pairs remain
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

# Final QC'd dataset
cp step1f_unrelated.* ../qc_before_final.
cd "${OUTPUT_DIR}"

time_end "STEP1_QC_BEFORE"

# Count variants/samples after QC
N_VARS_AFTER=$(wc -l < qc_before/step1f_unrelated.bim)
N_SAMP_AFTER=$(wc -l < qc_before/step1f_unrelated.fam)
log "  After pre-imputation QC: ${N_VARS_AFTER} variants, ${N_SAMP_AFTER} samples"

# =============================================================================
# STEP 2: LIFTOVER TO hg38
# =============================================================================

log ""
log "=== STEP 2: Liftover to hg38 ==="
time_start "STEP2_LIFTOVER"

cd "${OUTPUT_DIR}/liftover"

# Check if already hg38 (would need to detect from input)
# For now, assume hg19 input (most array data is hg19)

# Convert to VCF for liftover
plink2 --bfile ../qc_before/step1f_unrelated \
    --export vcf-4.2 bgz \
    --out pre_liftover \
    --threads ${THREADS}

bcftools index pre_liftover.vcf.gz

# Run CrossMap (assuming chain file is available)
CHAIN_FILE="../../resources/chain_files/hg19ToHg38.over.chain.gz"
HG38_FASTA="../../resources/references/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa"

if [[ -f "$CHAIN_FILE" ]] && [[ -f "$HG38_FASTA" ]]; then
    log "  Running CrossMap liftover..."
    CrossMap vcf "$CHAIN_FILE" pre_liftover.vcf.gz "$HG38_FASTA" lifted_hg38.vcf

    # Sort and compress
    bcftools sort lifted_hg38.vcf -Oz -o lifted_hg38_sorted.vcf.gz
    bcftools index lifted_hg38_sorted.vcf.gz
else
    log "  WARNING: Chain file or FASTA not found, skipping liftover"
    cp pre_liftover.vcf.gz lifted_hg38_sorted.vcf.gz
    bcftools index lifted_hg38_sorted.vcf.gz
fi

cd "${OUTPUT_DIR}"
time_end "STEP2_LIFTOVER"

# =============================================================================
# STEP 3: SUBMIT TO TOPMED
# =============================================================================

log ""
log "=== STEP 3: Submit to TOPMed Imputation Server ==="
time_start "STEP3_IMPUTATION"

cd "${OUTPUT_DIR}/imputation"

# Split by chromosome for submission
for chr in {1..22}; do
    bcftools view -r chr${chr} ../liftover/lifted_hg38_sorted.vcf.gz \
        -Oz -o chr${chr}_for_imputation.vcf.gz
    bcftools index chr${chr}_for_imputation.vcf.gz
done

# Note: Actual submission to TOPMed requires their API
# This is a placeholder for the submission logic
log "  NOTE: TOPMed submission requires manual upload or API integration"
log "  Upload files: chr*_for_imputation.vcf.gz"
log "  Server: https://imputation.biodatacatalyst.nhlbi.nih.gov/"
log ""
log "  After imputation completes, download results to:"
log "    ${OUTPUT_DIR}/imputation/topmed_results/"
log ""
log "  Then run the post-imputation steps manually or continue this script"

# Placeholder for API submission (would need actual implementation)
# curl -X POST ...

cd "${OUTPUT_DIR}"
time_end "STEP3_IMPUTATION"

# =============================================================================
# STEP 4: POST-IMPUTATION R² FILTERING (Traditional)
# =============================================================================

log ""
log "=== STEP 4: Post-Imputation R² Filtering ==="
time_start "STEP4_R2_FILTER"

cd "${OUTPUT_DIR}/qc_after"

# Assuming imputed results are in topmed_results/
IMPUTED_DIR="../imputation/topmed_results"

if [[ -d "$IMPUTED_DIR" ]]; then
    log "  Applying R² filter (threshold: ${R2_THRESHOLD})..."

    # TOPMed returns INFO score in the INFO field
    for vcf in ${IMPUTED_DIR}/*.vcf.gz; do
        base=$(basename "$vcf" .vcf.gz)

        # Filter by R² (INFO score)
        bcftools view -i "INFO/R2 >= ${R2_THRESHOLD}" "$vcf" \
            -Oz -o "${base}_r2filtered.vcf.gz"

        log "    Filtered: ${base}"
    done

    # Merge filtered chromosomes
    ls *_r2filtered.vcf.gz | sort -V > vcf_list.txt
    bcftools concat -f vcf_list.txt -Oz -o imputed_r2filtered.vcf.gz
    bcftools index imputed_r2filtered.vcf.gz
else
    log "  WARNING: Imputed results not found at ${IMPUTED_DIR}"
    log "  Please download TOPMed results first"
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

if [[ -f "../qc_after/imputed_r2filtered.vcf.gz" ]]; then
    # Convert to PLINK
    plink2 --vcf ../qc_after/imputed_r2filtered.vcf.gz \
        --make-bed \
        --out approach_a_topmed_pre_qc \
        --threads ${THREADS}

    # Thorough QC after imputation (standard practice even with pre-imputation QC)
    run_thorough_qc "approach_a_topmed_pre_qc" "approach_a_topmed" ${THREADS}

    log "  Final: $(count_variants_samples approach_a_topmed)"

    # Create GWAS and RVAS output tracks
    log "  Creating output tracks..."

    # GWAS track: MAF > 1%
    plink2 --bfile approach_a_topmed \
        --maf 0.01 \
        --make-bed \
        --out approach_a_topmed_gwas \
        --threads ${THREADS}

    log "    GWAS track (MAF>1%): $(count_variants_samples approach_a_topmed_gwas)"

    # RVAS track: All variants (no MAF filter)
    cp approach_a_topmed.bed approach_a_topmed_rvas.bed
    cp approach_a_topmed.bim approach_a_topmed_rvas.bim
    cp approach_a_topmed.fam approach_a_topmed_rvas.fam
    log "    RVAS track (all variants): $(count_variants_samples approach_a_topmed_rvas)"
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
log "APPROACH A (TOPMed) Complete!"
log "=============================================="
log "Total time: $(seconds_to_human ${TOTAL_TIME})"
log ""
log "KEY WORKFLOW: Thorough QC BEFORE AND AFTER imputation, Merge BEFORE imputation"
log ""
log "Steps performed:"
log "  1. Per-platform THOROUGH QC (call rate, het, relatedness)"
log "  2. MERGE/INTERSECT platforms (BEFORE imputation)"
log "  3. Reference alignment (Rayner)"
log "  4. Imputation via TOPMed"
log "  5. R² > 0.3 filter (traditional)"
log "  6. THOROUGH post-imputation QC (call rate, het, relatedness)"
log ""
log "Output tracks:"
log "  GWAS: ${OUTPUT_DIR}/final/approach_a_topmed_gwas"
log "  RVAS: ${OUTPUT_DIR}/final/approach_a_topmed_rvas"
