#!/bin/bash
################################################################################
# common_functions.sh
#
# Shared functions for benchmark approach scripts.
# Source this file in approach scripts: source "$(dirname "$0")/common_functions.sh"
#
################################################################################

# =============================================================================
# Logging Functions
# =============================================================================

setup_logging() {
    local output_dir="$1"
    mkdir -p "${output_dir}/logs"
    TIMING_LOG="${output_dir}/logs/timing.log"
    PIPELINE_LOG="${output_dir}/logs/pipeline.log"
    export TIMING_LOG PIPELINE_LOG
}

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${PIPELINE_LOG:-/dev/null}"
}

time_start() {
    echo "$1 START $(date +%s)" >> "${TIMING_LOG:-/dev/null}"
}

time_end() {
    echo "$1 END $(date +%s)" >> "${TIMING_LOG:-/dev/null}"
}

# =============================================================================
# QC Functions
# =============================================================================

# Standard stringent QC (Approach A style)
run_stringent_qc() {
    local input_prefix="$1"
    local output_prefix="$2"
    local threads="${3:-4}"
    local maf="${4:-0.01}"
    local hwe="${5:-1e-6}"
    local geno="${6:-0.02}"
    local mind="${7:-0.02}"

    log "Running stringent pre-imputation QC..."

    # Variant call rate
    plink2 --bfile "${input_prefix}" \
        --geno ${geno} \
        --make-bed \
        --out "${output_prefix}_step1_geno" \
        --threads ${threads}

    # Sample call rate
    plink2 --bfile "${output_prefix}_step1_geno" \
        --mind ${mind} \
        --make-bed \
        --out "${output_prefix}_step2_mind" \
        --threads ${threads}

    # MAF filter
    plink2 --bfile "${output_prefix}_step2_mind" \
        --maf ${maf} \
        --make-bed \
        --out "${output_prefix}_step3_maf" \
        --threads ${threads}

    # HWE filter
    plink2 --bfile "${output_prefix}_step3_maf" \
        --hwe ${hwe} midp \
        --make-bed \
        --out "${output_prefix}_step4_hwe" \
        --threads ${threads}

    # Heterozygosity filter
    plink2 --bfile "${output_prefix}_step4_hwe" \
        --het \
        --out "${output_prefix}_het" \
        --threads ${threads}

    Rscript - << 'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
het <- read.table(paste0(prefix, "_het.het"), header=TRUE)
het$HET_RATE <- (het$OBS_CT - het$O.HOM.) / het$OBS_CT
mean_het <- mean(het$HET_RATE, na.rm=TRUE)
sd_het <- sd(het$HET_RATE, na.rm=TRUE)
outliers <- het[abs(het$HET_RATE - mean_het) > 3 * sd_het, c("FID", "IID")]
write.table(outliers, paste0(prefix, "_het_outliers.txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
RSCRIPT

    plink2 --bfile "${output_prefix}_step4_hwe" \
        --remove "${output_prefix}_het_outliers.txt" \
        --make-bed \
        --out "${output_prefix}_step5_het" \
        --threads ${threads}

    # Relatedness filter
    run_relatedness_filter "${output_prefix}_step5_het" "${output_prefix}_qcd" ${threads}
}

# Minimal QC (Approach B / Southam 2011 style)
run_minimal_qc() {
    local input_prefix="$1"
    local output_prefix="$2"
    local threads="${3:-4}"
    local geno="${4:-0.10}"
    local mind="${5:-0.10}"

    log "Running minimal pre-imputation QC (Southam 2011 style)..."

    # Lenient variant call rate only
    plink2 --bfile "${input_prefix}" \
        --geno ${geno} \
        --make-bed \
        --out "${output_prefix}_step1_geno" \
        --threads ${threads}

    # Lenient sample call rate only
    plink2 --bfile "${output_prefix}_step1_geno" \
        --mind ${mind} \
        --make-bed \
        --out "${output_prefix}_qcd" \
        --threads ${threads}
}

# Relatedness filter using KING
run_relatedness_filter() {
    local input_prefix="$1"
    local output_prefix="$2"
    local threads="${3:-4}"
    local kinship_threshold="${4:-0.125}"

    log "Running relatedness filter..."

    # LD prune first
    plink2 --bfile "${input_prefix}" \
        --indep-pairwise 200 50 0.25 \
        --out "${output_prefix}_prune" \
        --threads ${threads}

    plink2 --bfile "${input_prefix}" \
        --extract "${output_prefix}_prune.prune.in" \
        --make-bed \
        --out "${output_prefix}_pruned" \
        --threads ${threads}

    # Calculate KING kinship
    plink2 --bfile "${output_prefix}_pruned" \
        --make-king-table \
        --out "${output_prefix}_king" \
        --threads ${threads}

    # Identify samples to remove
    Rscript - << RSCRIPT
king <- read.table("${output_prefix}_king.kin0", header=TRUE)
related <- king[king\$KINSHIP > ${kinship_threshold}, ]
if (nrow(related) > 0) {
    sample_counts <- table(c(related\$ID1, related\$ID2))
    to_remove <- names(sort(sample_counts, decreasing=TRUE))
    removed <- c()
    for (s in to_remove) {
        remaining <- related[!(related\$ID1 %in% removed | related\$ID2 %in% removed), ]
        if (nrow(remaining) == 0) break
        if (s %in% c(remaining\$ID1, remaining\$ID2)) removed <- c(removed, s)
    }
    write.table(data.frame(FID=removed, IID=removed), "${output_prefix}_remove.txt",
                row.names=FALSE, col.names=FALSE, quote=FALSE)
} else {
    file.create("${output_prefix}_remove.txt")
}
RSCRIPT

    plink2 --bfile "${input_prefix}" \
        --remove "${output_prefix}_remove.txt" \
        --make-bed \
        --out "${output_prefix}" \
        --threads ${threads}
}

# =============================================================================
# Merge/Intersect Functions
# =============================================================================

# Find common variants across multiple PLINK files
find_common_variants() {
    local output_file="$1"
    shift
    local inputs=("$@")

    log "Finding common variants across ${#inputs[@]} files..."

    # Extract variant IDs from first file
    cut -f2 "${inputs[0]}.bim" | sort > "${output_file}_tmp1.txt"

    # Intersect with remaining files
    for ((i=1; i<${#inputs[@]}; i++)); do
        cut -f2 "${inputs[$i]}.bim" | sort > "${output_file}_tmp2.txt"
        comm -12 "${output_file}_tmp1.txt" "${output_file}_tmp2.txt" > "${output_file}_tmp3.txt"
        mv "${output_file}_tmp3.txt" "${output_file}_tmp1.txt"
    done

    mv "${output_file}_tmp1.txt" "${output_file}"
    rm -f "${output_file}_tmp2.txt"

    n_common=$(wc -l < "${output_file}")
    log "  Found ${n_common} common variants"
}

# Merge multiple PLINK files (intersection)
merge_plink_intersect() {
    local output_prefix="$1"
    local common_variants="$2"
    shift 2
    local inputs=("$@")

    log "Merging ${#inputs[@]} files on common variants..."

    # Extract common variants from each file
    local merge_list="${output_prefix}_merge_list.txt"
    > "${merge_list}"

    for ((i=0; i<${#inputs[@]}; i++)); do
        plink2 --bfile "${inputs[$i]}" \
            --extract "${common_variants}" \
            --make-bed \
            --out "${output_prefix}_input${i}"
        echo "${output_prefix}_input${i}" >> "${merge_list}"
    done

    # Merge all
    first_file=$(head -1 "${merge_list}")
    tail -n +2 "${merge_list}" > "${output_prefix}_to_merge.txt"

    plink --bfile "${first_file}" \
        --merge-list "${output_prefix}_to_merge.txt" \
        --make-bed \
        --out "${output_prefix}"
}

# =============================================================================
# Liftover Functions
# =============================================================================

run_liftover_hg38() {
    local input_vcf="$1"
    local output_vcf="$2"
    local chain_file="${3:-../../resources/chain_files/hg19ToHg38.over.chain.gz}"
    local ref_fasta="${4:-../../resources/references/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa}"

    if [[ -f "$chain_file" ]] && [[ -f "$ref_fasta" ]]; then
        log "Running CrossMap liftover to hg38..."
        CrossMap vcf "$chain_file" "$input_vcf" "$ref_fasta" "${output_vcf%.gz}"
        bcftools sort "${output_vcf%.gz}" -Oz -o "${output_vcf}"
        bcftools index "${output_vcf}"
        rm -f "${output_vcf%.gz}"
    else
        log "WARNING: Chain file or reference FASTA not found, skipping liftover"
        cp "$input_vcf" "$output_vcf"
        bcftools index "$output_vcf"
    fi
}

# =============================================================================
# Imputation Server Submission Functions
# =============================================================================

# Submit to TOPMed via imputationbot
submit_topmed() {
    local input_dir="$1"
    local output_dir="$2"
    local token="$3"
    local password="${4:-}"

    log "Submitting to TOPMed Imputation Server..."

    if command -v imputationbot &> /dev/null; then
        imputationbot impute \
            --files "${input_dir}"/*.vcf.gz \
            --refpanel topmed-r2 \
            --build hg38 \
            --output "${output_dir}" \
            --token "${token}" \
            --password "${password}" \
            --wait
    else
        log "ERROR: imputationbot not installed"
        log "Install with: pip install imputationbot"
        return 1
    fi
}

# Submit to Michigan with 1000G panel via imputationbot
submit_michigan_1kg() {
    local input_dir="$1"
    local output_dir="$2"
    local token="$3"
    local password="${4:-}"

    log "Submitting to Michigan Imputation Server (1000G panel)..."

    if command -v imputationbot &> /dev/null; then
        imputationbot impute \
            --files "${input_dir}"/*.vcf.gz \
            --refpanel 1000g-phase3-v5 \
            --build hg38 \
            --output "${output_dir}" \
            --token "${token}" \
            --password "${password}" \
            --wait
    else
        log "ERROR: imputationbot not installed"
        return 1
    fi
}

# Submit to All of Us via terralab
submit_allofus() {
    local input_dir="$1"
    local output_dir="$2"

    log "Submitting to All of Us Imputation Server..."

    if command -v terralab &> /dev/null; then
        terralab submit imputation \
            --input "${input_dir}" \
            --output "${output_dir}" \
            --wait
    else
        log "ERROR: terralab not installed"
        log "See: https://github.com/DataBiosphere/terra-tools"
        return 1
    fi
}

# =============================================================================
# Post-Imputation Filtering
# =============================================================================

# Traditional R² filter
filter_r2() {
    local input_vcf="$1"
    local output_vcf="$2"
    local threshold="${3:-0.3}"

    log "Applying R² filter (threshold: ${threshold})..."
    bcftools view -i "INFO/R2 >= ${threshold}" "${input_vcf}" -Oz -o "${output_vcf}"
    bcftools index "${output_vcf}"
}

# Post-imputation QC (Approach B/D style - thorough)
run_post_imputation_qc() {
    local input_prefix="$1"
    local output_prefix="$2"
    local threads="${3:-4}"
    local maf="${4:-0.01}"
    local hwe="${5:-1e-6}"
    local geno="${6:-0.02}"
    local mind="${7:-0.02}"

    log "Running thorough post-imputation QC..."

    # Call rate filters
    plink2 --bfile "${input_prefix}" \
        --geno ${geno} \
        --mind ${mind} \
        --make-bed \
        --out "${output_prefix}_step1" \
        --threads ${threads}

    # MAF filter (optional for post-imputation)
    if [[ "${maf}" != "0" ]]; then
        plink2 --bfile "${output_prefix}_step1" \
            --maf ${maf} \
            --make-bed \
            --out "${output_prefix}_step2" \
            --threads ${threads}
    else
        cp "${output_prefix}_step1".* "${output_prefix}_step2".*
    fi

    # HWE filter (cautious for admixed)
    plink2 --bfile "${output_prefix}_step2" \
        --hwe ${hwe} midp \
        --make-bed \
        --out "${output_prefix}_step3" \
        --threads ${threads}

    # Relatedness
    run_relatedness_filter "${output_prefix}_step3" "${output_prefix}" ${threads}
}

# =============================================================================
# Utility Functions
# =============================================================================

count_variants_samples() {
    local prefix="$1"
    local n_vars=$(wc -l < "${prefix}.bim" 2>/dev/null || echo "0")
    local n_samp=$(wc -l < "${prefix}.fam" 2>/dev/null || echo "0")
    echo "${n_vars} variants, ${n_samp} samples"
}

seconds_to_human() {
    local seconds=$1
    local days=$((seconds / 86400))
    local hours=$(((seconds % 86400) / 3600))
    local minutes=$(((seconds % 3600) / 60))
    local secs=$((seconds % 60))

    if [[ $days -gt 0 ]]; then
        echo "${days}d ${hours}h ${minutes}m"
    elif [[ $hours -gt 0 ]]; then
        echo "${hours}h ${minutes}m ${secs}s"
    elif [[ $minutes -gt 0 ]]; then
        echo "${minutes}m ${secs}s"
    else
        echo "${secs}s"
    fi
}
