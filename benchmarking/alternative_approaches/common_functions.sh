#!/bin/bash
################################################################################
# common_functions.sh
#
# Shared functions for benchmark approach scripts.
# Source this file in approach scripts: source "$(dirname "$0")/common_functions.sh"
#
# STANDARDIZED SETTINGS FOR FAIR COMPARISON:
# - Call rate: 95% (geno=0.05, mind=0.05) - standard
# - HWE: SKIPPED by default (problematic for mixed cohorts)
# - MAF: 0.01 (1%) when applied
# - Relatedness: KING kinship > 0.125 (same for all)
# - R² filter: 0.3 for traditional approaches (A-D)
# - MagicalRsq-X: For our pipeline (E-F) only
# - Reference alignment: Rayner script for ALL approaches
#
################################################################################

# =============================================================================
# STANDARDIZED THRESHOLDS
# =============================================================================

# Call rate - 95% standard (geno/mind = 0.05)
STANDARD_GENO=0.05
STANDARD_MIND=0.05

# MAF - 1% when filtering
STANDARD_MAF=0.01

# HWE - SKIPPED by default for mixed cohorts
SKIP_HWE=true
HWE_PVALUE=1e-6  # Only used if SKIP_HWE=false

# Heterozygosity
HET_SD_THRESHOLD=3

# Relatedness - same for ALL approaches
KINSHIP_THRESHOLD=0.125  # ~2nd degree

# R² - traditional filter threshold
R2_THRESHOLD=0.3

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
# Reference Alignment (Rayner Script)
# =============================================================================

# Run Rayner-style reference alignment before imputation
# This is REQUIRED for ALL approaches for fair comparison
run_reference_alignment() {
    local input_prefix="$1"
    local output_prefix="$2"
    local ref_panel="${3:-HRC}"  # HRC, 1000G, or TOPMed
    local threads="${4:-4}"

    local script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    local helper_dir="${script_dir}/../../helper_scripts"

    log "Running reference alignment (Rayner-style)..."
    log "  Reference panel: ${ref_panel}"

    # Check for Rayner's tool or our equivalent
    if [[ -f "${helper_dir}/check_strand_topmed.pl" ]]; then
        log "  Using check_strand_topmed.pl..."

        # Create frequency file
        plink2 --bfile "${input_prefix}" \
            --freq \
            --out "${output_prefix}_freq" \
            --threads ${threads}

        # Run strand check
        perl "${helper_dir}/check_strand_topmed.pl" \
            -b "${input_prefix}.bim" \
            -f "${output_prefix}_freq.afreq" \
            -r "${helper_dir}/../resources/ref_panels/${ref_panel}_sites.txt.gz" \
            -h \
            -o "${output_prefix}_strand"

        # Apply fixes
        if [[ -f "${output_prefix}_strand-Exclude.txt" ]]; then
            plink2 --bfile "${input_prefix}" \
                --exclude "${output_prefix}_strand-Exclude.txt" \
                --make-bed \
                --out "${output_prefix}_step1" \
                --threads ${threads}
        else
            cp "${input_prefix}".* "${output_prefix}_step1".*
        fi

        # Flip strands if needed
        if [[ -f "${output_prefix}_strand-Flip.txt" ]] && [[ -s "${output_prefix}_strand-Flip.txt" ]]; then
            plink2 --bfile "${output_prefix}_step1" \
                --flip "${output_prefix}_strand-Flip.txt" \
                --make-bed \
                --out "${output_prefix}_aligned" \
                --threads ${threads}
        else
            mv "${output_prefix}_step1".* "${output_prefix}_aligned".*
        fi

    else
        log "  WARNING: Rayner script not found, using basic alignment..."
        # Fallback: basic checks
        plink2 --bfile "${input_prefix}" \
            --snps-only just-acgt \
            --rm-dup exclude-all \
            --make-bed \
            --out "${output_prefix}_aligned" \
            --threads ${threads}
    fi

    log "  Reference alignment complete"
}

# =============================================================================
# QC Functions - STANDARDIZED
# =============================================================================

# Thorough QC - CONSISTENT across all approaches
# Used BEFORE imputation for A, C
# Used AFTER imputation for B, D
run_thorough_qc() {
    local input_prefix="$1"
    local output_prefix="$2"
    local threads="${3:-4}"
    local apply_maf="${4:-true}"
    local skip_hwe="${5:-true}"  # Default: skip HWE for mixed cohorts

    log "Running thorough QC (standardized)..."
    log "  Call rate: ${STANDARD_GENO}/${STANDARD_MIND}"
    log "  MAF filter: ${apply_maf} (threshold: ${STANDARD_MAF})"
    log "  HWE filter: $([[ $skip_hwe == true ]] && echo 'SKIPPED' || echo 'APPLIED')"

    # Step 1: Variant call rate (95%)
    log "  Step 1: Variant call rate filter (>${STANDARD_GENO})..."
    plink2 --bfile "${input_prefix}" \
        --geno ${STANDARD_GENO} \
        --make-bed \
        --out "${output_prefix}_step1_geno" \
        --threads ${threads}

    # Step 2: Sample call rate (95%)
    log "  Step 2: Sample call rate filter (>${STANDARD_MIND})..."
    plink2 --bfile "${output_prefix}_step1_geno" \
        --mind ${STANDARD_MIND} \
        --make-bed \
        --out "${output_prefix}_step2_mind" \
        --threads ${threads}

    # Step 3: MAF filter (optional)
    if [[ "${apply_maf}" == "true" ]]; then
        log "  Step 3: MAF filter (>${STANDARD_MAF})..."
        plink2 --bfile "${output_prefix}_step2_mind" \
            --maf ${STANDARD_MAF} \
            --make-bed \
            --out "${output_prefix}_step3_maf" \
            --threads ${threads}
    else
        log "  Step 3: MAF filter SKIPPED"
        cp "${output_prefix}_step2_mind".{bed,bim,fam} ./ 2>/dev/null || true
        for ext in bed bim fam; do
            cp "${output_prefix}_step2_mind.${ext}" "${output_prefix}_step3_maf.${ext}" 2>/dev/null || true
        done
    fi

    # Step 4: HWE filter (SKIPPED by default for mixed cohorts)
    if [[ "${skip_hwe}" == "false" ]]; then
        log "  Step 4: HWE filter (p < ${HWE_PVALUE})..."
        plink2 --bfile "${output_prefix}_step3_maf" \
            --hwe ${HWE_PVALUE} midp \
            --make-bed \
            --out "${output_prefix}_step4_hwe" \
            --threads ${threads}
    else
        log "  Step 4: HWE filter SKIPPED (mixed cohort)"
        for ext in bed bim fam; do
            cp "${output_prefix}_step3_maf.${ext}" "${output_prefix}_step4_hwe.${ext}" 2>/dev/null || true
        done
    fi

    # Step 5: Heterozygosity filter
    log "  Step 5: Heterozygosity filter (±${HET_SD_THRESHOLD} SD)..."
    plink2 --bfile "${output_prefix}_step4_hwe" \
        --het \
        --out "${output_prefix}_het" \
        --threads ${threads}

    # Calculate het outliers
    Rscript --vanilla - "${output_prefix}" << 'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
het <- read.table(paste0(prefix, "_het.het"), header=TRUE)
het$HET_RATE <- (het$OBS_CT - het$O.HOM.) / het$OBS_CT
mean_het <- mean(het$HET_RATE, na.rm=TRUE)
sd_het <- sd(het$HET_RATE, na.rm=TRUE)
outliers <- het[abs(het$HET_RATE - mean_het) > 3 * sd_het, c("FID", "IID")]
write.table(outliers, paste0(prefix, "_het_outliers.txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
cat("Heterozygosity outliers:", nrow(outliers), "\n")
RSCRIPT

    if [[ -s "${output_prefix}_het_outliers.txt" ]]; then
        plink2 --bfile "${output_prefix}_step4_hwe" \
            --remove "${output_prefix}_het_outliers.txt" \
            --make-bed \
            --out "${output_prefix}_step5_het" \
            --threads ${threads}
    else
        for ext in bed bim fam; do
            cp "${output_prefix}_step4_hwe.${ext}" "${output_prefix}_step5_het.${ext}" 2>/dev/null || true
        done
    fi

    # Step 6: Relatedness filter (SAME for all approaches)
    log "  Step 6: Relatedness filter (KING > ${KINSHIP_THRESHOLD})..."
    run_relatedness_filter "${output_prefix}_step5_het" "${output_prefix}_qcd" ${threads}

    log "  Thorough QC complete: $(count_variants_samples ${output_prefix}_qcd)"
}

# Minimal QC - Only call rate (for pre-imputation in B, D approaches)
run_minimal_qc() {
    local input_prefix="$1"
    local output_prefix="$2"
    local threads="${3:-4}"

    log "Running minimal QC (call rate only)..."
    log "  Thresholds: geno=${STANDARD_GENO}, mind=${STANDARD_MIND}"

    plink2 --bfile "${input_prefix}" \
        --geno ${STANDARD_GENO} \
        --mind ${STANDARD_MIND} \
        --make-bed \
        --out "${output_prefix}_qcd" \
        --threads ${threads}

    log "  Minimal QC complete: $(count_variants_samples ${output_prefix}_qcd)"
}

# Relatedness filter - STANDARDIZED for all approaches
run_relatedness_filter() {
    local input_prefix="$1"
    local output_prefix="$2"
    local threads="${3:-4}"

    log "  Running relatedness filter (KING > ${KINSHIP_THRESHOLD})..."

    # LD prune first
    plink2 --bfile "${input_prefix}" \
        --indep-pairwise 200 50 0.25 \
        --out "${output_prefix}_prune" \
        --threads ${threads} 2>/dev/null

    plink2 --bfile "${input_prefix}" \
        --extract "${output_prefix}_prune.prune.in" \
        --make-bed \
        --out "${output_prefix}_pruned" \
        --threads ${threads} 2>/dev/null

    # Calculate KING kinship
    plink2 --bfile "${output_prefix}_pruned" \
        --make-king-table \
        --out "${output_prefix}_king" \
        --threads ${threads} 2>/dev/null

    # Identify samples to remove
    if [[ -f "${output_prefix}_king.kin0" ]]; then
        Rscript --vanilla - "${output_prefix}" "${KINSHIP_THRESHOLD}" << 'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
threshold <- as.numeric(args[2])
king <- read.table(paste0(prefix, "_king.kin0"), header=TRUE)
related <- king[king$KINSHIP > threshold, ]
if (nrow(related) > 0) {
    sample_counts <- table(c(related$ID1, related$ID2))
    to_remove <- names(sort(sample_counts, decreasing=TRUE))
    removed <- c()
    for (s in to_remove) {
        remaining <- related[!(related$ID1 %in% removed | related$ID2 %in% removed), ]
        if (nrow(remaining) == 0) break
        if (s %in% c(remaining$ID1, remaining$ID2)) removed <- c(removed, s)
    }
    write.table(data.frame(FID=removed, IID=removed), paste0(prefix, "_remove.txt"),
                row.names=FALSE, col.names=FALSE, quote=FALSE)
    cat("Related samples to remove:", length(removed), "\n")
} else {
    file.create(paste0(prefix, "_remove.txt"))
    cat("No related samples found\n")
}
RSCRIPT
    else
        touch "${output_prefix}_remove.txt"
    fi

    if [[ -s "${output_prefix}_remove.txt" ]]; then
        plink2 --bfile "${input_prefix}" \
            --remove "${output_prefix}_remove.txt" \
            --make-bed \
            --out "${output_prefix}" \
            --threads ${threads}
    else
        for ext in bed bim fam; do
            cp "${input_prefix}.${ext}" "${output_prefix}.${ext}" 2>/dev/null || true
        done
    fi
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

    local n_common=$(wc -l < "${output_file}")
    log "  Found ${n_common} common variants"
}

# Intersect merge across platforms
merge_plink_intersect() {
    local output_prefix="$1"
    shift
    local inputs=("$@")
    local threads="${THREADS:-4}"

    log "Intersect-merging ${#inputs[@]} files..."

    # Find common variants
    find_common_variants "${output_prefix}_common_vars.txt" "${inputs[@]}"

    # Extract common variants from each file
    for ((i=0; i<${#inputs[@]}; i++)); do
        plink2 --bfile "${inputs[$i]}" \
            --extract "${output_prefix}_common_vars.txt" \
            --make-bed \
            --out "${output_prefix}_input${i}" \
            --threads ${threads}
    done

    # Merge
    if [[ ${#inputs[@]} -eq 1 ]]; then
        for ext in bed bim fam; do
            cp "${output_prefix}_input0.${ext}" "${output_prefix}.${ext}"
        done
    else
        # Create merge list
        for ((i=1; i<${#inputs[@]}; i++)); do
            echo "${output_prefix}_input${i}"
        done > "${output_prefix}_merge_list.txt"

        plink2 --bfile "${output_prefix}_input0" \
            --pmerge-list "${output_prefix}_merge_list.txt" bfile \
            --make-bed \
            --out "${output_prefix}" \
            --threads ${threads}
    fi

    log "  Merged: $(count_variants_samples ${output_prefix})"
}

# Union merge (for our pipeline - within platform)
merge_plink_union() {
    local output_prefix="$1"
    shift
    local inputs=("$@")
    local threads="${THREADS:-4}"

    log "Union-merging ${#inputs[@]} files..."

    if [[ ${#inputs[@]} -eq 1 ]]; then
        for ext in bed bim fam; do
            cp "${inputs[0]}.${ext}" "${output_prefix}.${ext}"
        done
    else
        # Create merge list
        for ((i=1; i<${#inputs[@]}; i++)); do
            echo "${inputs[$i]}"
        done > "${output_prefix}_merge_list.txt"

        plink2 --bfile "${inputs[0]}" \
            --pmerge-list "${output_prefix}_merge_list.txt" bfile \
            --make-bed \
            --out "${output_prefix}" \
            --threads ${threads}
    fi

    log "  Union merged: $(count_variants_samples ${output_prefix})"
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

submit_topmed() {
    local input_dir="$1"
    local output_dir="$2"
    local token="${3:-$TOPMED_TOKEN}"
    local password="${4:-$TOPMED_PASSWORD}"

    log "Submitting to TOPMed Imputation Server..."

    if command -v imputationbot &> /dev/null && [[ -n "$token" ]]; then
        imputationbot impute \
            --files "${input_dir}"/*.vcf.gz \
            --refpanel topmed-r2 \
            --build hg38 \
            --output "${output_dir}" \
            --token "${token}" \
            ${password:+--password "${password}"} \
            --wait
    else
        log "NOTE: Manual submission required"
        log "  Upload: ${input_dir}/*.vcf.gz"
        log "  Server: https://imputation.biodatacatalyst.nhlbi.nih.gov/"
        return 1
    fi
}

submit_michigan_1kg() {
    local input_dir="$1"
    local output_dir="$2"
    local token="${3:-$MICHIGAN_TOKEN}"
    local password="${4:-$MICHIGAN_PASSWORD}"

    log "Submitting to Michigan Imputation Server (1000G panel)..."

    if command -v imputationbot &> /dev/null && [[ -n "$token" ]]; then
        imputationbot impute \
            --files "${input_dir}"/*.vcf.gz \
            --refpanel 1000g-phase3-v5 \
            --build hg38 \
            --output "${output_dir}" \
            --token "${token}" \
            ${password:+--password "${password}"} \
            --wait
    else
        log "NOTE: Manual submission required"
        log "  Upload: ${input_dir}/*.vcf.gz"
        log "  Server: https://imputationserver.sph.umich.edu/"
        return 1
    fi
}

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
        log "NOTE: Manual submission required via terralab"
        log "  Install: pip install terralab-cli"
        return 1
    fi
}

# =============================================================================
# Post-Imputation Filtering
# =============================================================================

# Traditional R² filter (for approaches A-D)
filter_r2() {
    local input_vcf="$1"
    local output_vcf="$2"
    local threshold="${3:-$R2_THRESHOLD}"

    log "Applying traditional R² filter (threshold: ${threshold})..."
    bcftools view -i "INFO/R2 >= ${threshold}" "${input_vcf}" -Oz -o "${output_vcf}"
    bcftools index "${output_vcf}"
}

# MagicalRsq-X filter (for our pipeline E-F only)
filter_magicalrsq() {
    local input_vcf="$1"
    local output_vcf="$2"
    local ancestry="${3:-mixed}"
    local threshold="${4:-0.3}"

    local script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    local helper_dir="${script_dir}/../../helper_scripts"

    log "Applying MagicalRsq-X filter (ancestry-calibrated)..."
    log "  Ancestry: ${ancestry}, threshold: ${threshold}"

    if [[ -f "${helper_dir}/magicalrsq_filter.R" ]]; then
        Rscript "${helper_dir}/magicalrsq_filter.R" \
            --vcf "${input_vcf}" \
            --ancestry "${ancestry}" \
            --threshold "${threshold}" \
            --output "${output_vcf}"
    else
        log "WARNING: MagicalRsq-X script not found, falling back to R² filter"
        filter_r2 "${input_vcf}" "${output_vcf}" "${threshold}"
    fi
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

# Prepare VCFs for imputation (split by chromosome)
prepare_vcfs_for_imputation() {
    local input_prefix="$1"
    local output_dir="$2"
    local threads="${3:-4}"

    log "Preparing VCFs for imputation..."
    mkdir -p "${output_dir}"

    # Convert to VCF
    plink2 --bfile "${input_prefix}" \
        --export vcf-4.2 bgz \
        --out "${output_dir}/all_chrs" \
        --threads ${threads}

    bcftools index "${output_dir}/all_chrs.vcf.gz"

    # Split by chromosome
    for chr in {1..22}; do
        bcftools view -r chr${chr} "${output_dir}/all_chrs.vcf.gz" \
            -Oz -o "${output_dir}/chr${chr}.vcf.gz"
        bcftools index "${output_dir}/chr${chr}.vcf.gz"
    done

    log "  VCFs prepared in ${output_dir}/"
}
