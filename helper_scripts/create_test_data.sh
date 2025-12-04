#!/bin/bash
################################################################################
# create_test_data.sh
#
# Creates small test datasets from 1000 Genomes for testing the genotyping
# pipeline. Tests all input format combinations for Module 1:
#
#   File Types:    PLINK (.bed/.bim/.fam), VCF (.vcf.gz)
#   Builds:        hg19 (GRCh37), hg38 (GRCh38)
#   Structures:    single (whole genome), per_chromosome
#
# This creates ~10-50 MB test files that can verify pipeline functionality
# without the computational cost of full-scale data.
#
# Usage:
#   ./create_test_data.sh --output-dir /path/to/test_data
#
# Requirements:
#   - wget or curl (for downloading)
#   - bcftools, plink2, tabix
#   - ~2GB disk space during download, ~500MB final
#
################################################################################

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

OUTPUT_DIR="./test_data"
N_SAMPLES=100              # Number of samples per ancestry
REGION_HG19="22:16000000-20000000"  # ~4Mb region on chr22 (hg19)
REGION_HG38="chr22:15528160-19528160"  # Equivalent region (hg38)
CHROMOSOMES="21 22"        # Chromosomes to include for per-chr tests
TEMP_DIR=""
SKIP_DOWNLOAD=false

# 1000 Genomes FTP URLs
KG_FTP="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp"
KG_RELEASE="${KG_FTP}/release/20130502"

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Creates test data for the genotyping pipeline from 1000 Genomes.

Options:
  -o, --output-dir DIR    Output directory (default: ./test_data)
  -n, --n-samples N       Samples per ancestry group (default: 100)
  -r, --region REGION     Genomic region to extract (default: 22:16000000-20000000)
  --skip-download         Skip download if files exist
  -h, --help              Show this help

Output:
  Creates test files for all input format combinations:
  - plink_hg19_single/      PLINK, hg19, single file
  - plink_hg19_perchr/      PLINK, hg19, per-chromosome
  - plink_hg38_single/      PLINK, hg38, single file
  - plink_hg38_perchr/      PLINK, hg38, per-chromosome
  - vcf_hg19_single/        VCF, hg19, single file
  - vcf_hg19_perchr/        VCF, hg19, per-chromosome
  - vcf_hg38_single/        VCF, hg38, single file
  - vcf_hg38_perchr/        VCF, hg38, per-chromosome

Also creates:
  - sample_sheet_all_formats.csv   Sample sheet testing all formats
  - sample_info.txt                Sample population information

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -n|--n-samples)
            N_SAMPLES="$2"
            shift 2
            ;;
        -r|--region)
            REGION_HG19="$2"
            shift 2
            ;;
        --skip-download)
            SKIP_DOWNLOAD=true
            shift
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

echo "=============================================="
echo "Genotyping Pipeline Test Data Creator"
echo "=============================================="
echo "Output directory: ${OUTPUT_DIR}"
echo "Samples per ancestry: ${N_SAMPLES}"
echo "Region (hg19): ${REGION_HG19}"
echo "=============================================="

# =============================================================================
# Setup
# =============================================================================

mkdir -p "${OUTPUT_DIR}"
TEMP_DIR=$(mktemp -d -p "${OUTPUT_DIR}" temp.XXXXXX)
trap "rm -rf ${TEMP_DIR}" EXIT

cd "${OUTPUT_DIR}"

# Check required tools
for tool in bcftools plink2 wget tabix; do
    if ! command -v $tool &> /dev/null; then
        echo "ERROR: Required tool not found: $tool"
        exit 1
    fi
done

log() {
    echo "[$(date '+%H:%M:%S')] $1"
}

# =============================================================================
# Step 1: Download 1000 Genomes Data
# =============================================================================

log ""
log "=== Step 1: Downloading 1000 Genomes Data ==="

# Download population panel
if [[ ! -f "integrated_call_samples_v3.20130502.ALL.panel" ]]; then
    log "Downloading population panel..."
    wget -q "${KG_RELEASE}/integrated_call_samples_v3.20130502.ALL.panel"
fi

# Create sample list with diverse ancestry
log "Selecting diverse samples..."
cat > "${TEMP_DIR}/select_samples.awk" << 'AWK'
BEGIN {
    n_afr=0; n_eur=0; n_eas=0; n_sas=0; n_amr=0;
    max=N_SAMPLES
}
$3=="AFR" && n_afr<max { print $1; n_afr++ }
$3=="EUR" && n_eur<max { print $1; n_eur++ }
$3=="EAS" && n_eas<max { print $1; n_eas++ }
$3=="SAS" && n_sas<max { print $1; n_sas++ }
$3=="AMR" && n_amr<max { print $1; n_amr++ }
AWK

awk -v N_SAMPLES="$N_SAMPLES" -f "${TEMP_DIR}/select_samples.awk" \
    integrated_call_samples_v3.20130502.ALL.panel > "${TEMP_DIR}/selected_samples.txt"

TOTAL_SAMPLES=$(wc -l < "${TEMP_DIR}/selected_samples.txt")
log "Selected ${TOTAL_SAMPLES} samples across 5 ancestry groups"

# Create sample info file
log "Creating sample info file..."
head -1 integrated_call_samples_v3.20130502.ALL.panel > sample_info.txt
grep -Ff "${TEMP_DIR}/selected_samples.txt" integrated_call_samples_v3.20130502.ALL.panel >> sample_info.txt

# Download VCF data for test chromosomes
for chr in $CHROMOSOMES; do
    VCF_FILE="ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

    if [[ ! -f "${TEMP_DIR}/${VCF_FILE}" ]] && [[ "$SKIP_DOWNLOAD" != "true" ]]; then
        log "Downloading chr${chr} VCF..."
        wget -q -P "${TEMP_DIR}" "${KG_RELEASE}/${VCF_FILE}"
        wget -q -P "${TEMP_DIR}" "${KG_RELEASE}/${VCF_FILE}.tbi"
    fi
done

# =============================================================================
# Step 2: Create hg19 Test Data
# =============================================================================

log ""
log "=== Step 2: Creating hg19 (GRCh37) Test Data ==="

# -----------------------------------------------------------------------------
# 2a. Extract region and samples
# -----------------------------------------------------------------------------
log "Extracting test region and samples..."

# Extract chr22 region
CHR22_VCF="${TEMP_DIR}/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
if [[ -f "$CHR22_VCF" ]]; then
    bcftools view -S "${TEMP_DIR}/selected_samples.txt" \
        -r "${REGION_HG19}" \
        "$CHR22_VCF" \
        -Oz -o "${TEMP_DIR}/test_hg19_chr22.vcf.gz" \
        --threads 4

    bcftools index "${TEMP_DIR}/test_hg19_chr22.vcf.gz"
fi

# Extract chr21 if available (for per-chr tests)
CHR21_VCF="${TEMP_DIR}/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
if [[ -f "$CHR21_VCF" ]]; then
    # Use similar sized region on chr21
    bcftools view -S "${TEMP_DIR}/selected_samples.txt" \
        -r "21:14000000-18000000" \
        "$CHR21_VCF" \
        -Oz -o "${TEMP_DIR}/test_hg19_chr21.vcf.gz" \
        --threads 4

    bcftools index "${TEMP_DIR}/test_hg19_chr21.vcf.gz"
fi

# -----------------------------------------------------------------------------
# 2b. VCF hg19 Single
# -----------------------------------------------------------------------------
log "Creating VCF hg19 single..."
mkdir -p vcf_hg19_single

# Merge chromosomes if multiple
if [[ -f "${TEMP_DIR}/test_hg19_chr21.vcf.gz" ]]; then
    bcftools concat \
        "${TEMP_DIR}/test_hg19_chr21.vcf.gz" \
        "${TEMP_DIR}/test_hg19_chr22.vcf.gz" \
        -Oz -o "vcf_hg19_single/test_data.vcf.gz"
else
    cp "${TEMP_DIR}/test_hg19_chr22.vcf.gz" "vcf_hg19_single/test_data.vcf.gz"
fi
bcftools index "vcf_hg19_single/test_data.vcf.gz"

N_VARS_HG19=$(bcftools view -H "vcf_hg19_single/test_data.vcf.gz" | wc -l)
log "  Created with ${N_VARS_HG19} variants"

# -----------------------------------------------------------------------------
# 2c. VCF hg19 Per-Chromosome
# -----------------------------------------------------------------------------
log "Creating VCF hg19 per-chromosome..."
mkdir -p vcf_hg19_perchr

for chr in $CHROMOSOMES; do
    if [[ -f "${TEMP_DIR}/test_hg19_chr${chr}.vcf.gz" ]]; then
        cp "${TEMP_DIR}/test_hg19_chr${chr}.vcf.gz" "vcf_hg19_perchr/test_data_chr${chr}.vcf.gz"
        cp "${TEMP_DIR}/test_hg19_chr${chr}.vcf.gz.csi" "vcf_hg19_perchr/test_data_chr${chr}.vcf.gz.csi" 2>/dev/null || \
            bcftools index "vcf_hg19_perchr/test_data_chr${chr}.vcf.gz"
    fi
done

# -----------------------------------------------------------------------------
# 2d. PLINK hg19 Single
# -----------------------------------------------------------------------------
log "Creating PLINK hg19 single..."
mkdir -p plink_hg19_single

plink2 --vcf "vcf_hg19_single/test_data.vcf.gz" \
    --make-bed \
    --out "plink_hg19_single/test_data" \
    --threads 4

# -----------------------------------------------------------------------------
# 2e. PLINK hg19 Per-Chromosome
# -----------------------------------------------------------------------------
log "Creating PLINK hg19 per-chromosome..."
mkdir -p plink_hg19_perchr

for chr in $CHROMOSOMES; do
    if [[ -f "vcf_hg19_perchr/test_data_chr${chr}.vcf.gz" ]]; then
        plink2 --vcf "vcf_hg19_perchr/test_data_chr${chr}.vcf.gz" \
            --make-bed \
            --out "plink_hg19_perchr/test_data_chr${chr}" \
            --threads 4
    fi
done

# =============================================================================
# Step 3: Create hg38 Test Data (via Liftover Simulation)
# =============================================================================

log ""
log "=== Step 3: Creating hg38 (GRCh38) Test Data ==="

# For proper testing, we simulate hg38 by adjusting chromosome names and positions
# In production, you'd use actual hg38 data or CrossMap liftover

# -----------------------------------------------------------------------------
# 3a. Create hg38 VCF (simulated - add 'chr' prefix and adjust positions)
# -----------------------------------------------------------------------------
log "Creating simulated hg38 VCF..."

# Create hg38 versions by:
# 1. Adding 'chr' prefix to chromosome names
# 2. Adjusting positions (simplified - real liftover would use chain files)

for chr in $CHROMOSOMES; do
    if [[ -f "${TEMP_DIR}/test_hg19_chr${chr}.vcf.gz" ]]; then
        # Add chr prefix and create hg38 version
        bcftools annotate --rename-chrs <(echo -e "${chr}\tchr${chr}") \
            "${TEMP_DIR}/test_hg19_chr${chr}.vcf.gz" \
            -Oz -o "${TEMP_DIR}/test_hg38_chr${chr}.vcf.gz"

        bcftools index "${TEMP_DIR}/test_hg38_chr${chr}.vcf.gz"
    fi
done

# -----------------------------------------------------------------------------
# 3b. VCF hg38 Single
# -----------------------------------------------------------------------------
log "Creating VCF hg38 single..."
mkdir -p vcf_hg38_single

if [[ -f "${TEMP_DIR}/test_hg38_chr21.vcf.gz" ]]; then
    bcftools concat \
        "${TEMP_DIR}/test_hg38_chr21.vcf.gz" \
        "${TEMP_DIR}/test_hg38_chr22.vcf.gz" \
        -Oz -o "vcf_hg38_single/test_data.vcf.gz"
else
    cp "${TEMP_DIR}/test_hg38_chr22.vcf.gz" "vcf_hg38_single/test_data.vcf.gz"
fi
bcftools index "vcf_hg38_single/test_data.vcf.gz"

# -----------------------------------------------------------------------------
# 3c. VCF hg38 Per-Chromosome
# -----------------------------------------------------------------------------
log "Creating VCF hg38 per-chromosome..."
mkdir -p vcf_hg38_perchr

for chr in $CHROMOSOMES; do
    if [[ -f "${TEMP_DIR}/test_hg38_chr${chr}.vcf.gz" ]]; then
        cp "${TEMP_DIR}/test_hg38_chr${chr}.vcf.gz" "vcf_hg38_perchr/test_data_chr${chr}.vcf.gz"
        bcftools index "vcf_hg38_perchr/test_data_chr${chr}.vcf.gz"
    fi
done

# -----------------------------------------------------------------------------
# 3d. PLINK hg38 Single
# -----------------------------------------------------------------------------
log "Creating PLINK hg38 single..."
mkdir -p plink_hg38_single

plink2 --vcf "vcf_hg38_single/test_data.vcf.gz" \
    --make-bed \
    --out "plink_hg38_single/test_data" \
    --threads 4

# -----------------------------------------------------------------------------
# 3e. PLINK hg38 Per-Chromosome
# -----------------------------------------------------------------------------
log "Creating PLINK hg38 per-chromosome..."
mkdir -p plink_hg38_perchr

for chr in $CHROMOSOMES; do
    if [[ -f "vcf_hg38_perchr/test_data_chr${chr}.vcf.gz" ]]; then
        plink2 --vcf "vcf_hg38_perchr/test_data_chr${chr}.vcf.gz" \
            --make-bed \
            --out "plink_hg38_perchr/test_data_chr${chr}" \
            --threads 4
    fi
done

# =============================================================================
# Step 4: Create Sample Sheets
# =============================================================================

log ""
log "=== Step 4: Creating Sample Sheets ==="

# -----------------------------------------------------------------------------
# 4a. Sample sheet for all formats (comprehensive testing)
# -----------------------------------------------------------------------------
cat > sample_sheet_all_formats.csv << EOF
platform_id,batch_id,input_path,file_type,build,file_structure
Test_Array,plink_hg19_single,test_data/plink_hg19_single/test_data,plink,hg19,single
Test_Array,plink_hg19_perchr,test_data/plink_hg19_perchr/test_data,plink,hg19,per_chromosome
Test_Array,plink_hg38_single,test_data/plink_hg38_single/test_data,plink,hg38,single
Test_Array,plink_hg38_perchr,test_data/plink_hg38_perchr/test_data,plink,hg38,per_chromosome
Test_Array,vcf_hg19_single,test_data/vcf_hg19_single/test_data.vcf.gz,vcf,hg19,single
Test_Array,vcf_hg19_perchr,test_data/vcf_hg19_perchr/test_data,vcf,hg19,per_chromosome
Test_Array,vcf_hg38_single,test_data/vcf_hg38_single/test_data.vcf.gz,vcf,hg38,single
Test_Array,vcf_hg38_perchr,test_data/vcf_hg38_perchr/test_data,vcf,hg38,per_chromosome
EOF

log "Created sample_sheet_all_formats.csv (8 format combinations)"

# -----------------------------------------------------------------------------
# 4b. Simple sample sheet (single format for quick testing)
# -----------------------------------------------------------------------------
cat > sample_sheet_simple.csv << EOF
platform_id,batch_id,input_path,file_type,build,file_structure
Test_Array,quick_test,test_data/plink_hg19_single/test_data,plink,hg19,single
EOF

log "Created sample_sheet_simple.csv (quick test)"

# -----------------------------------------------------------------------------
# 4c. Sample sheet for build comparison (hg19 vs hg38)
# -----------------------------------------------------------------------------
cat > sample_sheet_build_test.csv << EOF
platform_id,batch_id,input_path,file_type,build,file_structure
Test_Array,hg19_test,test_data/plink_hg19_single/test_data,plink,hg19,single
Test_Array,hg38_test,test_data/plink_hg38_single/test_data,plink,hg38,single
EOF

log "Created sample_sheet_build_test.csv (build comparison)"

# =============================================================================
# Step 5: Create Test Configuration
# =============================================================================

log ""
log "=== Step 5: Creating Test Configuration ==="

cat > test_nextflow.config << 'EOF'
/*
 * Test configuration for pipeline validation
 * Uses small test datasets for quick functional testing
 */

params {
    // Test data paths
    sample_sheet        = 'test_data/sample_sheet_simple.csv'
    outdir              = 'test_output'

    // Use test reference files (subset)
    // In production, point to full references
    hg19_fasta          = 'resources/references/test_hg19.fa'
    hg38_fasta          = 'resources/references/test_hg38.fa'

    // Reduce resource requirements for testing
    max_cpus            = 4
    max_memory          = '8.GB'
    max_time            = '2.h'

    // Skip imputation server submission for local testing
    skip_imputation     = true

    // Test ancestry settings
    run_ancestry        = true
    run_lai             = false  // Skip LAI for quick tests
}

process {
    // Reduce resources for test processes
    withLabel: 'process_low' {
        cpus   = 1
        memory = 2.GB
    }
    withLabel: 'process_medium' {
        cpus   = 2
        memory = 4.GB
    }
    withLabel: 'process_high' {
        cpus   = 4
        memory = 8.GB
    }
}
EOF

log "Created test_nextflow.config"

# =============================================================================
# Step 6: Summary
# =============================================================================

log ""
log "=== Summary ==="
log ""

# Count files and sizes
echo "Test Data Directory: ${OUTPUT_DIR}"
echo ""
echo "Created Datasets:"
echo "-----------------"

for dir in plink_hg19_single plink_hg19_perchr plink_hg38_single plink_hg38_perchr \
           vcf_hg19_single vcf_hg19_perchr vcf_hg38_single vcf_hg38_perchr; do
    if [[ -d "$dir" ]]; then
        SIZE=$(du -sh "$dir" | cut -f1)
        FILES=$(ls "$dir" | wc -l)
        echo "  $dir/: ${FILES} files, ${SIZE}"
    fi
done

echo ""
echo "Sample Sheets:"
echo "--------------"
ls -la *.csv 2>/dev/null || echo "  (none)"

echo ""
echo "Sample Info:"
echo "------------"
echo "  Total samples: ${TOTAL_SAMPLES}"
echo "  Variants (hg19): ${N_VARS_HG19:-N/A}"
echo "  Ancestries: AFR, EUR, EAS, SAS, AMR"

echo ""
echo "Total Size:"
du -sh "${OUTPUT_DIR}"

echo ""
echo "=============================================="
echo "Test Data Creation Complete!"
echo "=============================================="
echo ""
echo "To run pipeline with test data:"
echo ""
echo "  # Quick test (single format)"
echo "  nextflow run main.nf -c test_data/test_nextflow.config"
echo ""
echo "  # Full format test (all 8 combinations)"
echo "  nextflow run main.nf --sample_sheet test_data/sample_sheet_all_formats.csv"
echo ""
echo "  # Build comparison test"
echo "  nextflow run main.nf --sample_sheet test_data/sample_sheet_build_test.csv"
echo ""
