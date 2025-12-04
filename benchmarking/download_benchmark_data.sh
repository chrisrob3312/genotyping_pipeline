#!/bin/bash
################################################################################
# download_benchmark_data.sh
#
# Downloads publicly available data for benchmarking:
#   1. 1000 Genomes Omni 2.5M array genotypes (input for imputation)
#   2. 1000 Genomes 30x WGS data (truth for validation)
#   3. Sample metadata (population/ancestry info)
#
# This data allows benchmarking across multiple ancestry groups.
#
# Usage:
#   ./download_benchmark_data.sh -o /path/to/output
#
# Requirements:
#   - wget, curl
#   - bcftools, plink2
#   - ~100GB disk space (full dataset)
#   - ~20GB disk space (subset for quick testing)
#
################################################################################

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

OUTPUT_DIR="${OUTPUT_DIR:-./benchmarking/test_data}"
SUBSET_ONLY="${SUBSET_ONLY:-false}"  # Download only chr22 for quick testing
THREADS="${THREADS:-4}"

# 1000 Genomes FTP
IGSR_FTP="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp"

# Specific data URLs
OMNI_GENOTYPES="${IGSR_FTP}/release/20130502/supporting/hd_genotype_chip"
WGS_30X_BASE="ftp://ftp.sra.ebi.ac.uk/vol1/ERZ134/ERZ1346188"  # 1KG 30x
SAMPLE_INFO="${IGSR_FTP}/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Downloads benchmark data (1000 Genomes genotypes + WGS truth).

Options:
  -o, --output-dir DIR   Output directory (default: ./benchmarking/test_data)
  -s, --subset           Download only chr22 for quick testing (~5GB)
  -t, --threads N        Threads for processing (default: 4)
  -h, --help             Show this help

Data Downloaded:
----------------

1. Genotype Array Data (Omni 2.5M):
   - 2504 samples from 1000 Genomes Phase 3
   - 5 superpopulations: AFR, EUR, EAS, SAS, AMR
   - ~2.5M variants per sample

2. WGS Truth Data (30x):
   - Matched samples from 1000 Genomes
   - High-coverage for accurate truth genotypes
   - Used for concordance calculation

3. Sample Metadata:
   - Population assignments
   - Superpopulation (for GRAF-anc comparison)
   - Sex, family info

Disk Space Required:
  --subset:     ~5GB  (chr22 only, quick testing)
  Full dataset: ~100GB (all autosomes)

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -s|--subset)
            SUBSET_ONLY=true
            shift
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
            print_usage
            exit 1
            ;;
    esac
done

echo "=============================================="
echo "Benchmark Data Download"
echo "=============================================="
echo "Output directory: ${OUTPUT_DIR}"
echo "Subset only (chr22): ${SUBSET_ONLY}"
echo "Threads: ${THREADS}"
echo ""

mkdir -p "${OUTPUT_DIR}"/{genotypes,wgs_truth,sample_info}
cd "${OUTPUT_DIR}"

# =============================================================================
# 1. SAMPLE METADATA
# =============================================================================

echo ""
echo "=== 1/5: Downloading Sample Metadata ==="
echo ""

cd sample_info

if [ ! -f "integrated_call_samples.panel" ]; then
    echo "  Downloading 1000 Genomes sample panel..."
    wget -q -c "${SAMPLE_INFO}" -O integrated_call_samples.panel && \
        echo "  Downloaded sample panel" || \
        echo "  WARNING: Could not download sample panel"
fi

# Download superpopulation definitions
if [ ! -f "superpopulations.txt" ]; then
    cat > superpopulations.txt << 'SUPERPOPS'
pop	superpop	description
CHB	EAS	Han Chinese in Beijing, China
JPT	EAS	Japanese in Tokyo, Japan
CHS	EAS	Southern Han Chinese
CDX	EAS	Chinese Dai in Xishuangbanna, China
KHV	EAS	Kinh in Ho Chi Minh City, Vietnam
CEU	EUR	Utah Residents with Northern and Western European Ancestry
TSI	EUR	Toscani in Italia
FIN	EUR	Finnish in Finland
GBR	EUR	British in England and Scotland
IBS	EUR	Iberian Population in Spain
YRI	AFR	Yoruba in Ibadan, Nigeria
LWK	AFR	Luhya in Webuye, Kenya
GWD	AFR	Gambian in Western Divisions in the Gambia
MSL	AFR	Mende in Sierra Leone
ESN	AFR	Esan in Nigeria
ASW	AFR	Americans of African Ancestry in SW USA
ACB	AFR	African Caribbeans in Barbados
MXL	AMR	Mexican Ancestry from Los Angeles USA
PUR	AMR	Puerto Ricans from Puerto Rico
CLM	AMR	Colombians from Medellin, Colombia
PEL	AMR	Peruvians from Lima, Peru
GIH	SAS	Gujarati Indian from Houston, Texas
PJL	SAS	Punjabi from Lahore, Pakistan
BEB	SAS	Bengali from Bangladesh
STU	SAS	Sri Lankan Tamil from the UK
ITU	SAS	Indian Telugu from the UK
SUPERPOPS
    echo "  Created superpopulations.txt"
fi

# Create sample list per superpopulation (for subsetting)
echo "  Creating per-superpopulation sample lists..."
for pop in AFR EUR EAS SAS AMR; do
    awk -v pop="$pop" '$3 == pop {print $1}' integrated_call_samples.panel > "samples_${pop}.txt"
    n=$(wc -l < "samples_${pop}.txt")
    echo "    ${pop}: ${n} samples"
done

cd "${OUTPUT_DIR}"

# =============================================================================
# 2. GENOTYPE ARRAY DATA (OMNI 2.5M)
# =============================================================================

echo ""
echo "=== 2/5: Downloading Genotype Array Data ==="
echo ""

cd genotypes

# 1000 Genomes Omni data is available as PLINK files
OMNI_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes"

if [ "$SUBSET_ONLY" = "true" ]; then
    # Download only chr22 for quick testing
    CHROMOSOMES="22"
    echo "  Downloading chr22 only (subset mode)..."
else
    CHROMOSOMES=$(seq 1 22)
    echo "  Downloading all autosomes..."
fi

# Check if we need VCF or can get PLINK directly
# The 1KG Omni data is typically available as VCF
for chr in $CHROMOSOMES; do
    vcf_file="ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.chr${chr}.vcf.gz"

    if [ ! -f "${vcf_file}" ]; then
        echo "  Downloading chr${chr}..."
        wget -q -c "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/${vcf_file}" || \
            echo "  WARNING: Could not download chr${chr}"

        if [ -f "${vcf_file}" ]; then
            # Index the VCF
            bcftools index "${vcf_file}" 2>/dev/null || true
        fi
    else
        echo "  chr${chr} exists, skipping..."
    fi
done

# Merge chromosomes into single PLINK file
if [ ! -f "1kg_omni.bed" ]; then
    echo "  Converting to PLINK format..."

    # Create list of VCFs
    ls ALL.chip.omni*.vcf.gz 2>/dev/null | sort -V > vcf_list.txt

    if [ -s vcf_list.txt ]; then
        # Concatenate and convert
        bcftools concat -f vcf_list.txt -Oz -o 1kg_omni_merged.vcf.gz --threads "${THREADS}"

        plink2 --vcf 1kg_omni_merged.vcf.gz \
            --make-bed \
            --out 1kg_omni \
            --threads "${THREADS}"

        rm -f 1kg_omni_merged.vcf.gz vcf_list.txt
        echo "  Created 1kg_omni PLINK files"
    else
        echo "  WARNING: No VCF files found to merge"
    fi
fi

cd "${OUTPUT_DIR}"

# =============================================================================
# 3. WGS TRUTH DATA (30x)
# =============================================================================

echo ""
echo "=== 3/5: Downloading WGS Truth Data ==="
echo ""

cd wgs_truth

# 1000 Genomes 30x data is available from multiple sources
# Using the IGSR hosted version

# Base URL for 30x data
WGS_BASE="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot"

if [ "$SUBSET_ONLY" = "true" ]; then
    CHROMOSOMES="22"
    echo "  Downloading chr22 WGS only (subset mode)..."
else
    # For full benchmarking, still use chr22 initially (it's enough for validation)
    # Full WGS is very large (~1TB)
    CHROMOSOMES="22 21 20"  # Representative chromosomes
    echo "  Downloading representative chromosomes (20, 21, 22)..."
    echo "  (Full WGS is ~1TB - using subset for validation)"
fi

for chr in $CHROMOSOMES; do
    vcf_file="CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz"

    if [ ! -f "${vcf_file}" ]; then
        echo "  Downloading chr${chr} WGS..."
        wget -q -c "${WGS_BASE}/${vcf_file}" 2>/dev/null || \
            # Try alternative URL
            wget -q -c "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/${vcf_file}" || \
            echo "  WARNING: Could not download chr${chr} WGS"

        # Download index if available
        if [ -f "${vcf_file}" ]; then
            wget -q -c "${WGS_BASE}/${vcf_file}.tbi" 2>/dev/null || \
                bcftools index -t "${vcf_file}" 2>/dev/null || true
        fi
    else
        echo "  chr${chr} WGS exists, skipping..."
    fi
done

cd "${OUTPUT_DIR}"

# =============================================================================
# 4. CREATE BENCHMARK SAMPLE SHEETS
# =============================================================================

echo ""
echo "=== 4/5: Creating Benchmark Sample Sheets ==="
echo ""

# Create sample sheet for full benchmark
cat > benchmark_sample_sheet.csv << EOF
platform_id,batch_id,input_path,file_type,build,file_structure
1KG_Omni,AFR_samples,${OUTPUT_DIR}/genotypes/1kg_omni,plink,hg19,merged_batch
1KG_Omni,EUR_samples,${OUTPUT_DIR}/genotypes/1kg_omni,plink,hg19,merged_batch
1KG_Omni,EAS_samples,${OUTPUT_DIR}/genotypes/1kg_omni,plink,hg19,merged_batch
1KG_Omni,SAS_samples,${OUTPUT_DIR}/genotypes/1kg_omni,plink,hg19,merged_batch
1KG_Omni,AMR_samples,${OUTPUT_DIR}/genotypes/1kg_omni,plink,hg19,merged_batch
EOF

echo "  Created benchmark_sample_sheet.csv"

# Create sample sheet for quick test (subset)
cat > benchmark_sample_sheet_quick.csv << EOF
platform_id,batch_id,input_path,file_type,build,file_structure
1KG_Omni,ALL_samples,${OUTPUT_DIR}/genotypes/1kg_omni,plink,hg19,merged_batch
EOF

echo "  Created benchmark_sample_sheet_quick.csv"

cd "${OUTPUT_DIR}"

# =============================================================================
# 5. SIMULATE PHENOTYPES FOR GWAS/PRS TESTING
# =============================================================================

echo ""
echo "=== 5/5: Simulating Phenotypes for GWAS/PRS Testing ==="
echo ""

mkdir -p phenotypes

# Check if WGS truth VCF exists for phenotype simulation
WGS_VCF=$(ls wgs_truth/CCDG_*_chr22*.vcf.gz 2>/dev/null | head -1 || true)

if [ -n "$WGS_VCF" ] && [ -f "$WGS_VCF" ]; then
    echo "  Using WGS truth for phenotype simulation: $(basename $WGS_VCF)"

    # Check if R and required packages are available
    if command -v Rscript &> /dev/null; then
        SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

        # Simulate height (continuous trait, ~50% heritable)
        echo "  Simulating HEIGHT phenotype..."
        Rscript "${SCRIPT_DIR}/bench-helper-scripts/simulate_phenotypes.R" \
            --vcf "$WGS_VCF" \
            --trait height \
            --h2 0.5 \
            --output phenotypes/height_simulated \
            2>/dev/null || echo "    (Skipped - missing R packages)"

        # Simulate T2D (binary trait, tests ancestry-specific effects)
        echo "  Simulating T2D phenotype..."
        Rscript "${SCRIPT_DIR}/bench-helper-scripts/simulate_phenotypes.R" \
            --vcf "$WGS_VCF" \
            --trait t2d \
            --h2 0.4 \
            --prevalence 0.1 \
            --output phenotypes/t2d_simulated \
            2>/dev/null || echo "    (Skipped - missing R packages)"

        # Simulate LDL (tests rare variant detection)
        echo "  Simulating LDL phenotype..."
        Rscript "${SCRIPT_DIR}/bench-helper-scripts/simulate_phenotypes.R" \
            --vcf "$WGS_VCF" \
            --trait ldl \
            --h2 0.3 \
            --output phenotypes/ldl_simulated \
            2>/dev/null || echo "    (Skipped - missing R packages)"

        if [ -f "phenotypes/height_simulated.pheno" ]; then
            echo "  Created simulated phenotypes:"
            echo "    - phenotypes/height_simulated.pheno (continuous)"
            echo "    - phenotypes/t2d_simulated.pheno (binary)"
            echo "    - phenotypes/ldl_simulated.pheno (continuous)"
            echo "    - *_causal_variants.txt (for hit recovery testing)"
        fi
    else
        echo "  Rscript not found - skipping phenotype simulation"
        echo "  Run manually: Rscript bench-helper-scripts/simulate_phenotypes.R --help"
    fi
else
    echo "  WGS truth not found - skipping phenotype simulation"
    echo "  Download WGS data first, then run:"
    echo "    Rscript bench-helper-scripts/simulate_phenotypes.R --vcf <wgs.vcf.gz> --trait height"
fi

cd "${OUTPUT_DIR}"

# =============================================================================
# CREATE SUMMARY
# =============================================================================

cat > BENCHMARK_DATA_README.txt << EOF
================================================================================
Benchmark Test Data
================================================================================

Downloaded: $(date)

Directory Structure:
--------------------

genotypes/
├── 1kg_omni.bed/bim/fam          # 1000 Genomes Omni 2.5M array data
└── ALL.chip.omni*.vcf.gz         # Per-chromosome VCFs (if downloaded)

wgs_truth/
└── CCDG_*_chr*.vcf.gz            # 30x WGS truth data

sample_info/
├── integrated_call_samples.panel  # Sample metadata
├── superpopulations.txt           # Population definitions
└── samples_*.txt                  # Per-superpopulation sample lists

phenotypes/                        # Simulated phenotypes for GWAS/PRS testing
├── height_simulated.pheno         # Height (continuous, h2=0.5)
├── t2d_simulated.pheno            # Type 2 Diabetes (binary, h2=0.4)
├── ldl_simulated.pheno            # LDL cholesterol (continuous, h2=0.3)
└── *_causal_variants.txt          # Known causal variants for hit recovery


Sample Counts by Superpopulation:
---------------------------------

$(for pop in AFR EUR EAS SAS AMR; do
    if [ -f "sample_info/samples_${pop}.txt" ]; then
        n=$(wc -l < "sample_info/samples_${pop}.txt")
        echo "  ${pop}: ${n} samples"
    fi
done)


Usage:
------

1. Run GRAF-anc to verify ancestry classifications:

   nextflow run main.nf \\
       --sample_sheet benchmark_sample_sheet.csv \\
       --skip_modules "1,2,3,4,5,6" \\
       --outdir benchmarking/results/ancestry_check

2. Run benchmark approaches (see benchmarking/alternative_approaches/)

3. Calculate concordance with WGS truth:

   Rscript benchmarking/bench-helper-scripts/calculate_concordance.R \\
       --imputed results/approach_a/final.vcf.gz \\
       --truth wgs_truth/CCDG_*_chr22.vcf.gz \\
       --output results/concordance_approach_a.txt

4. Run GWAS with simulated phenotypes:

   plink2 --bfile imputed_data \\
       --pheno phenotypes/height_simulated.pheno \\
       --glm \\
       --out gwas_height

5. Check causal variant recovery:

   Rscript benchmarking/bench-helper-scripts/check_hit_recovery.R \\
       --gwas gwas_height.glm.linear \\
       --causal phenotypes/height_simulated_causal_variants.txt \\
       --output hit_recovery_height


Notes:
------

- WGS truth data is limited to chr20-22 by default (full WGS is ~1TB)
- For publication-quality benchmarks, consider downloading more chromosomes
- Sample overlap between array and WGS should be verified before concordance

================================================================================
EOF

echo ""
echo "=============================================="
echo "Benchmark Data Download Complete!"
echo "=============================================="
echo ""
echo "Downloaded to: ${OUTPUT_DIR}"
echo ""
echo "Next steps:"
echo "  1. Run GRAF-anc to classify samples by ancestry"
echo "  2. Run alternative pipeline approaches"
echo "  3. Calculate concordance metrics"
echo ""
echo "See BENCHMARK_DATA_README.txt for details."
echo ""
