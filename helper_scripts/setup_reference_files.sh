#!/bin/bash
################################################################################
# setup_reference_files.sh
#
# Downloads and organizes reference files needed for the genotyping pipeline:
#   - Chain files (liftover between hg19/hg38)
#   - TOPMed/HRC reference files for imputation QC (Rayner/McCarthy format)
#   - Genetic maps for phasing and LAI
#   - Reference genome FASTA files
#
# Run on HPC, then copy formatted data to pipeline resources/ directory
#
# Usage:
#   ./setup_reference_files.sh --output-dir /path/to/references
#
# Requirements:
#   - wget or curl
#   - gunzip, unzip
#   - ~50GB disk space for all references
#
################################################################################

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

# Default output directory
OUTPUT_DIR="${OUTPUT_DIR:-./pipeline_references}"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [--output-dir DIR]"
            echo ""
            echo "Downloads reference files for the genotyping pipeline:"
            echo "  - Chain files (hg19 <-> hg38)"
            echo "  - TOPMed/HRC strand files (Rayner format)"
            echo "  - Genetic maps"
            echo "  - Reference FASTA files"
            echo ""
            echo "Options:"
            echo "  -o, --output-dir DIR   Output directory (default: ./pipeline_references)"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

echo "=============================================="
echo "Genotyping Pipeline Reference File Setup"
echo "=============================================="
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Create directory structure
mkdir -p "${OUTPUT_DIR}"/{chain_files,strand_files,genetic_maps,fasta,checksums}

cd "${OUTPUT_DIR}"

# =============================================================================
# 1. CHAIN FILES (Liftover)
# =============================================================================

echo ""
echo "=== Downloading Chain Files ==="
echo ""

cd chain_files

# hg19 to hg38
if [ ! -f hg19ToHg38.over.chain.gz ]; then
    echo "Downloading hg19 to hg38 chain file..."
    wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
else
    echo "hg19ToHg38.over.chain.gz already exists, skipping..."
fi

# hg38 to hg19
if [ ! -f hg38ToHg19.over.chain.gz ]; then
    echo "Downloading hg38 to hg19 chain file..."
    wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
else
    echo "hg38ToHg19.over.chain.gz already exists, skipping..."
fi

# GRCh37 to GRCh38 (Ensembl - sometimes needed for VCF liftover)
if [ ! -f GRCh37_to_GRCh38.chain.gz ]; then
    echo "Downloading GRCh37 to GRCh38 chain file (Ensembl)..."
    wget -c https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz || \
        echo "Warning: Ensembl chain file not available, using UCSC version"
fi

cd "${OUTPUT_DIR}"

# =============================================================================
# 2. STRAND/REFERENCE FILES (Rayner/McCarthy Format for Imputation QC)
# =============================================================================

echo ""
echo "=== Downloading Strand Files (Rayner/McCarthy Format) ==="
echo ""

cd strand_files

# -----------------------------------------------------------------------------
# TOPMed Reference (hg38) - Primary for TOPMed imputation
# -----------------------------------------------------------------------------
echo ""
echo "--- TOPMed Reference Files (hg38) ---"

mkdir -p topmed_hg38
cd topmed_hg38

# TOPMed reference allele frequencies and strand info
# From: https://www.well.ox.ac.uk/~wrayner/tools/
if [ ! -f PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz ]; then
    echo "Downloading TOPMed freeze5 hg38 reference..."
    wget -c https://www.well.ox.ac.uk/~wrayner/tools/PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz || {
        echo "Direct download failed, trying alternative..."
        # Alternative: Create from TOPMed VCF if direct download fails
        echo "NOTE: You may need to generate this from TOPMed reference VCF"
        echo "See: https://www.well.ox.ac.uk/~wrayner/tools/"
    }
fi

cd "${OUTPUT_DIR}/strand_files"

# -----------------------------------------------------------------------------
# HRC Reference (hg19/GRCh37) - For HRC imputation
# -----------------------------------------------------------------------------
echo ""
echo "--- HRC Reference Files (hg19/GRCh37) ---"

mkdir -p hrc_hg19
cd hrc_hg19

if [ ! -f HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz ]; then
    echo "Downloading HRC r1.1 reference..."
    wget -c ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz || {
        echo "Primary source failed, trying alternative..."
        wget -c https://www.well.ox.ac.uk/~wrayner/tools/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz || \
            echo "Warning: Could not download HRC reference"
    }
fi

cd "${OUTPUT_DIR}/strand_files"

# -----------------------------------------------------------------------------
# 1000 Genomes Reference (hg19) - Fallback reference
# -----------------------------------------------------------------------------
echo ""
echo "--- 1000 Genomes Reference Files (hg19) ---"

mkdir -p 1kg_hg19
cd 1kg_hg19

if [ ! -f 1000GP_Phase3_combined.legend.gz ]; then
    echo "Downloading 1000G Phase 3 legend files..."
    # These are used by some tools for strand checking
    for chr in {1..22}; do
        wget -c https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/1000GP_Phase3_chr${chr}.legend.gz || true
    done
fi

cd "${OUTPUT_DIR}/strand_files"

# -----------------------------------------------------------------------------
# Illumina/Affymetrix Strand Files
# -----------------------------------------------------------------------------
echo ""
echo "--- Array-Specific Strand Files ---"

mkdir -p array_strand_files
cd array_strand_files

# Download Will Rayner's strand files for common arrays
# https://www.well.ox.ac.uk/~wrayner/strand/

STRAND_BASE="https://www.well.ox.ac.uk/~wrayner/strand"

# Common Illumina arrays
declare -a ILLUMINA_ARRAYS=(
    "HumanOmni2.5-8v1-Multi"
    "HumanOmni2.5-8v1_A"
    "HumanOmniExpress-24v1-0_A"
    "HumanOmniExpress-24v1-1_A"
    "Multi-EthnicGlobal_A1"
    "InfiniumOmni5-4v1-2_A1"
    "HumanCoreExome-24v1-0_A"
    "HumanCore-24v1-0_A"
    "GSA-24v3-0_A1"
    "GSA-24v2-0_A1"
)

# Common Affymetrix arrays
declare -a AFFY_ARRAYS=(
    "Axiom_PMRA.na35"
    "Axiom_UKB_WCSG"
    "Axiom_GW_ASI_SNP"
    "Axiom_GW_AFR_SNP"
    "Axiom_GW_LAT_SNP"
    "Axiom_GW_EU_SNP"
    "GenomeWideSNP_6"
)

echo "Downloading Illumina strand files..."
for array in "${ILLUMINA_ARRAYS[@]}"; do
    if [ ! -f "${array}-b37-strand.zip" ] && [ ! -f "${array}-b37.strand" ]; then
        wget -c "${STRAND_BASE}/${array}-b37-strand.zip" 2>/dev/null || \
            echo "  Note: ${array} strand file not available"
    fi
done

echo "Downloading Affymetrix strand files..."
for array in "${AFFY_ARRAYS[@]}"; do
    if [ ! -f "${array}-b37-strand.zip" ] && [ ! -f "${array}.strand" ]; then
        wget -c "${STRAND_BASE}/${array}-b37-strand.zip" 2>/dev/null || \
            echo "  Note: ${array} strand file not available"
    fi
done

# Unzip strand files
for zip in *.zip; do
    [ -f "$zip" ] && unzip -o "$zip" && rm "$zip"
done 2>/dev/null || true

cd "${OUTPUT_DIR}"

# =============================================================================
# 3. GENETIC MAPS (for Phasing and LAI)
# =============================================================================

echo ""
echo "=== Downloading Genetic Maps ==="
echo ""

cd genetic_maps

# -----------------------------------------------------------------------------
# SHAPEIT4/Eagle genetic maps (GRCh37/hg19)
# -----------------------------------------------------------------------------
echo "--- GRCh37/hg19 Genetic Maps ---"

mkdir -p hg19
cd hg19

if [ ! -f genetic_map_chr1_combined_b37.txt ]; then
    echo "Downloading SHAPEIT4 genetic maps (b37)..."

    # Try SHAPEIT4 GitHub first
    wget -c https://github.com/odelaneau/shapeit4/raw/master/maps/genetic_maps.b37.tar.gz 2>/dev/null && \
        tar -xzf genetic_maps.b37.tar.gz && rm genetic_maps.b37.tar.gz || {

        # Fallback to Eagle maps
        echo "Trying Eagle genetic maps..."
        wget -c https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg19_withX.txt.gz && \
            gunzip genetic_map_hg19_withX.txt.gz
    }
fi

# Also get HapMap genetic maps (used by some older tools)
if [ ! -f genetic_map_GRCh37_chr1.txt ]; then
    echo "Downloading HapMap genetic maps..."
    wget -c https://github.com/odelaneau/shapeit4/raw/master/maps/genetic_maps.b37.tar.gz 2>/dev/null || true
fi

cd "${OUTPUT_DIR}/genetic_maps"

# -----------------------------------------------------------------------------
# GRCh38/hg38 Genetic Maps
# -----------------------------------------------------------------------------
echo "--- GRCh38/hg38 Genetic Maps ---"

mkdir -p hg38
cd hg38

if [ ! -f genetic_map_chr1_combined_b38.txt ]; then
    echo "Downloading SHAPEIT4 genetic maps (b38)..."

    wget -c https://github.com/odelaneau/shapeit4/raw/master/maps/genetic_maps.b38.tar.gz 2>/dev/null && \
        tar -xzf genetic_maps.b38.tar.gz && rm genetic_maps.b38.tar.gz || {

        echo "Trying Eagle hg38 maps..."
        wget -c https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz && \
            gunzip genetic_map_hg38_withX.txt.gz
    }
fi

cd "${OUTPUT_DIR}"

# =============================================================================
# 4. REFERENCE FASTA FILES
# =============================================================================

echo ""
echo "=== Downloading Reference FASTA Files ==="
echo ""

cd fasta

# -----------------------------------------------------------------------------
# GRCh37/hg19 Reference
# -----------------------------------------------------------------------------
echo "--- GRCh37/hg19 Reference ---"

mkdir -p hg19
cd hg19

if [ ! -f human_g1k_v37.fasta ]; then
    echo "Downloading GRCh37 reference (1000 Genomes version)..."
    wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz && \
        gunzip human_g1k_v37.fasta.gz

    # Create index
    if command -v samtools &> /dev/null; then
        echo "Creating FASTA index..."
        samtools faidx human_g1k_v37.fasta
    fi
fi

cd "${OUTPUT_DIR}/fasta"

# -----------------------------------------------------------------------------
# GRCh38/hg38 Reference
# -----------------------------------------------------------------------------
echo "--- GRCh38/hg38 Reference ---"

mkdir -p hg38
cd hg38

if [ ! -f GRCh38_full_analysis_set_plus_decoy_hla.fa ]; then
    echo "Downloading GRCh38 reference..."

    # Try 1000 Genomes version first
    wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa 2>/dev/null || {

        # Alternative: NCBI
        echo "Trying NCBI reference..."
        wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz && \
            gunzip -c GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > GRCh38_full_analysis_set_plus_decoy_hla.fa
    }

    # Create index
    if command -v samtools &> /dev/null; then
        echo "Creating FASTA index..."
        samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa
    fi
fi

cd "${OUTPUT_DIR}"

# =============================================================================
# 5. DOWNLOAD CHECKING TOOLS (McCarthy/Rayner)
# =============================================================================

echo ""
echo "=== Downloading QC Tools ==="
echo ""

mkdir -p tools
cd tools

# HRC/1000G Imputation preparation checking tool
if [ ! -f HRC-1000G-check-bim-v4.3.0.zip ]; then
    echo "Downloading McCarthy/Rayner checking tool..."
    wget -c https://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip && \
        unzip -o HRC-1000G-check-bim-v4.3.0.zip
fi

cd "${OUTPUT_DIR}"

# =============================================================================
# 6. CREATE SUMMARY
# =============================================================================

echo ""
echo "=== Creating Summary ==="
echo ""

cat > REFERENCE_FILES_README.txt << 'EOF'
Pipeline Reference Files
========================

This directory contains reference files for the genotyping pipeline.

Directory Structure:
--------------------

chain_files/
├── hg19ToHg38.over.chain.gz      # Liftover hg19 → hg38
├── hg38ToHg19.over.chain.gz      # Liftover hg38 → hg19
└── GRCh37_to_GRCh38.chain.gz     # Ensembl chain file

strand_files/
├── topmed_hg38/                   # TOPMed reference (primary for imputation)
│   └── PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz
├── hrc_hg19/                      # HRC reference (GRCh37)
│   └── HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
├── 1kg_hg19/                      # 1000 Genomes legend files
│   └── 1000GP_Phase3_chr*.legend.gz
└── array_strand_files/            # Array-specific strand files
    └── *.strand

genetic_maps/
├── hg19/                          # GRCh37 genetic maps
│   └── genetic_map_*_combined_b37.txt
└── hg38/                          # GRCh38 genetic maps
    └── genetic_map_*_combined_b38.txt

fasta/
├── hg19/
│   ├── human_g1k_v37.fasta       # GRCh37 reference
│   └── human_g1k_v37.fasta.fai
└── hg38/
    ├── GRCh38_full_analysis_set_plus_decoy_hla.fa
    └── GRCh38_full_analysis_set_plus_decoy_hla.fa.fai

tools/
└── HRC-1000G-check-bim-v4.3.0/   # McCarthy/Rayner checking tool


Copy to Pipeline:
-----------------

# Copy to your pipeline resources directory:
cp -r chain_files/* /path/to/pipeline/resources/chain_files/
cp -r strand_files/* /path/to/pipeline/resources/strand_files/
cp -r genetic_maps/* /path/to/pipeline/resources/genetic_maps/
cp -r fasta/* /path/to/pipeline/resources/references/


References:
-----------
- TOPMed: https://imputation.biodatacatalyst.nhlbi.nih.gov/
- HRC: http://www.haplotype-reference-consortium.org/
- Rayner tools: https://www.well.ox.ac.uk/~wrayner/tools/
- Strand files: https://www.well.ox.ac.uk/~wrayner/strand/
EOF

echo ""
echo "=============================================="
echo "Reference File Setup Complete!"
echo "=============================================="
echo ""
echo "Files downloaded to: ${OUTPUT_DIR}"
echo ""
echo "Summary saved to: ${OUTPUT_DIR}/REFERENCE_FILES_README.txt"
echo ""
echo "Next steps:"
echo "  1. Copy files to pipeline resources/ directory"
echo "  2. Run format_reference_panel.sh on your WGS reference data"
echo "  3. Run create_test_data.sh to generate test datasets"
echo ""

# List what was downloaded
echo "Downloaded files:"
find "${OUTPUT_DIR}" -type f -name "*.gz" -o -name "*.fa" -o -name "*.fasta" -o -name "*.strand" -o -name "*.tab" 2>/dev/null | head -30
echo ""
echo "Total size:"
du -sh "${OUTPUT_DIR}"
