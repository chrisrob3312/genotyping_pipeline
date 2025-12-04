#!/bin/bash
################################################################################
# setup_reference_files.sh
#
# Downloads and organizes reference files needed for the genotyping pipeline:
#   - Chain files (liftover between hg19/hg38)
#   - BRAVO Freeze 10 VCFs → Rayner format for imputation QC
#   - Genetic maps from TractorWorkflow (SHAPEIT5/RFMix2 formats)
#   - Reference genome FASTA files
#   - HRC reference for hg19 fallback
#
# NOTE: LAI reference panel (HGDP-1KG + MX Biobank) is user-provided
#       Use format_reference_panel.sh to convert your custom reference
#
# Usage:
#   ./setup_reference_files.sh --output-dir /path/to/references
#
# Requirements:
#   - wget, curl
#   - bcftools (for BRAVO VCF processing)
#   - ~100GB disk space for all references
#
################################################################################

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

OUTPUT_DIR="${OUTPUT_DIR:-./pipeline_references}"
THREADS="${THREADS:-4}"
SKIP_BRAVO="${SKIP_BRAVO:-false}"

# BRAVO Freeze 10 base URL (requires visiting https://bravo.sph.umich.edu/vcfs.html for signed URLs)
BRAVO_BASE="https://bravo.sph.umich.edu"

# TractorWorkflow reference data
TRACTOR_REPO="https://github.com/Atkinson-Lab/TractorWorkflow"
TRACTOR_TESTDATA="https://github.com/Atkinson-Lab/Tractor-tutorial/raw/refs/heads/main/test_data.zip"

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Downloads reference files for the genotyping pipeline.

Options:
  -o, --output-dir DIR   Output directory (default: ./pipeline_references)
  -t, --threads N        Threads for processing (default: 4)
  --skip-bravo           Skip BRAVO download (use if already downloaded)
  -h, --help             Show this help

Reference Sources:
  - Chain files:    UCSC Genome Browser
  - TOPMed (hg38):  BRAVO Freeze 10 (https://bravo.sph.umich.edu/vcfs.html)
  - HRC (hg19):     Will Rayner's site
  - Genetic maps:   TractorWorkflow / SHAPEIT5
  - FASTA:          1000 Genomes / NCBI

NOTE: LAI reference (HGDP-1KG + custom samples) must be formatted separately
      using format_reference_panel.sh with your WGS reference data.

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --skip-bravo)
            SKIP_BRAVO=true
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
echo "Genotyping Pipeline Reference File Setup"
echo "=============================================="
echo "Output directory: ${OUTPUT_DIR}"
echo "Threads: ${THREADS}"
echo ""

mkdir -p "${OUTPUT_DIR}"/{chain_files,bravo_vcfs,rayner_refs,genetic_maps,fasta,tractor_refs}

cd "${OUTPUT_DIR}"

# =============================================================================
# 1. CHAIN FILES (Liftover)
# =============================================================================

echo ""
echo "=== 1. Downloading Chain Files ==="
echo ""

cd chain_files

for chain in "hg19ToHg38" "hg38ToHg19"; do
    if [ ! -f "${chain}.over.chain.gz" ]; then
        build=$(echo $chain | cut -d'T' -f1)
        echo "Downloading ${chain}.over.chain.gz..."
        wget -q -c "https://hgdownload.soe.ucsc.edu/goldenPath/${build}/liftOver/${chain}.over.chain.gz" || \
            echo "Warning: Could not download ${chain}"
    else
        echo "${chain}.over.chain.gz exists, skipping..."
    fi
done

cd "${OUTPUT_DIR}"

# =============================================================================
# 2. BRAVO FREEZE 10 (TOPMed hg38)
# =============================================================================

echo ""
echo "=== 2. BRAVO Freeze 10 (TOPMed hg38 Reference) ==="
echo ""

if [ "$SKIP_BRAVO" = "true" ]; then
    echo "Skipping BRAVO download (--skip-bravo flag set)"
else
    cat << 'BRAVO_INSTRUCTIONS'
--------------------------------------------------------------------------------
BRAVO Freeze 10 VCF Download Instructions
--------------------------------------------------------------------------------

The BRAVO VCF files require signed URLs that expire after 15 minutes.
You must download them manually from:

    https://bravo.sph.umich.edu/vcfs.html

Steps:
1. Visit the URL above in your browser
2. Right-click each chromosome link and copy the URL
3. Download using wget/curl (URLs are time-limited)

Example for chromosome 1:
    wget -O bravo_vcfs/chr1.bravo.pub.vcf.gz "PASTE_SIGNED_URL_HERE"

After downloading all chromosomes (1-22, X), run:
    ./helper_scripts/convert_bravo_to_rayner.sh \
        --input-dir bravo_vcfs/ \
        --output rayner_refs/PASS.Variants.TOPMed_freeze10_hg38.tab.gz

Alternative: Use the existing Freeze 5 reference from Will Rayner's site:
    wget -P rayner_refs/ https://www.chg.ox.ac.uk/~wrayner/tools/PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz

--------------------------------------------------------------------------------
BRAVO_INSTRUCTIONS

    # Create download helper script
    cat > bravo_vcfs/download_bravo.sh << 'DOWNLOAD_SCRIPT'
#!/bin/bash
# Download BRAVO Freeze 10 VCFs
# Run this script after getting signed URLs from https://bravo.sph.umich.edu/vcfs.html

echo "Paste the signed URLs for each chromosome when prompted"
echo "URLs expire after 15 minutes, so download quickly"
echo ""

for chr in {1..22} X; do
    if [ ! -f "chr${chr}.bravo.pub.vcf.gz" ]; then
        echo "Enter URL for chromosome ${chr}:"
        read -r url
        if [ -n "$url" ]; then
            wget -O "chr${chr}.bravo.pub.vcf.gz" "$url"
        fi
    else
        echo "chr${chr}.bravo.pub.vcf.gz exists, skipping..."
    fi
done

echo "Download complete. Now run convert_bravo_to_rayner.sh"
DOWNLOAD_SCRIPT
    chmod +x bravo_vcfs/download_bravo.sh
fi

# Download Freeze 5 as fallback/alternative
cd rayner_refs
if [ ! -f "PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz" ]; then
    echo "Downloading TOPMed Freeze 5 reference (fallback)..."
    wget -q -c "https://www.chg.ox.ac.uk/~wrayner/tools/PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz" || \
        echo "Warning: Could not download Freeze 5 reference"
fi

cd "${OUTPUT_DIR}"

# =============================================================================
# 3. HRC REFERENCE (hg19/GRCh37)
# =============================================================================

echo ""
echo "=== 3. HRC Reference (hg19/GRCh37) ==="
echo ""

cd rayner_refs

if [ ! -f "HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz" ]; then
    echo "Downloading HRC r1.1 reference..."
    wget -q -c "ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz" 2>/dev/null || \
    wget -q -c "https://www.chg.ox.ac.uk/~wrayner/tools/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz" || \
        echo "Warning: Could not download HRC reference"
else
    echo "HRC reference exists, skipping..."
fi

cd "${OUTPUT_DIR}"

# =============================================================================
# 4. GENETIC MAPS (from TractorWorkflow)
# =============================================================================

echo ""
echo "=== 4. Genetic Maps (TractorWorkflow / SHAPEIT5) ==="
echo ""

cd genetic_maps

# Download TractorWorkflow test data which includes genetic maps
if [ ! -d "tractor_maps" ]; then
    echo "Downloading TractorWorkflow reference files..."

    wget -q -c "$TRACTOR_TESTDATA" -O tractor_test_data.zip && \
        unzip -q -o tractor_test_data.zip && \
        mv test_data/references tractor_maps && \
        rm -rf test_data tractor_test_data.zip || \
        echo "Warning: Could not download TractorWorkflow data"
fi

# Also download SHAPEIT5 genetic maps for all chromosomes
mkdir -p shapeit5_b37 shapeit5_b38

echo "Downloading SHAPEIT5 genetic maps (b37)..."
wget -q -c "https://github.com/odelaneau/shapeit5/raw/main/maps/genetic_maps.b37.tar.gz" -O shapeit5_b37.tar.gz && \
    tar -xzf shapeit5_b37.tar.gz -C shapeit5_b37 --strip-components=1 && \
    rm shapeit5_b37.tar.gz || \
    echo "Warning: Could not download b37 maps"

echo "Downloading SHAPEIT5 genetic maps (b38)..."
wget -q -c "https://github.com/odelaneau/shapeit5/raw/main/maps/genetic_maps.b38.tar.gz" -O shapeit5_b38.tar.gz && \
    tar -xzf shapeit5_b38.tar.gz -C shapeit5_b38 --strip-components=1 && \
    rm shapeit5_b38.tar.gz || \
    echo "Warning: Could not download b38 maps"

cd "${OUTPUT_DIR}"

# =============================================================================
# 5. REFERENCE FASTA FILES
# =============================================================================

echo ""
echo "=== 5. Reference FASTA Files ==="
echo ""

cd fasta

# GRCh37/hg19
mkdir -p hg19
cd hg19
if [ ! -f "human_g1k_v37.fasta" ]; then
    echo "Downloading GRCh37 reference..."
    wget -q -c "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz" && \
        gunzip human_g1k_v37.fasta.gz || \
        echo "Warning: Could not download GRCh37 reference"

    if command -v samtools &> /dev/null && [ -f human_g1k_v37.fasta ]; then
        echo "Creating FASTA index..."
        samtools faidx human_g1k_v37.fasta
    fi
fi
cd ..

# GRCh38/hg38
mkdir -p hg38
cd hg38
if [ ! -f "GRCh38_full_analysis_set_plus_decoy_hla.fa" ]; then
    echo "Downloading GRCh38 reference..."
    wget -q -c "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa" || \
        echo "Warning: Could not download GRCh38 reference"

    if command -v samtools &> /dev/null && [ -f GRCh38_full_analysis_set_plus_decoy_hla.fa ]; then
        echo "Creating FASTA index..."
        samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa
    fi
fi
cd "${OUTPUT_DIR}"

# =============================================================================
# 6. CREATE SUMMARY
# =============================================================================

echo ""
echo "=== Creating Summary ==="

cat > README.txt << 'EOF'
================================================================================
Pipeline Reference Files
================================================================================

This directory contains reference files for the genotyping pipeline.

Directory Structure:
--------------------

chain_files/
├── hg19ToHg38.over.chain.gz      # Liftover hg19 → hg38
└── hg38ToHg19.over.chain.gz      # Liftover hg38 → hg19

bravo_vcfs/
├── download_bravo.sh             # Helper script for BRAVO download
└── chr*.bravo.pub.vcf.gz         # BRAVO Freeze 10 VCFs (manual download)

rayner_refs/
├── PASS.Variants.TOPMed_freeze10_hg38.tab.gz  # TOPMed Fr10 (from BRAVO)
├── PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz  # TOPMed Fr5 (fallback)
└── HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz      # HRC for hg19

genetic_maps/
├── shapeit5_b37/                 # SHAPEIT5 genetic maps (hg19)
├── shapeit5_b38/                 # SHAPEIT5 genetic maps (hg38)
└── tractor_maps/                 # TractorWorkflow format maps
    ├── chr*.b37.gmap.gz          # SHAPEIT5 format
    └── chr*.genetic_map.modified.txt  # RFMix2/GNomix format

fasta/
├── hg19/human_g1k_v37.fasta      # GRCh37 reference
└── hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa  # GRCh38 reference


BRAVO Freeze 10 Setup:
----------------------

1. Visit https://bravo.sph.umich.edu/vcfs.html
2. Download all chromosome VCFs to bravo_vcfs/
3. Run: ./helper_scripts/convert_bravo_to_rayner.sh \
        --input-dir bravo_vcfs/ \
        --output rayner_refs/PASS.Variants.TOPMed_freeze10_hg38.tab.gz


LAI Reference Panel (User-Provided):
------------------------------------

Your HGDP-1KG + MX Biobank reference is NOT included here.
Format it using:

    ./helper_scripts/format_reference_panel.sh \
        --input /path/to/your_reference.vcf.gz \
        --sample-map /path/to/sample_populations.txt \
        --genetic-map-dir genetic_maps/shapeit5_b38/ \
        --output-dir ancestry_references/


Copy to Pipeline:
-----------------

cp -r chain_files/* /path/to/pipeline/resources/chain_files/
cp -r rayner_refs/* /path/to/pipeline/resources/rayner/
cp -r genetic_maps/* /path/to/pipeline/resources/genetic_maps/
cp -r fasta/* /path/to/pipeline/resources/references/

================================================================================
EOF

echo ""
echo "=============================================="
echo "Reference File Setup Complete!"
echo "=============================================="
echo ""
echo "Downloaded files:"
find "${OUTPUT_DIR}" -type f \( -name "*.gz" -o -name "*.fa" -o -name "*.fasta" \) 2>/dev/null | head -20
echo ""
echo "Total size:"
du -sh "${OUTPUT_DIR}"
echo ""
echo "IMPORTANT: For TOPMed Freeze 10 (BRAVO), follow instructions in:"
echo "  ${OUTPUT_DIR}/README.txt"
echo ""
echo "Next steps:"
echo "  1. Download BRAVO Freeze 10 VCFs (see bravo_vcfs/download_bravo.sh)"
echo "  2. Convert BRAVO to Rayner format (see convert_bravo_to_rayner.sh)"
echo "  3. Format your LAI reference (HGDP-1KG + MX Biobank) with format_reference_panel.sh"
echo ""
