#!/bin/bash
################################################################################
# download_core_references.sh
#
# Downloads ALL publicly available reference files for Modules 1-7:
#   - Chain files (liftover between hg19/hg38)
#   - TOPMed Freeze 10 reference (BRAVO - required for imputation QC)
#   - HRC reference for hg19 fallback
#   - Reference genome FASTA files (hg19/hg38)
#   - Genetic maps (SHAPEIT5 b37/b38 + TractorWorkflow format)
#
# IMPORTANT: These are ALL publicly available - NO custom data required!
#
# The ONLY thing that requires your custom HGDP-1KG + MX Biobank data is
# the LAI reference panel, which you prepare separately using:
#   ./format_reference_panel.sh
#
# Usage:
#   ./download_core_references.sh --output-dir /path/to/references
#
# Requirements:
#   - wget, curl, unzip
#   - bcftools (for BRAVO Freeze 10 conversion)
#   - samtools (for FASTA indexing, optional)
#   - ~60GB disk space
#
################################################################################

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

OUTPUT_DIR="${OUTPUT_DIR:-./pipeline_references}"
SKIP_FASTA="${SKIP_FASTA:-false}"
SKIP_GENETIC_MAPS="${SKIP_GENETIC_MAPS:-false}"

# TractorWorkflow reference data (genetic maps only - no custom panel needed)
TRACTOR_TESTDATA="https://github.com/Atkinson-Lab/Tractor-tutorial/raw/refs/heads/main/test_data.zip"

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Downloads ALL publicly available reference files for Modules 1-7.

Options:
  -o, --output-dir DIR   Output directory (default: ./pipeline_references)
  --skip-fasta           Skip FASTA download (if you already have them)
  --skip-genetic-maps    Skip genetic maps download
  -h, --help             Show this help

What This Script Downloads (ALL PUBLIC - no custom data needed):
----------------------------------------------------------------

For Modules 1-6 (Imputation Pipeline):
  [x] Chain files (hg19 <-> hg38 liftover)
  [x] TOPMed Freeze 10 reference (BRAVO) - REQUIRED
  [x] HRC r1.1 reference (hg19)
  [x] Reference FASTA files (hg19, hg38)

For Module 7 (Ancestry - genetic maps only):
  [x] SHAPEIT5 genetic maps (b37 and b38)
  [x] TractorWorkflow genetic maps (RFMix2/GNomix format)

What You Must Prepare Separately (requires YOUR data):
------------------------------------------------------
  [ ] LAI reference panel (HGDP-1KG + MX Biobank)
      -> Use format_reference_panel.sh after obtaining your WGS reference

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --skip-fasta)
            SKIP_FASTA=true
            shift
            ;;
        --skip-genetic-maps)
            SKIP_GENETIC_MAPS=true
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
echo "Pipeline Reference Files Download"
echo "=============================================="
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "All files are PUBLICLY available."
echo "NO custom reference panel needed for this script."
echo ""

mkdir -p "${OUTPUT_DIR}"/{chain_files,bravo_vcfs,rayner_refs,fasta,genetic_maps}

cd "${OUTPUT_DIR}"

# =============================================================================
# 1. CHAIN FILES (Liftover)
# =============================================================================

echo ""
echo "=== 1/5: Downloading Chain Files ==="
echo ""

cd chain_files

for chain in "hg19ToHg38" "hg38ToHg19"; do
    if [ ! -f "${chain}.over.chain.gz" ]; then
        build=$(echo $chain | cut -d'T' -f1)
        echo "  Downloading ${chain}.over.chain.gz..."
        wget -q -c "https://hgdownload.soe.ucsc.edu/goldenPath/${build}/liftOver/${chain}.over.chain.gz" || \
            echo "  WARNING: Could not download ${chain}"
    else
        echo "  ${chain}.over.chain.gz exists, skipping..."
    fi
done

cd "${OUTPUT_DIR}"

# =============================================================================
# 2. TOPMED FREEZE 10 REFERENCE (BRAVO) - REQUIRED
# =============================================================================

echo ""
echo "=== 2/5: TOPMed Freeze 10 Reference (BRAVO) ==="
echo ""

cat << 'BRAVO_INSTRUCTIONS'
--------------------------------------------------------------------------------
BRAVO Freeze 10 VCF Download (REQUIRED for imputation QC)
--------------------------------------------------------------------------------

The BRAVO Freeze 10 VCFs require signed URLs that expire after 15 minutes.
You must download them manually:

STEP 1: Visit https://bravo.sph.umich.edu/vcfs.html

STEP 2: Download all chromosome VCFs (chr1-22) to this directory:
        ${OUTPUT_DIR}/bravo_vcfs/

STEP 3: Run the conversion script:
        ./helper_scripts/convert_bravo_to_rayner.sh \
            --input-dir ${OUTPUT_DIR}/bravo_vcfs/ \
            --output ${OUTPUT_DIR}/rayner_refs/PASS.Variants.TOPMed_freeze10_hg38.tab.gz

--------------------------------------------------------------------------------
BRAVO_INSTRUCTIONS

# Create interactive download helper script
cat > bravo_vcfs/download_bravo_freeze10.sh << 'DOWNLOAD_SCRIPT'
#!/bin/bash
################################################################################
# BRAVO Freeze 10 Download Helper
#
# URLs expire after 15 minutes - download quickly!
# Visit: https://bravo.sph.umich.edu/vcfs.html
################################################################################

echo "=============================================="
echo "BRAVO Freeze 10 VCF Download"
echo "=============================================="
echo ""
echo "1. Open https://bravo.sph.umich.edu/vcfs.html in your browser"
echo "2. For each chromosome, right-click the download link and copy URL"
echo "3. Paste the URL when prompted below"
echo ""
echo "URLs expire in 15 minutes - work quickly!"
echo ""

for chr in {1..22}; do
    if [ ! -f "chr${chr}.bravo_freeze10.vcf.gz" ]; then
        echo "----------------------------------------------"
        echo "Chromosome ${chr}: Paste the signed URL and press Enter"
        echo "(or press Enter to skip)"
        read -r url
        if [ -n "$url" ]; then
            echo "Downloading chr${chr}..."
            wget -O "chr${chr}.bravo_freeze10.vcf.gz" "$url" && \
                echo "  Downloaded chr${chr} successfully" || \
                echo "  WARNING: Failed to download chr${chr}"
        else
            echo "  Skipped chr${chr}"
        fi
    else
        echo "chr${chr}.bravo_freeze10.vcf.gz exists, skipping..."
    fi
done

echo ""
echo "=============================================="
echo "Download complete!"
echo "=============================================="
echo ""
echo "Next step: Convert to Rayner format:"
echo "  ../helper_scripts/convert_bravo_to_rayner.sh \\"
echo "      --input-dir ./ \\"
echo "      --output ../rayner_refs/PASS.Variants.TOPMed_freeze10_hg38.tab.gz"
echo ""
DOWNLOAD_SCRIPT
chmod +x bravo_vcfs/download_bravo_freeze10.sh

echo "  Created: bravo_vcfs/download_bravo_freeze10.sh"
echo "  Run this script after visiting https://bravo.sph.umich.edu/vcfs.html"

# =============================================================================
# 3. HRC REFERENCE (hg19/GRCh37)
# =============================================================================

echo ""
echo "=== 3/5: HRC Reference (hg19/GRCh37) ==="
echo ""

cd rayner_refs

if [ ! -f "HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz" ]; then
    echo "  Downloading HRC r1.1 reference..."
    wget -q -c "https://www.chg.ox.ac.uk/~wrayner/tools/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz" && \
        echo "  Downloaded HRC reference successfully" || \
        echo "  WARNING: Could not download HRC reference"
else
    echo "  HRC reference exists, skipping..."
fi

cd "${OUTPUT_DIR}"

# =============================================================================
# 4. REFERENCE FASTA FILES
# =============================================================================

echo ""
echo "=== 4/5: Reference FASTA Files ==="
echo ""

if [ "$SKIP_FASTA" = "true" ]; then
    echo "  Skipping FASTA download (--skip-fasta flag set)"
else
    cd fasta

    # GRCh37/hg19
    mkdir -p hg19
    cd hg19
    if [ ! -f "human_g1k_v37.fasta" ]; then
        echo "  Downloading GRCh37 reference (~3GB)..."
        wget -q -c "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz" && \
            gunzip -k human_g1k_v37.fasta.gz 2>/dev/null && \
            echo "  Downloaded GRCh37 reference" || \
            echo "  WARNING: Could not download GRCh37 reference"

        if command -v samtools &> /dev/null && [ -f human_g1k_v37.fasta ]; then
            echo "  Creating FASTA index..."
            samtools faidx human_g1k_v37.fasta
        fi
    else
        echo "  GRCh37 reference exists, skipping..."
    fi
    cd ..

    # GRCh38/hg38
    mkdir -p hg38
    cd hg38
    if [ ! -f "GRCh38_full_analysis_set_plus_decoy_hla.fa" ]; then
        echo "  Downloading GRCh38 reference (~3GB)..."
        wget -q -c "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa" && \
            echo "  Downloaded GRCh38 reference" || \
            echo "  WARNING: Could not download GRCh38 reference"

        if command -v samtools &> /dev/null && [ -f GRCh38_full_analysis_set_plus_decoy_hla.fa ]; then
            echo "  Creating FASTA index..."
            samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa
        fi
    else
        echo "  GRCh38 reference exists, skipping..."
    fi

    cd "${OUTPUT_DIR}"
fi

# =============================================================================
# 5. GENETIC MAPS (Public - no custom panel needed!)
# =============================================================================

echo ""
echo "=== 5/5: Genetic Maps (SHAPEIT5 + TractorWorkflow) ==="
echo ""

if [ "$SKIP_GENETIC_MAPS" = "true" ]; then
    echo "  Skipping genetic maps download (--skip-genetic-maps flag set)"
else
    cd genetic_maps

    # Create subdirectories
    mkdir -p shapeit5_b37 shapeit5_b38 tractor_format

    # Download SHAPEIT5 genetic maps (b37/hg19)
    if [ ! -f "shapeit5_b37/chr1.b37.gmap.gz" ]; then
        echo "  Downloading SHAPEIT5 genetic maps (b37/hg19)..."
        wget -q -c "https://github.com/odelaneau/shapeit5/raw/main/maps/genetic_maps.b37.tar.gz" -O shapeit5_b37.tar.gz && \
            tar -xzf shapeit5_b37.tar.gz -C shapeit5_b37 --strip-components=1 2>/dev/null && \
            rm -f shapeit5_b37.tar.gz && \
            echo "  Downloaded b37 genetic maps" || \
            echo "  WARNING: Could not download b37 maps"
    else
        echo "  SHAPEIT5 b37 maps exist, skipping..."
    fi

    # Download SHAPEIT5 genetic maps (b38/hg38)
    if [ ! -f "shapeit5_b38/chr1.b38.gmap.gz" ]; then
        echo "  Downloading SHAPEIT5 genetic maps (b38/hg38)..."
        wget -q -c "https://github.com/odelaneau/shapeit5/raw/main/maps/genetic_maps.b38.tar.gz" -O shapeit5_b38.tar.gz && \
            tar -xzf shapeit5_b38.tar.gz -C shapeit5_b38 --strip-components=1 2>/dev/null && \
            rm -f shapeit5_b38.tar.gz && \
            echo "  Downloaded b38 genetic maps" || \
            echo "  WARNING: Could not download b38 maps"
    else
        echo "  SHAPEIT5 b38 maps exist, skipping..."
    fi

    # Download TractorWorkflow genetic maps (RFMix2/GNomix format)
    if [ ! -d "tractor_format/references" ]; then
        echo "  Downloading TractorWorkflow genetic maps..."
        wget -q -c "$TRACTOR_TESTDATA" -O tractor_test_data.zip && \
            unzip -q -o tractor_test_data.zip && \
            mv test_data/references tractor_format/ && \
            rm -rf test_data tractor_test_data.zip && \
            echo "  Downloaded TractorWorkflow maps" || \
            echo "  WARNING: Could not download TractorWorkflow data"
    else
        echo "  TractorWorkflow maps exist, skipping..."
    fi

    cd "${OUTPUT_DIR}"
fi

# =============================================================================
# CREATE SUMMARY
# =============================================================================

echo ""
echo "=== Creating Summary ==="

cat > REFERENCES_README.txt << 'EOF'
================================================================================
Pipeline Reference Files (Modules 1-7)
================================================================================

ALL files in this directory are PUBLICLY AVAILABLE.
NO custom reference panel is needed for these downloads.

Directory Structure:
--------------------

chain_files/
├── hg19ToHg38.over.chain.gz      # Liftover hg19 → hg38
└── hg38ToHg19.over.chain.gz      # Liftover hg38 → hg19

bravo_vcfs/
├── download_bravo_freeze10.sh    # Helper script for BRAVO download
└── chr*.bravo_freeze10.vcf.gz    # BRAVO Freeze 10 VCFs (manual download)

rayner_refs/
├── PASS.Variants.TOPMed_freeze10_hg38.tab.gz  # TOPMed Fr10 (from BRAVO)
└── HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz      # HRC for hg19

fasta/
├── hg19/human_g1k_v37.fasta      # GRCh37 reference
└── hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa  # GRCh38 reference

genetic_maps/
├── shapeit5_b37/                 # SHAPEIT5 genetic maps (hg19)
├── shapeit5_b38/                 # SHAPEIT5 genetic maps (hg38)
└── tractor_format/               # RFMix2/GNomix format maps


BRAVO Freeze 10 Setup (REQUIRED):
---------------------------------

1. Run: cd bravo_vcfs && ./download_bravo_freeze10.sh
2. Follow prompts to download each chromosome
3. Convert: ./helper_scripts/convert_bravo_to_rayner.sh ...


What Still Requires YOUR Custom Data:
-------------------------------------

ONLY the LAI reference panel (for RFMix, FLARE, G-NOMIX) requires your
HGDP-1KG + MX Biobank WGS data. Prepare it using:

    ./helper_scripts/format_reference_panel.sh \
        --input /path/to/your_combined_reference.vcf.gz \
        --sample-map /path/to/sample_populations.txt \
        --genetic-map-dir genetic_maps/shapeit5_b38/ \
        --output-dir lai_reference/


Copy to Pipeline:
-----------------

# All core references
cp chain_files/* resources/chain_files/
cp rayner_refs/* resources/rayner/
cp -r fasta/* resources/references/
cp -r genetic_maps/* resources/genetic_maps/

# LAI reference (after you prepare it with your custom data)
cp -r lai_reference/* resources/ancestry_references/lai_reference/

================================================================================
EOF

echo ""
echo "=============================================="
echo "Reference Download Complete!"
echo "=============================================="
echo ""
echo "Downloaded files:"
find "${OUTPUT_DIR}" -type f \( -name "*.gz" -o -name "*.fa" -o -name "*.fasta" \) 2>/dev/null | head -15
echo ""
echo "Total size:"
du -sh "${OUTPUT_DIR}"
echo ""
echo "=============================================="
echo "IMPORTANT: BRAVO Freeze 10 (TOPMed) still needs manual download!"
echo "=============================================="
echo ""
echo "Run: cd ${OUTPUT_DIR}/bravo_vcfs && ./download_bravo_freeze10.sh"
echo ""
echo "After downloading BRAVO, convert it:"
echo "  ./helper_scripts/convert_bravo_to_rayner.sh \\"
echo "      --input-dir ${OUTPUT_DIR}/bravo_vcfs/ \\"
echo "      --output ${OUTPUT_DIR}/rayner_refs/PASS.Variants.TOPMed_freeze10_hg38.tab.gz"
echo ""
