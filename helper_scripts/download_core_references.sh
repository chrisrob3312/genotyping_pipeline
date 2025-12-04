#!/bin/bash
################################################################################
# download_core_references.sh
#
# Downloads and organizes CORE reference files for Modules 1-6 (Pre-Ancestry):
#   - Chain files (liftover between hg19/hg38)
#   - TOPMed Freeze 5/10 reference (Rayner format for imputation QC)
#   - HRC reference for hg19 fallback
#   - Reference genome FASTA files (hg19/hg38)
#
# These references are publicly available and can be downloaded immediately.
# Run this script BEFORE running Modules 1-6.
#
# For Module 7 (Ancestry), use download_ancestry_references.sh instead
# (requires HGDP-1KG + MX Biobank reference panel to be prepared first).
#
# Usage:
#   ./download_core_references.sh --output-dir /path/to/references
#
# Requirements:
#   - wget, curl
#   - samtools (for FASTA indexing, optional)
#   - ~50GB disk space
#
################################################################################

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

OUTPUT_DIR="${OUTPUT_DIR:-./pipeline_references}"
SKIP_BRAVO="${SKIP_BRAVO:-false}"
SKIP_FASTA="${SKIP_FASTA:-false}"

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Downloads CORE reference files for Modules 1-6 (Pre-Ancestry).

Options:
  -o, --output-dir DIR   Output directory (default: ./pipeline_references)
  --skip-bravo           Skip BRAVO/TOPMed download instructions
  --skip-fasta           Skip FASTA download (if you already have them)
  -h, --help             Show this help

Reference Sources:
  - Chain files:    UCSC Genome Browser
  - TOPMed (hg38):  BRAVO Freeze 10 or Freeze 5 fallback
  - HRC (hg19):     Wellcome Sanger / Will Rayner's site
  - FASTA:          1000 Genomes / NCBI

What This Script Downloads:
  [x] Chain files (hg19 <-> hg38 liftover)
  [x] TOPMed Freeze 5 reference (fallback, auto-download)
  [x] HRC r1.1 reference (hg19)
  [x] Reference FASTA files (hg19, hg38)
  [x] BRAVO Freeze 10 download helper (manual step required)

For Module 7 (Ancestry), run download_ancestry_references.sh separately.

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --skip-bravo)
            SKIP_BRAVO=true
            shift
            ;;
        --skip-fasta)
            SKIP_FASTA=true
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
echo "Core Reference Files (Modules 1-6)"
echo "=============================================="
echo "Output directory: ${OUTPUT_DIR}"
echo ""

mkdir -p "${OUTPUT_DIR}"/{chain_files,bravo_vcfs,rayner_refs,fasta}

cd "${OUTPUT_DIR}"

# =============================================================================
# 1. CHAIN FILES (Liftover)
# =============================================================================

echo ""
echo "=== 1/4: Downloading Chain Files ==="
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
# 2. TOPMED REFERENCE (hg38)
# =============================================================================

echo ""
echo "=== 2/4: TOPMed Reference (hg38) ==="
echo ""

cd rayner_refs

# Download Freeze 5 (automatic, always available)
if [ ! -f "PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz" ]; then
    echo "  Downloading TOPMed Freeze 5 reference..."
    wget -q -c "https://www.chg.ox.ac.uk/~wrayner/tools/PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz" && \
        echo "  Downloaded TOPMed Freeze 5 successfully" || \
        echo "  WARNING: Could not download Freeze 5 reference"
else
    echo "  TOPMed Freeze 5 exists, skipping..."
fi

cd "${OUTPUT_DIR}"

# BRAVO Freeze 10 instructions (manual download required)
if [ "$SKIP_BRAVO" != "true" ]; then
    cat << 'BRAVO_INSTRUCTIONS'

--------------------------------------------------------------------------------
OPTIONAL: BRAVO Freeze 10 (Recommended for Better Coverage)
--------------------------------------------------------------------------------

The BRAVO VCF files require signed URLs that expire after 15 minutes.
You must download them manually from:

    https://bravo.sph.umich.edu/vcfs.html

Steps:
1. Visit the URL above in your browser
2. Download all chromosome VCFs (chr1-22, X) to bravo_vcfs/
3. Run the conversion script:

   ./helper_scripts/convert_bravo_to_rayner.sh \
       --input-dir bravo_vcfs/ \
       --output rayner_refs/PASS.Variants.TOPMed_freeze10_hg38.tab.gz

NOTE: Freeze 5 (already downloaded) is sufficient for most use cases.
      Freeze 10 provides better rare variant coverage.

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

# =============================================================================
# 3. HRC REFERENCE (hg19/GRCh37)
# =============================================================================

echo ""
echo "=== 3/4: HRC Reference (hg19/GRCh37) ==="
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
echo "=== 4/4: Reference FASTA Files ==="
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
            gunzip -k human_g1k_v37.fasta.gz && \
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
# CREATE SUMMARY
# =============================================================================

echo ""
echo "=== Creating Summary ==="

cat > CORE_REFERENCES_README.txt << 'EOF'
================================================================================
Core Pipeline Reference Files (Modules 1-6)
================================================================================

This directory contains reference files for Modules 1-6 (Pre-Ancestry).

Directory Structure:
--------------------

chain_files/
├── hg19ToHg38.over.chain.gz      # Liftover hg19 → hg38
└── hg38ToHg19.over.chain.gz      # Liftover hg38 → hg19

rayner_refs/
├── PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz   # TOPMed Fr5 (default)
├── PASS.Variants.TOPMed_freeze10_hg38.tab.gz        # TOPMed Fr10 (optional)
└── HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz            # HRC for hg19

fasta/
├── hg19/human_g1k_v37.fasta      # GRCh37 reference
└── hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa  # GRCh38 reference


Copy to Pipeline:
-----------------

cp chain_files/* /path/to/pipeline/resources/chain_files/
cp rayner_refs/* /path/to/pipeline/resources/rayner/
cp -r fasta/* /path/to/pipeline/resources/references/


Modules Supported:
------------------

  Module 1 (Pre-Imputation QC):  chain_files/, rayner_refs/, fasta/
  Module 2 (Imputation):         (uses imputation server references)
  Module 3 (Post-Imputation QC): rayner_refs/
  Module 4 (Platform Merging):   fasta/
  Module 5 (Re-Imputation):      (uses imputation server references)
  Module 6 (Post-Merge QC):      rayner_refs/

For Module 7 (Ancestry), run download_ancestry_references.sh separately.

================================================================================
EOF

echo ""
echo "=============================================="
echo "Core Reference Download Complete!"
echo "=============================================="
echo ""
echo "Downloaded files:"
find "${OUTPUT_DIR}" -type f \( -name "*.gz" -o -name "*.fa" -o -name "*.fasta" \) 2>/dev/null | head -10
echo ""
echo "Total size:"
du -sh "${OUTPUT_DIR}"
echo ""
echo "These references are ready for Modules 1-6."
echo ""
echo "Next steps:"
echo "  1. (Optional) Download BRAVO Freeze 10 for better coverage"
echo "  2. Copy files to pipeline resources directory"
echo "  3. For Module 7, run download_ancestry_references.sh"
echo ""
