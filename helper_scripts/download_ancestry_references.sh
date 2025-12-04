#!/bin/bash
################################################################################
# download_ancestry_references.sh
#
# Downloads and organizes reference files for Module 7 (Ancestry Estimation):
#   - Genetic maps (SHAPEIT5 format for LAI tools)
#   - TractorWorkflow reference data
#   - Instructions for HGDP-1KG + MX Biobank reference panel preparation
#
# IMPORTANT: This script downloads publicly available genetic maps, but the
# LAI reference panel (HGDP-1KG + MX Biobank) must be prepared separately
# using format_reference_panel.sh with YOUR WGS reference data.
#
# Usage:
#   ./download_ancestry_references.sh --output-dir /path/to/references
#
# Requirements:
#   - wget, unzip
#   - ~10GB disk space for genetic maps
#   - User-provided HGDP-1KG + MX Biobank reference panel
#
################################################################################

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

OUTPUT_DIR="${OUTPUT_DIR:-./pipeline_references}"

# TractorWorkflow reference data
TRACTOR_TESTDATA="https://github.com/Atkinson-Lab/Tractor-tutorial/raw/refs/heads/main/test_data.zip"

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Downloads reference files for Module 7 (Ancestry Estimation).

Options:
  -o, --output-dir DIR   Output directory (default: ./pipeline_references)
  -h, --help             Show this help

What This Script Downloads:
  [x] SHAPEIT5 genetic maps (b37 and b38)
  [x] TractorWorkflow reference data (RFMix2/GNomix format maps)

What You Must Provide (NOT Downloaded):
  [ ] HGDP-1KG high-coverage WGS reference panel
  [ ] MX Biobank samples (or other Latino reference samples)
  [ ] Sample-to-population mapping file

After downloading, use format_reference_panel.sh to prepare your LAI reference.

Reference Sources:
  - Genetic maps:       SHAPEIT5 / TractorWorkflow
  - HGDP-1KG WGS:       NYGC (https://www.internationalgenome.org/data-portal/)
  - MX Biobank:         User-provided / collaborator data

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output-dir)
            OUTPUT_DIR="$2"
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
echo "Ancestry Reference Files (Module 7)"
echo "=============================================="
echo "Output directory: ${OUTPUT_DIR}"
echo ""

mkdir -p "${OUTPUT_DIR}"/{genetic_maps,lai_reference,tractor_refs}

cd "${OUTPUT_DIR}"

# =============================================================================
# 1. GENETIC MAPS (SHAPEIT5)
# =============================================================================

echo ""
echo "=== 1/3: Downloading Genetic Maps ==="
echo ""

cd genetic_maps

# Create subdirectories
mkdir -p shapeit5_b37 shapeit5_b38

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

cd "${OUTPUT_DIR}"

# =============================================================================
# 2. TRACTORWORKFLOW REFERENCE DATA
# =============================================================================

echo ""
echo "=== 2/3: Downloading TractorWorkflow Reference Data ==="
echo ""

cd tractor_refs

if [ ! -d "references" ]; then
    echo "  Downloading TractorWorkflow test data..."
    wget -q -c "$TRACTOR_TESTDATA" -O tractor_test_data.zip && \
        unzip -q -o tractor_test_data.zip && \
        mv test_data/references . && \
        rm -rf test_data tractor_test_data.zip && \
        echo "  Downloaded TractorWorkflow references" || \
        echo "  WARNING: Could not download TractorWorkflow data"
else
    echo "  TractorWorkflow references exist, skipping..."
fi

cd "${OUTPUT_DIR}"

# =============================================================================
# 3. LAI REFERENCE PANEL INSTRUCTIONS
# =============================================================================

echo ""
echo "=== 3/3: LAI Reference Panel Instructions ==="
echo ""

cat << 'LAI_INSTRUCTIONS'
--------------------------------------------------------------------------------
HGDP-1KG + MX Biobank Reference Panel Setup
--------------------------------------------------------------------------------

The Local Ancestry Inference (LAI) reference panel is NOT automatically
downloaded because it requires user-provided data (MX Biobank or other
Latino reference samples).

Required Components:
--------------------

1. HGDP-1KG High-Coverage WGS (~4000 samples)
   - Download from: https://www.internationalgenome.org/data-portal/
   - Files: phased VCF files per chromosome

2. MX Biobank or Latino Reference Samples
   - User-provided WGS data
   - Must be phased (or will be phased during processing)

3. Sample-to-Population Mapping File
   - Tab-separated: SampleID <TAB> Population <TAB> Superpopulation
   - Example:
       NA19700    YRI    AFR
       NA12878    CEU    EUR
       MXB00001   MXL    AMR

Preparation Steps:
------------------

1. Combine your reference samples into a single phased VCF:

   bcftools merge --force-samples \
       hgdp_1kg_phased.vcf.gz \
       mxbiobank_phased.vcf.gz \
       -Oz -o combined_reference.vcf.gz

2. Create sample-population mapping file

3. Run the formatting script:

   ./helper_scripts/format_reference_panel.sh \
       --input combined_reference.vcf.gz \
       --sample-map sample_populations.txt \
       --genetic-map-dir genetic_maps/shapeit5_b38/ \
       --output-dir lai_reference/

This will create:
  - lai_reference/           Shared VCF + sample maps for RFMix v2, FLARE, G-NOMIX
  - rfmix1/                  Legacy RFMix v1 format
  - admixture/               LD-pruned PLINK for ADMIXTURE

--------------------------------------------------------------------------------
LAI_INSTRUCTIONS

# Create placeholder files to remind user
mkdir -p lai_reference
cat > lai_reference/README.txt << 'README'
================================================================================
LAI Reference Panel (User-Provided)
================================================================================

This directory should contain your formatted LAI reference panel.

Required files after running format_reference_panel.sh:
  - reference_haplotypes.vcf.gz      # Phased reference VCF
  - reference_haplotypes.vcf.gz.tbi  # Index
  - rfmix2_sample_map.txt            # RFMix v2 sample map
  - flare_panels.txt                 # FLARE panel file
  - gnomix_sample_map.tsv            # G-NOMIX sample map

To prepare these files, run:

  ./helper_scripts/format_reference_panel.sh \
      --input /path/to/your_combined_reference.vcf.gz \
      --sample-map /path/to/sample_populations.txt \
      --genetic-map-dir ../genetic_maps/shapeit5_b38/ \
      --output-dir ./

See download_ancestry_references.sh for detailed instructions.

================================================================================
README

# =============================================================================
# CREATE SUMMARY
# =============================================================================

echo ""
echo "=== Creating Summary ==="

cat > ANCESTRY_REFERENCES_README.txt << 'EOF'
================================================================================
Ancestry Reference Files (Module 7)
================================================================================

This directory contains reference files for Module 7 (Ancestry Estimation).

Directory Structure:
--------------------

genetic_maps/
├── shapeit5_b37/                 # SHAPEIT5 genetic maps (hg19)
│   └── chr*.b37.gmap.gz
├── shapeit5_b38/                 # SHAPEIT5 genetic maps (hg38)
│   └── chr*.b38.gmap.gz
└── tractor_maps/                 # RFMix2/GNomix format (from TractorWorkflow)
    ├── chr*.genetic_map.modified.txt
    └── ...

tractor_refs/
└── references/                   # TractorWorkflow test reference data
    ├── genetic_maps/
    └── sample_maps/

lai_reference/                    # YOUR REFERENCE PANEL (user-provided)
├── reference_haplotypes.vcf.gz   # Combined HGDP-1KG + MX Biobank
├── rfmix2_sample_map.txt         # RFMix v2 format
├── flare_panels.txt              # FLARE format
└── gnomix_sample_map.tsv         # G-NOMIX format


What's Downloaded vs User-Provided:
-----------------------------------

Downloaded (ready to use):
  [x] SHAPEIT5 genetic maps (b37 and b38)
  [x] TractorWorkflow reference data

User-Provided (you must prepare):
  [ ] HGDP-1KG high-coverage WGS
  [ ] MX Biobank / Latino reference samples
  [ ] Combined + formatted LAI reference panel


Copy to Pipeline:
-----------------

# Genetic maps (downloaded)
cp -r genetic_maps/* /path/to/pipeline/resources/genetic_maps/

# LAI reference (after you prepare it)
cp -r lai_reference/* /path/to/pipeline/resources/ancestry_references/lai_reference/


Module 7 Tools & Their Reference Needs:
---------------------------------------

GRAF-anc:     Built-in reference (no external files needed)
ADMIXTURE:    LD-pruned PLINK from lai_reference/ (format_reference_panel.sh)
RFMix v2:     lai_reference/ + genetic_maps/
FLARE:        lai_reference/ + genetic_maps/
G-NOMIX:      lai_reference/ + genetic_maps/
RFMix v1:     rfmix1/ format (legacy, format_reference_panel.sh)

================================================================================
EOF

echo ""
echo "=============================================="
echo "Ancestry Reference Download Complete!"
echo "=============================================="
echo ""
echo "Downloaded files:"
find "${OUTPUT_DIR}/genetic_maps" -type f -name "*.gz" 2>/dev/null | head -10
echo ""
echo "Total size:"
du -sh "${OUTPUT_DIR}/genetic_maps" "${OUTPUT_DIR}/tractor_refs" 2>/dev/null || true
echo ""
echo "IMPORTANT: LAI reference panel (HGDP-1KG + MX Biobank) NOT included!"
echo ""
echo "Next steps:"
echo "  1. Obtain HGDP-1KG phased WGS from NYGC"
echo "  2. Combine with your MX Biobank / Latino reference samples"
echo "  3. Create sample-population mapping file"
echo "  4. Run format_reference_panel.sh to prepare LAI reference"
echo ""
echo "See ANCESTRY_REFERENCES_README.txt for detailed instructions."
echo ""
