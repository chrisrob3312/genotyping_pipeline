#!/bin/bash
################################################################################
# setup_reference_files.sh
#
# WRAPPER SCRIPT - Calls the appropriate reference download script based on
# which modules you need to run.
#
# This script has been split into two parts for flexibility:
#
#   1. download_core_references.sh     - Modules 1-6 (can run immediately)
#   2. download_ancestry_references.sh - Module 7 (requires LAI reference panel)
#
# Usage:
#   ./setup_reference_files.sh [--all | --core | --ancestry] [OPTIONS]
#
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${OUTPUT_DIR:-./pipeline_references}"
MODE="all"

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [MODE] [OPTIONS]

Wrapper script to download reference files for the genotyping pipeline.

Modes:
  --all        Download both core and ancestry references (default)
  --core       Download only core references (Modules 1-6)
  --ancestry   Download only ancestry references (Module 7)

Options:
  -o, --output-dir DIR   Output directory (default: ./pipeline_references)
  -h, --help             Show this help

Reference Scripts:
------------------

Core References (Modules 1-6):
  ./download_core_references.sh

  Downloads:
    - Chain files (hg19 <-> hg38 liftover)
    - TOPMed Freeze 5 reference (hg38)
    - HRC reference (hg19)
    - Reference FASTA files

  These are publicly available and can be downloaded immediately.

Ancestry References (Module 7):
  ./download_ancestry_references.sh

  Downloads:
    - SHAPEIT5 genetic maps
    - TractorWorkflow reference data

  IMPORTANT: The LAI reference panel (HGDP-1KG + MX Biobank) is NOT
  downloaded automatically. You must prepare it using:

    ./format_reference_panel.sh

Recommended Workflow:
---------------------

1. Run core references first:
   ./setup_reference_files.sh --core -o /path/to/references

2. Test Modules 1-6 with the core references

3. When ready for Module 7, download ancestry references:
   ./setup_reference_files.sh --ancestry -o /path/to/references

4. Prepare your LAI reference panel:
   ./format_reference_panel.sh --input your_reference.vcf.gz ...

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --all)
            MODE="all"
            shift
            ;;
        --core)
            MODE="core"
            shift
            ;;
        --ancestry)
            MODE="ancestry"
            shift
            ;;
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
echo "Genotyping Pipeline Reference Setup"
echo "=============================================="
echo "Mode:             ${MODE}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# =============================================================================
# Run Appropriate Scripts
# =============================================================================

if [[ "$MODE" == "all" ]] || [[ "$MODE" == "core" ]]; then
    echo ""
    echo ">>> Running Core References Download (Modules 1-6) <<<"
    echo ""

    "${SCRIPT_DIR}/download_core_references.sh" --output-dir "${OUTPUT_DIR}"
fi

if [[ "$MODE" == "all" ]] || [[ "$MODE" == "ancestry" ]]; then
    echo ""
    echo ">>> Running Ancestry References Download (Module 7) <<<"
    echo ""

    "${SCRIPT_DIR}/download_ancestry_references.sh" --output-dir "${OUTPUT_DIR}"
fi

# =============================================================================
# Summary
# =============================================================================

echo ""
echo "=============================================="
echo "Reference Setup Complete!"
echo "=============================================="
echo ""
echo "Downloaded references are in: ${OUTPUT_DIR}"
echo ""

if [[ "$MODE" == "all" ]] || [[ "$MODE" == "ancestry" ]]; then
    echo "REMINDER: For Module 7, you still need to prepare your LAI reference panel."
    echo "Run: ./format_reference_panel.sh --help for instructions."
    echo ""
fi

echo "To copy references to your pipeline:"
echo ""
echo "  # Core references (Modules 1-6)"
echo "  cp ${OUTPUT_DIR}/chain_files/* resources/chain_files/"
echo "  cp ${OUTPUT_DIR}/rayner_refs/* resources/rayner/"
echo "  cp -r ${OUTPUT_DIR}/fasta/* resources/references/"
echo ""

if [[ "$MODE" == "all" ]] || [[ "$MODE" == "ancestry" ]]; then
    echo "  # Ancestry references (Module 7)"
    echo "  cp -r ${OUTPUT_DIR}/genetic_maps/* resources/genetic_maps/"
    echo "  cp -r ${OUTPUT_DIR}/lai_reference/* resources/ancestry_references/lai_reference/"
    echo ""
fi
