#!/bin/bash
################################################################################
# setup_reference_files.sh
#
# Wrapper script to set up all reference files for the genotyping pipeline.
#
# Downloads ALL publicly available references (Modules 1-7):
#   ./download_core_references.sh
#
# For LAI reference panel (requires YOUR custom HGDP-1KG + MX Biobank data):
#   ./format_reference_panel.sh
#
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${OUTPUT_DIR:-./pipeline_references}"

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Downloads all publicly available reference files for the genotyping pipeline.

Options:
  -o, --output-dir DIR   Output directory (default: ./pipeline_references)
  --skip-fasta           Skip FASTA download
  --skip-genetic-maps    Skip genetic maps download
  -h, --help             Show this help

What Gets Downloaded (ALL PUBLIC - no custom data needed):
----------------------------------------------------------

For Modules 1-6 (Imputation):
  - Chain files (hg19 <-> hg38 liftover)
  - TOPMed Freeze 10 reference (BRAVO) + helper script
  - HRC reference (hg19)
  - Reference FASTA files

For Module 7 (Ancestry):
  - SHAPEIT5 genetic maps (b37 and b38)
  - TractorWorkflow genetic maps

What Requires YOUR Custom Data (not downloaded):
------------------------------------------------

  LAI Reference Panel (for RFMix, FLARE, G-NOMIX)
  -> Requires HGDP-1KG + MX Biobank WGS data
  -> Prepare using: ./format_reference_panel.sh

Workflow:
---------

1. Download public references:
   ./setup_reference_files.sh -o /path/to/references

2. Download BRAVO Freeze 10 (manual - signed URLs):
   cd /path/to/references/bravo_vcfs
   ./download_bravo_freeze10.sh

3. Convert BRAVO to Rayner format:
   ./convert_bravo_to_rayner.sh ...

4. (When ready) Prepare LAI reference panel:
   ./format_reference_panel.sh --input your_reference.vcf.gz ...

EOF
}

# Pass all arguments to download_core_references.sh
ARGS=()
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            print_usage
            exit 0
            ;;
        *)
            ARGS+=("$1")
            shift
            ;;
    esac
done

echo "=============================================="
echo "Genotyping Pipeline Reference Setup"
echo "=============================================="
echo ""

# Run the core references download
"${SCRIPT_DIR}/download_core_references.sh" "${ARGS[@]}"

echo ""
echo "=============================================="
echo "Setup Complete!"
echo "=============================================="
echo ""
echo "All public references have been downloaded."
echo ""
echo "Remaining manual steps:"
echo "  1. Download BRAVO Freeze 10 (see bravo_vcfs/download_bravo_freeze10.sh)"
echo "  2. Convert BRAVO to Rayner format"
echo "  3. (Optional) Prepare LAI reference panel with format_reference_panel.sh"
echo ""
