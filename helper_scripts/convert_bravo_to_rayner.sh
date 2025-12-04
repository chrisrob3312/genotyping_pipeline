#!/bin/bash
################################################################################
# convert_bravo_to_rayner.sh
#
# Converts BRAVO Freeze 10 VCF files to Rayner/McCarthy tab format for use with
# the HRC-1000G-check-bim.pl script.
#
# Input:  BRAVO VCF files (chr1-22, X) from https://bravo.sph.umich.edu/vcfs.html
# Output: PASS.Variants.TOPMed_freeze10_hg38.tab.gz (Rayner format)
#
# Rayner tab format columns:
#   1. CHR      - Chromosome (1-22, X)
#   2. POS      - Position (1-based)
#   3. ID       - Variant ID (rsID or chr:pos:ref:alt)
#   4. REF      - Reference allele
#   5. ALT      - Alternate allele
#   6. AF       - Allele frequency
#
# Usage:
#   ./convert_bravo_to_rayner.sh \
#       --input-dir /path/to/bravo_vcfs/ \
#       --output /path/to/PASS.Variants.TOPMed_freeze10_hg38.tab.gz
#
# Requirements:
#   - bcftools
#   - bgzip, tabix
#   - ~50GB disk space for intermediate files
#
################################################################################

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

INPUT_DIR=""
OUTPUT_FILE=""
THREADS=4
KEEP_TEMP=false
CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Converts BRAVO Freeze 10 VCFs to Rayner tab format.

Required:
  -i, --input-dir DIR     Directory containing BRAVO VCF files (chr*.bravo.pub.vcf.gz)
  -o, --output FILE       Output file path (will be gzipped)

Optional:
  -t, --threads N         Number of threads (default: 4)
  -c, --chromosomes LIST  Chromosomes to process (default: 1-22 X)
  --keep-temp             Keep intermediate files
  -h, --help              Show this help

Input files expected:
  chr1.bravo.pub.vcf.gz
  chr2.bravo.pub.vcf.gz
  ...
  chr22.bravo.pub.vcf.gz
  chrX.bravo.pub.vcf.gz

Output format (Rayner tab):
  CHR  POS  ID  REF  ALT  AF

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -c|--chromosomes)
            CHROMOSOMES="$2"
            shift 2
            ;;
        --keep-temp)
            KEEP_TEMP=true
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

# Validate arguments
if [[ -z "$INPUT_DIR" ]] || [[ -z "$OUTPUT_FILE" ]]; then
    echo "ERROR: --input-dir and --output are required"
    print_usage
    exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

# Check for required tools
for tool in bcftools bgzip; do
    if ! command -v $tool &> /dev/null; then
        echo "ERROR: Required tool not found: $tool"
        exit 1
    fi
done

echo "=============================================="
echo "BRAVO to Rayner Format Converter"
echo "=============================================="
echo "Input directory: $INPUT_DIR"
echo "Output file: $OUTPUT_FILE"
echo "Threads: $THREADS"
echo "Chromosomes: $CHROMOSOMES"
echo ""

# Create output directory
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

# Create temp directory
TEMP_DIR=$(mktemp -d)
if [[ "$KEEP_TEMP" != "true" ]]; then
    trap "rm -rf $TEMP_DIR" EXIT
fi

# =============================================================================
# Process Each Chromosome
# =============================================================================

echo "=== Processing BRAVO VCFs ==="
echo ""

# Create header
echo -e "#CHR\tPOS\tID\tREF\tALT\tAF" > "${TEMP_DIR}/header.txt"

for chr in $CHROMOSOMES; do
    # Find VCF file (handle both chr1 and 1 naming)
    VCF=""
    for pattern in "chr${chr}.bravo.pub.vcf.gz" "${chr}.bravo.pub.vcf.gz" "chr${chr}.vcf.gz"; do
        if [[ -f "${INPUT_DIR}/${pattern}" ]]; then
            VCF="${INPUT_DIR}/${pattern}"
            break
        fi
    done

    if [[ -z "$VCF" ]]; then
        echo "WARNING: VCF not found for chromosome ${chr}, skipping..."
        continue
    fi

    echo "Processing chromosome ${chr}..."

    # Extract required fields from VCF
    # BRAVO VCF INFO field contains AN (allele number) and AC (allele count)
    # AF = AC / AN
    bcftools query \
        -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\n' \
        "$VCF" \
        2>/dev/null | \
    awk -F'\t' -v OFS='\t' '
    {
        # Clean chromosome name (remove chr prefix if present)
        chr = $1
        gsub(/^chr/, "", chr)

        # Handle missing AF
        af = $6
        if (af == "." || af == "") af = "0"

        # Handle multi-allelic (take first AF if comma-separated)
        split(af, afs, ",")
        af = afs[1]

        # Print in Rayner format
        print chr, $2, $3, $4, $5, af
    }' > "${TEMP_DIR}/chr${chr}.tab"

    # Count variants
    n_vars=$(wc -l < "${TEMP_DIR}/chr${chr}.tab")
    echo "  Extracted ${n_vars} variants from chr${chr}"

done

# =============================================================================
# Combine All Chromosomes
# =============================================================================

echo ""
echo "=== Combining all chromosomes ==="

# Combine all chromosome files
cat "${TEMP_DIR}/header.txt" > "${TEMP_DIR}/combined.tab"

for chr in $CHROMOSOMES; do
    if [[ -f "${TEMP_DIR}/chr${chr}.tab" ]]; then
        cat "${TEMP_DIR}/chr${chr}.tab" >> "${TEMP_DIR}/combined.tab"
    fi
done

# Count total variants
TOTAL_VARS=$(tail -n +2 "${TEMP_DIR}/combined.tab" | wc -l)
echo "Total variants: ${TOTAL_VARS}"

# =============================================================================
# Compress and Finalize
# =============================================================================

echo ""
echo "=== Compressing output ==="

# Compress with bgzip
bgzip -c "${TEMP_DIR}/combined.tab" > "$OUTPUT_FILE"

# Create tabix index (optional, for quick lookup)
if command -v tabix &> /dev/null; then
    tabix -s 1 -b 2 -e 2 -S 1 "$OUTPUT_FILE" 2>/dev/null || true
fi

# =============================================================================
# Summary
# =============================================================================

echo ""
echo "=============================================="
echo "Conversion Complete!"
echo "=============================================="
echo ""
echo "Output file: $OUTPUT_FILE"
echo "Total variants: ${TOTAL_VARS}"
echo "File size: $(ls -lh "$OUTPUT_FILE" | awk '{print $5}')"
echo ""
echo "Usage with Rayner script:"
echo "  perl HRC-1000G-check-bim.pl \\"
echo "      -b your_data.bim \\"
echo "      -f your_data.frq \\"
echo "      -r $OUTPUT_FILE \\"
echo "      -h"
echo ""

if [[ "$KEEP_TEMP" == "true" ]]; then
    echo "Temporary files kept in: $TEMP_DIR"
fi
