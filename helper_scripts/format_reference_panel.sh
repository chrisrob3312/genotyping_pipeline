#!/bin/bash
################################################################################
# format_reference_panel.sh
#
# Converts a phased WGS reference panel (PLINK or VCF) to all formats required
# by the genotyping pipeline's LAI and ancestry tools:
#
#   - RFMix v2:   Phased VCF + sample map (tab-separated)
#   - FLARE:      Phased VCF + panel file (space-separated)
#   - RFMix v1:   .alleles, .classes, .snp_locations per chromosome
#   - G-NOMIX:    Phased VCF + sample map (for training)
#   - ADMIXTURE:  LD-pruned PLINK binary files
#   - GRAF-anc:   Uses built-in reference (no conversion needed)
#
# Usage:
#   ./format_reference_panel.sh \
#       --input /path/to/reference.vcf.gz \
#       --sample-map /path/to/sample_populations.txt \
#       --genetic-map-dir /path/to/genetic_maps/ \
#       --output-dir /path/to/formatted_references
#
# Input sample map format (tab-separated):
#   SampleID    Population    Superpopulation(optional)
#   NA19700     YRI           AFR
#   NA12878     CEU           EUR
#
################################################################################

set -euo pipefail

# =============================================================================
# Default Configuration
# =============================================================================

INPUT_FILE=""
SAMPLE_MAP=""
GENETIC_MAP_DIR=""
OUTPUT_DIR="./formatted_references"
THREADS=4
CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
INPUT_TYPE="vcf"  # vcf or plink
LD_PRUNE_WINDOW=200
LD_PRUNE_STEP=50
LD_PRUNE_R2=0.1   # Strict pruning for ADMIXTURE

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Required:
  -i, --input FILE          Input reference panel (VCF.gz or PLINK prefix)
  -s, --sample-map FILE     Sample to population mapping (tab-separated)
  -g, --genetic-map-dir DIR Directory containing genetic maps per chromosome

Optional:
  -o, --output-dir DIR      Output directory (default: ./formatted_references)
  -t, --threads N           Number of threads (default: 4)
  -c, --chromosomes LIST    Chromosomes to process (default: 1-22)
  --plink                   Input is PLINK format (default: VCF)
  -h, --help                Show this help message

Sample Map Format:
  Tab-separated file with columns: SampleID, Population, [Superpopulation]

  Example:
    NA19700    YRI    AFR
    NA12878    CEU    EUR
    HG00403    CHS    EAS

Output Formats Created:
  rfmix2/           RFMix v2 format (VCF + sample_map.txt)
  flare/            FLARE format (VCF + panels.txt)
  rfmix1/           RFMix v1 format (.alleles, .classes, .snp_locations)
  gnomix/           G-NOMIX format (VCF + sample_map.tsv)
  admixture/        ADMIXTURE format (LD-pruned PLINK)

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_FILE="$2"
            shift 2
            ;;
        -s|--sample-map)
            SAMPLE_MAP="$2"
            shift 2
            ;;
        -g|--genetic-map-dir)
            GENETIC_MAP_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
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
        --plink)
            INPUT_TYPE="plink"
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

# Validate required arguments
if [[ -z "$INPUT_FILE" ]] || [[ -z "$SAMPLE_MAP" ]] || [[ -z "$GENETIC_MAP_DIR" ]]; then
    echo "ERROR: Missing required arguments"
    print_usage
    exit 1
fi

# Check input files exist
if [[ ! -f "$SAMPLE_MAP" ]]; then
    echo "ERROR: Sample map not found: $SAMPLE_MAP"
    exit 1
fi

if [[ "$INPUT_TYPE" == "vcf" ]] && [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERROR: Input VCF not found: $INPUT_FILE"
    exit 1
fi

if [[ "$INPUT_TYPE" == "plink" ]] && [[ ! -f "${INPUT_FILE}.bed" ]]; then
    echo "ERROR: Input PLINK files not found: ${INPUT_FILE}.bed/bim/fam"
    exit 1
fi

echo "=============================================="
echo "Reference Panel Formatting Script"
echo "=============================================="
echo "Input file:      $INPUT_FILE"
echo "Input type:      $INPUT_TYPE"
echo "Sample map:      $SAMPLE_MAP"
echo "Genetic maps:    $GENETIC_MAP_DIR"
echo "Output dir:      $OUTPUT_DIR"
echo "Threads:         $THREADS"
echo "Chromosomes:     $CHROMOSOMES"
echo "=============================================="

# =============================================================================
# Create Output Directories
# =============================================================================

mkdir -p "${OUTPUT_DIR}"/{rfmix2,flare,rfmix1,gnomix,admixture,logs}

LOGDIR="${OUTPUT_DIR}/logs"

# =============================================================================
# Helper Functions
# =============================================================================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${LOGDIR}/format_reference.log"
}

check_tool() {
    if ! command -v "$1" &> /dev/null; then
        echo "ERROR: Required tool not found: $1"
        exit 1
    fi
}

# Check required tools
check_tool bcftools
check_tool plink2
check_tool python3

# =============================================================================
# Step 0: Prepare Input VCF
# =============================================================================

log "Step 0: Preparing input data..."

# If input is PLINK, convert to VCF first
if [[ "$INPUT_TYPE" == "plink" ]]; then
    log "Converting PLINK to VCF..."

    plink2 --bfile "$INPUT_FILE" \
        --export vcf-4.2 bgz \
        --out "${OUTPUT_DIR}/temp_input" \
        --threads "$THREADS"

    bcftools index "${OUTPUT_DIR}/temp_input.vcf.gz"
    INPUT_VCF="${OUTPUT_DIR}/temp_input.vcf.gz"
else
    INPUT_VCF="$INPUT_FILE"
fi

# Get list of samples in VCF
bcftools query -l "$INPUT_VCF" > "${OUTPUT_DIR}/vcf_samples.txt"
VCF_SAMPLES=$(wc -l < "${OUTPUT_DIR}/vcf_samples.txt")
log "Found ${VCF_SAMPLES} samples in input VCF"

# Validate sample map against VCF
log "Validating sample map..."
awk 'NR>0 {print $1}' "$SAMPLE_MAP" | sort > "${OUTPUT_DIR}/map_samples.txt"
comm -23 <(sort "${OUTPUT_DIR}/vcf_samples.txt") <(sort "${OUTPUT_DIR}/map_samples.txt") > "${OUTPUT_DIR}/samples_not_in_map.txt"

NOT_IN_MAP=$(wc -l < "${OUTPUT_DIR}/samples_not_in_map.txt")
if [[ $NOT_IN_MAP -gt 0 ]]; then
    log "WARNING: ${NOT_IN_MAP} samples in VCF not found in sample map"
    head "${OUTPUT_DIR}/samples_not_in_map.txt"
fi

# =============================================================================
# Step 1: RFMix v2 Format
# =============================================================================

log ""
log "=== Step 1: Creating RFMix v2 Format ==="

RFMIX2_DIR="${OUTPUT_DIR}/rfmix2"

# Copy/link VCF (RFMix v2 uses VCF directly)
log "Linking reference VCF..."
if [[ "$INPUT_VCF" != "${RFMIX2_DIR}/reference_haplotypes.vcf.gz" ]]; then
    ln -sf "$(readlink -f "$INPUT_VCF")" "${RFMIX2_DIR}/reference_haplotypes.vcf.gz"

    # Copy or create index
    if [[ -f "${INPUT_VCF}.tbi" ]]; then
        ln -sf "$(readlink -f "${INPUT_VCF}.tbi")" "${RFMIX2_DIR}/reference_haplotypes.vcf.gz.tbi"
    elif [[ -f "${INPUT_VCF}.csi" ]]; then
        ln -sf "$(readlink -f "${INPUT_VCF}.csi")" "${RFMIX2_DIR}/reference_haplotypes.vcf.gz.csi"
    else
        bcftools index "${RFMIX2_DIR}/reference_haplotypes.vcf.gz"
    fi
fi

# Create sample map (tab-separated: SampleID<TAB>Population)
log "Creating RFMix v2 sample map..."
awk -F'\t' 'NR>0 && NF>=2 {print $1"\t"$2}' "$SAMPLE_MAP" > "${RFMIX2_DIR}/reference_sample_map.txt"

# Count populations
POPS=$(awk '{print $2}' "${RFMIX2_DIR}/reference_sample_map.txt" | sort -u | wc -l)
log "RFMix v2 format complete: ${POPS} populations"

# =============================================================================
# Step 2: FLARE Format
# =============================================================================

log ""
log "=== Step 2: Creating FLARE Format ==="

FLARE_DIR="${OUTPUT_DIR}/flare"

# Link VCF (FLARE also uses VCF)
ln -sf "$(readlink -f "$INPUT_VCF")" "${FLARE_DIR}/reference_haplotypes.vcf.gz"

if [[ -f "${INPUT_VCF}.tbi" ]]; then
    ln -sf "$(readlink -f "${INPUT_VCF}.tbi")" "${FLARE_DIR}/reference_haplotypes.vcf.gz.tbi"
elif [[ -f "${INPUT_VCF}.csi" ]]; then
    ln -sf "$(readlink -f "${INPUT_VCF}.csi")" "${FLARE_DIR}/reference_haplotypes.vcf.gz.csi"
fi

# Create panel file (space-separated: SampleID<SPACE>Population)
log "Creating FLARE panel file..."
awk -F'\t' 'NR>0 && NF>=2 {print $1" "$2}' "$SAMPLE_MAP" > "${FLARE_DIR}/flare_panels.txt"

log "FLARE format complete"

# =============================================================================
# Step 3: RFMix v1 Format (Per Chromosome)
# =============================================================================

log ""
log "=== Step 3: Creating RFMix v1 Format ==="

RFMIX1_DIR="${OUTPUT_DIR}/rfmix1"

# Check if conversion script exists
CONVERT_SCRIPT="$(dirname "$0")/convert_vcf_to_rfmix1.py"
if [[ ! -f "$CONVERT_SCRIPT" ]]; then
    log "ERROR: Conversion script not found: $CONVERT_SCRIPT"
    log "Skipping RFMix v1 format..."
else
    for chr in $CHROMOSOMES; do
        log "Processing chromosome ${chr}..."

        # Find genetic map for this chromosome
        GENETIC_MAP=""
        for pattern in "${GENETIC_MAP_DIR}/chr${chr}.map" \
                       "${GENETIC_MAP_DIR}/genetic_map_chr${chr}"*.txt \
                       "${GENETIC_MAP_DIR}/chr${chr}"*.gmap; do
            if [[ -f "$pattern" ]]; then
                GENETIC_MAP="$pattern"
                break
            fi
        done

        if [[ -z "$GENETIC_MAP" ]]; then
            log "WARNING: No genetic map found for chr${chr}, skipping..."
            continue
        fi

        # Extract chromosome from VCF
        CHR_VCF="${OUTPUT_DIR}/temp_chr${chr}.vcf.gz"
        bcftools view -r "chr${chr},${chr}" "$INPUT_VCF" -Oz -o "$CHR_VCF" --threads "$THREADS"
        bcftools index "$CHR_VCF"

        # Convert to RFMix v1 format
        python3 "$CONVERT_SCRIPT" \
            --vcf "$CHR_VCF" \
            --sample-map "$SAMPLE_MAP" \
            --genetic-map "$GENETIC_MAP" \
            --output-prefix "${RFMIX1_DIR}/rfmix1_reference" \
            --chromosome "$chr" \
            --verbose 2>&1 | tee -a "${LOGDIR}/rfmix1_chr${chr}.log"

        # Clean up temp file
        rm -f "$CHR_VCF" "${CHR_VCF}.csi"

    done

    log "RFMix v1 format complete"
fi

# =============================================================================
# Step 4: G-NOMIX Format
# =============================================================================

log ""
log "=== Step 4: Creating G-NOMIX Format ==="

GNOMIX_DIR="${OUTPUT_DIR}/gnomix"

# Link VCF
ln -sf "$(readlink -f "$INPUT_VCF")" "${GNOMIX_DIR}/reference_haplotypes.vcf.gz"

if [[ -f "${INPUT_VCF}.tbi" ]]; then
    ln -sf "$(readlink -f "${INPUT_VCF}.tbi")" "${GNOMIX_DIR}/reference_haplotypes.vcf.gz.tbi"
elif [[ -f "${INPUT_VCF}.csi" ]]; then
    ln -sf "$(readlink -f "${INPUT_VCF}.csi")" "${GNOMIX_DIR}/reference_haplotypes.vcf.gz.csi"
fi

# Create sample map (tab-separated, with header for G-NOMIX)
log "Creating G-NOMIX sample map..."
echo -e "Sample\tPopulation" > "${GNOMIX_DIR}/sample_map.tsv"
awk -F'\t' 'NR>0 && NF>=2 {print $1"\t"$2}' "$SAMPLE_MAP" >> "${GNOMIX_DIR}/sample_map.tsv"

log "G-NOMIX format complete"
log "NOTE: To use G-NOMIX, either download pre-trained models or train on this reference"

# =============================================================================
# Step 5: ADMIXTURE Format (LD-Pruned PLINK)
# =============================================================================

log ""
log "=== Step 5: Creating ADMIXTURE Format ==="

ADMIX_DIR="${OUTPUT_DIR}/admixture"

log "Converting VCF to PLINK..."
plink2 --vcf "$INPUT_VCF" \
    --make-bed \
    --out "${ADMIX_DIR}/reference_full" \
    --threads "$THREADS"

log "LD pruning for ADMIXTURE (window=${LD_PRUNE_WINDOW}, step=${LD_PRUNE_STEP}, r2=${LD_PRUNE_R2})..."
plink2 --bfile "${ADMIX_DIR}/reference_full" \
    --indep-pairwise $LD_PRUNE_WINDOW $LD_PRUNE_STEP $LD_PRUNE_R2 \
    --out "${ADMIX_DIR}/prune_list" \
    --threads "$THREADS"

plink2 --bfile "${ADMIX_DIR}/reference_full" \
    --extract "${ADMIX_DIR}/prune_list.prune.in" \
    --make-bed \
    --out "${ADMIX_DIR}/reference_ldpruned" \
    --threads "$THREADS"

# Count variants
FULL_VARS=$(wc -l < "${ADMIX_DIR}/reference_full.bim")
PRUNED_VARS=$(wc -l < "${ADMIX_DIR}/reference_ldpruned.bim")
log "LD pruning: ${FULL_VARS} → ${PRUNED_VARS} variants"

# Create .pop file for supervised ADMIXTURE (if desired)
log "Creating population file for supervised ADMIXTURE..."
# Match sample order in FAM file
awk 'NR==FNR {pop[$1]=$2; next} {print pop[$2]}' \
    "$SAMPLE_MAP" \
    "${ADMIX_DIR}/reference_ldpruned.fam" \
    > "${ADMIX_DIR}/reference_ldpruned.pop"

log "ADMIXTURE format complete"

# =============================================================================
# Step 6: Copy Genetic Maps
# =============================================================================

log ""
log "=== Step 6: Organizing Genetic Maps ==="

mkdir -p "${OUTPUT_DIR}/genetic_maps"
cp -r "${GENETIC_MAP_DIR}"/* "${OUTPUT_DIR}/genetic_maps/" 2>/dev/null || \
    ln -sf "$(readlink -f "$GENETIC_MAP_DIR")"/* "${OUTPUT_DIR}/genetic_maps/"

log "Genetic maps copied"

# =============================================================================
# Cleanup
# =============================================================================

log ""
log "=== Cleanup ==="

if [[ -f "${OUTPUT_DIR}/temp_input.vcf.gz" ]]; then
    rm -f "${OUTPUT_DIR}/temp_input.vcf.gz" "${OUTPUT_DIR}/temp_input.vcf.gz.csi"
fi

rm -f "${OUTPUT_DIR}"/*.txt 2>/dev/null || true

# =============================================================================
# Summary
# =============================================================================

log ""
log "=============================================="
log "Reference Panel Formatting Complete!"
log "=============================================="

cat << EOF | tee -a "${LOGDIR}/format_reference.log"

Output Directory: ${OUTPUT_DIR}

Created Formats:
----------------

rfmix2/
├── reference_haplotypes.vcf.gz     # Phased VCF
├── reference_haplotypes.vcf.gz.tbi # Index
└── reference_sample_map.txt        # Tab-separated: SampleID<TAB>Population

flare/
├── reference_haplotypes.vcf.gz     # Phased VCF
└── flare_panels.txt                # Space-separated: SampleID<SPACE>Population

rfmix1/
├── rfmix1_reference_chr1.alleles   # Binary alleles (one haplotype per line)
├── rfmix1_reference_chr1.classes   # Population labels (space-separated)
├── rfmix1_reference_chr1.snp_locations  # Genetic positions (cM)
└── ... (per chromosome)

gnomix/
├── reference_haplotypes.vcf.gz     # Phased VCF
└── sample_map.tsv                  # Tab-separated with header

admixture/
├── reference_full.bed/bim/fam      # Full reference
├── reference_ldpruned.bed/bim/fam  # LD-pruned for ADMIXTURE
└── reference_ldpruned.pop          # Population labels for supervised mode

genetic_maps/
└── ... (per chromosome)


Copy to Pipeline:
-----------------

cp -r ${OUTPUT_DIR}/rfmix2/* /path/to/pipeline/resources/ancestry_references/
cp -r ${OUTPUT_DIR}/flare/* /path/to/pipeline/resources/ancestry_references/
cp -r ${OUTPUT_DIR}/rfmix1/* /path/to/pipeline/resources/ancestry_references/
cp -r ${OUTPUT_DIR}/gnomix/* /path/to/pipeline/resources/ancestry_references/
cp -r ${OUTPUT_DIR}/admixture/* /path/to/pipeline/resources/ancestry_references/
cp -r ${OUTPUT_DIR}/genetic_maps/* /path/to/pipeline/resources/genetic_maps/

EOF

log ""
log "Done!"
