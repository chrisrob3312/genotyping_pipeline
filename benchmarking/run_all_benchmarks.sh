#!/bin/bash
################################################################################
# run_all_benchmarks.sh
#
# Master script to run all benchmark approaches for comparison.
# Runs traditional approaches (A-D) and our pipeline variants (E-F).
#
# Usage:
#   ./run_all_benchmarks.sh \
#       --input /path/to/test_data \
#       --output /path/to/results \
#       --topmed-token YOUR_TOKEN \
#       --approaches "a,b,e,f"
#
# Requirements:
#   - Nextflow installed
#   - API tokens for imputation servers
#   - Test data downloaded (see download_benchmark_data.sh)
#
################################################################################

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

INPUT_DIR=""
OUTPUT_DIR="./benchmark_results"
TOPMED_TOKEN=""
TOPMED_PASSWORD=""
MICHIGAN_TOKEN=""
APPROACHES="all"  # Comma-separated: a,b,c,d,e,f or "all" or "traditional" or "ours"
THREADS=4
DRY_RUN=false

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Run benchmark comparisons across imputation approaches.

Required:
  -i, --input DIR            Input directory with test data
  --topmed-token TOKEN       TOPMed API token (for approaches A, C, E, F)

Recommended:
  --topmed-password PASS     TOPMed password for auto-decrypt
  --michigan-token TOKEN     Michigan API token (for approaches B, D)

Optional:
  -o, --output DIR           Output directory (default: ./benchmark_results)
  --approaches LIST          Approaches to run (default: all)
                             Options: a,b,c,d,e,f or "all", "traditional", "ours"
  -t, --threads N            Threads per job (default: 4)
  --dry-run                  Print commands without executing
  -h, --help                 Show this help

Approaches:
  TRADITIONAL (QC before imputation):
    A - Traditional QC + TOPMed
    B - Traditional QC + Michigan (HRC)
    C - Intersect-first + QC + TOPMed
    D - Per-platform QC + Intersect + Michigan

  OUR PIPELINE (QC after imputation):
    E - 1-step: Union → Impute → Merge → MagicalRsq-X → QC
    F - 2-step: Union → Impute → Merge → Re-impute → MagicalRsq-X → QC

Server Variants (for E and F):
    topmed  - TOPMed imputation server
    allofus - All of Us AnVIL server

Example:
  # Run all approaches with TOPMed
  $0 -i test_data --topmed-token TOKEN --approaches all

  # Run only our pipeline variants
  $0 -i test_data --topmed-token TOKEN --approaches ours

  # Run specific approaches
  $0 -i test_data --topmed-token TOKEN --michigan-token TOKEN2 --approaches a,b,e

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --topmed-token)
            TOPMED_TOKEN="$2"
            shift 2
            ;;
        --topmed-password)
            TOPMED_PASSWORD="$2"
            shift 2
            ;;
        --michigan-token)
            MICHIGAN_TOKEN="$2"
            shift 2
            ;;
        --approaches)
            APPROACHES="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            print_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

if [[ -z "$INPUT_DIR" ]]; then
    echo "ERROR: --input is required"
    print_usage
    exit 1
fi

# =============================================================================
# Setup
# =============================================================================

mkdir -p "${OUTPUT_DIR}/logs"
LOG_FILE="${OUTPUT_DIR}/logs/benchmark_master.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

run_cmd() {
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  [DRY-RUN] $1"
    else
        log "  Running: $1"
        eval "$1" 2>&1 | tee -a "$LOG_FILE"
    fi
}

# Expand approach aliases
case "$APPROACHES" in
    all)
        APPROACH_LIST="a b c d e_topmed e_allofus f_topmed f_allofus"
        ;;
    traditional)
        APPROACH_LIST="a b c d"
        ;;
    ours)
        APPROACH_LIST="e_topmed e_allofus f_topmed f_allofus"
        ;;
    *)
        APPROACH_LIST=$(echo "$APPROACHES" | tr ',' ' ')
        ;;
esac

log "=============================================="
log "BENCHMARK COMPARISON SUITE"
log "=============================================="
log "Input: ${INPUT_DIR}"
log "Output: ${OUTPUT_DIR}"
log "Approaches: ${APPROACH_LIST}"
log "Dry run: ${DRY_RUN}"
log ""

# Find input files
PLINK_PREFIX=$(find "$INPUT_DIR" -name "*.bed" | head -1 | sed 's/.bed$//')
if [[ -z "$PLINK_PREFIX" ]]; then
    log "ERROR: No PLINK files found in ${INPUT_DIR}"
    exit 1
fi
log "Using input: ${PLINK_PREFIX}"

# =============================================================================
# Run Benchmarks
# =============================================================================

BENCHMARK_START=$(date +%s)

for approach in $APPROACH_LIST; do
    log ""
    log "=============================================="
    log "Running Approach: ${approach}"
    log "=============================================="

    APPROACH_START=$(date +%s)

    case "$approach" in
        a|A)
            if [[ -z "$TOPMED_TOKEN" ]]; then
                log "  SKIPPING: No TOPMed token provided"
                continue
            fi
            run_cmd "${SCRIPT_DIR}/alternative_approaches/approach_a_topmed.sh \
                --input '${PLINK_PREFIX}' \
                --output '${OUTPUT_DIR}/approach_a_topmed' \
                --topmed-token '${TOPMED_TOKEN}' \
                --topmed-password '${TOPMED_PASSWORD}' \
                --threads ${THREADS}"
            ;;

        b|B)
            if [[ -z "$MICHIGAN_TOKEN" ]]; then
                log "  SKIPPING: No Michigan token provided"
                continue
            fi
            run_cmd "${SCRIPT_DIR}/alternative_approaches/approach_b_michigan.sh \
                --input '${PLINK_PREFIX}' \
                --output '${OUTPUT_DIR}/approach_b_michigan' \
                --michigan-token '${MICHIGAN_TOKEN}' \
                --threads ${THREADS}"
            ;;

        c|C)
            if [[ -z "$TOPMED_TOKEN" ]]; then
                log "  SKIPPING: No TOPMed token provided"
                continue
            fi
            run_cmd "${SCRIPT_DIR}/alternative_approaches/approach_c_intersect_first.sh \
                --inputs '${PLINK_PREFIX}' \
                --output '${OUTPUT_DIR}/approach_c_intersect' \
                --topmed-token '${TOPMED_TOKEN}' \
                --topmed-password '${TOPMED_PASSWORD}' \
                --threads ${THREADS}"
            ;;

        d|D)
            if [[ -z "$MICHIGAN_TOKEN" ]]; then
                log "  SKIPPING: No Michigan token provided"
                continue
            fi
            run_cmd "${SCRIPT_DIR}/alternative_approaches/approach_d_qcbefore_michigan.sh \
                --inputs '${PLINK_PREFIX}' \
                --output '${OUTPUT_DIR}/approach_d_qcbefore_michigan' \
                --michigan-token '${MICHIGAN_TOKEN}' \
                --threads ${THREADS}"
            ;;

        e|E|e_topmed)
            if [[ -z "$TOPMED_TOKEN" ]]; then
                log "  SKIPPING: No TOPMed token provided"
                continue
            fi
            run_cmd "nextflow run '${PIPELINE_DIR}/main.nf' \
                -c '${SCRIPT_DIR}/our_pipeline_variants/ours_1step_topmed.config' \
                --sample_sheet '${INPUT_DIR}/benchmark_sample_sheet.csv' \
                --topmed_api_token '${TOPMED_TOKEN}' \
                --topmed_password '${TOPMED_PASSWORD}' \
                --outdir '${OUTPUT_DIR}/ours_1step_topmed' \
                -profile benchmark"
            ;;

        e_allofus)
            run_cmd "nextflow run '${PIPELINE_DIR}/main.nf' \
                -c '${SCRIPT_DIR}/our_pipeline_variants/ours_1step_allofus.config' \
                --sample_sheet '${INPUT_DIR}/benchmark_sample_sheet.csv' \
                --outdir '${OUTPUT_DIR}/ours_1step_allofus' \
                -profile benchmark"
            ;;

        f|F|f_topmed)
            if [[ -z "$TOPMED_TOKEN" ]]; then
                log "  SKIPPING: No TOPMed token provided"
                continue
            fi
            run_cmd "nextflow run '${PIPELINE_DIR}/main.nf' \
                -c '${SCRIPT_DIR}/our_pipeline_variants/ours_2step_topmed.config' \
                --sample_sheet '${INPUT_DIR}/benchmark_sample_sheet.csv' \
                --topmed_api_token '${TOPMED_TOKEN}' \
                --topmed_password '${TOPMED_PASSWORD}' \
                --outdir '${OUTPUT_DIR}/ours_2step_topmed' \
                -profile benchmark"
            ;;

        f_allofus)
            run_cmd "nextflow run '${PIPELINE_DIR}/main.nf' \
                -c '${SCRIPT_DIR}/our_pipeline_variants/ours_2step_allofus.config' \
                --sample_sheet '${INPUT_DIR}/benchmark_sample_sheet.csv' \
                --outdir '${OUTPUT_DIR}/ours_2step_allofus' \
                -profile benchmark"
            ;;

        *)
            log "  WARNING: Unknown approach '${approach}', skipping"
            ;;
    esac

    APPROACH_END=$(date +%s)
    APPROACH_TIME=$((APPROACH_END - APPROACH_START))
    log "  Completed in ${APPROACH_TIME} seconds"
done

# =============================================================================
# Post-Processing
# =============================================================================

BENCHMARK_END=$(date +%s)
TOTAL_TIME=$((BENCHMARK_END - BENCHMARK_START))

log ""
log "=============================================="
log "BENCHMARK SUITE COMPLETE"
log "=============================================="
log "Total time: ${TOTAL_TIME} seconds"
log ""

# Extract timing metrics
if [[ "$DRY_RUN" != "true" ]]; then
    log "Extracting timing metrics..."
    run_cmd "${SCRIPT_DIR}/bench-helper-scripts/extract_timing_metrics.sh \
        --results-dir '${OUTPUT_DIR}' \
        --output '${OUTPUT_DIR}/timing_comparison.txt'"
fi

log ""
log "Next steps:"
log "  1. Download imputed results from servers (if not auto-downloaded)"
log "  2. Run concordance analysis:"
log "     Rscript bench-helper-scripts/calculate_concordance.R ..."
log "  3. Compare approaches:"
log "     Rscript bench-helper-scripts/compare_approaches.R ..."
log "  4. Generate publication figures:"
log "     Rscript bench-helper-scripts/generate_publication_figures.R ..."
log ""
log "Results in: ${OUTPUT_DIR}/"
log "Log file: ${LOG_FILE}"
log ""
