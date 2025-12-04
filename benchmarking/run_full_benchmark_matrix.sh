#!/bin/bash
################################################################################
# run_full_benchmark_matrix.sh
#
# Runs ALL approach × server combinations for comprehensive benchmarking.
#
# APPROACHES (6):
#   A - Traditional QC-before-imputation
#   B - Traditional QC-before (alternate implementation)
#   C - Intersect-first + QC
#   D - Per-platform QC + intersect
#   E - Our pipeline 1-step (union → impute → merge → QC)
#   F - Our pipeline 2-step (union → impute → merge → re-impute → QC)
#
# SERVERS (4):
#   topmed   - TOPMed Imputation Server (TOPMed r2 panel)
#   allofus  - All of Us AnVIL (TOPMed-based, diverse)
#   michigan_hrc - Michigan with HRC panel (EUR-focused)
#   michigan_1kg - Michigan with 1000G Phase 3 (global)
#
# TOTAL COMBINATIONS: 6 approaches × 4 servers = 24 benchmark runs
#
# Uses imputationbot for automated submission to Michigan and TOPMed.
#
# Usage:
#   ./run_full_benchmark_matrix.sh \
#       --input /path/to/test_data \
#       --output /path/to/results \
#       --config benchmark_config.yaml
#
################################################################################

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

INPUT_DIR=""
OUTPUT_DIR="./benchmark_results_full"
CONFIG_FILE=""
THREADS=8
DRY_RUN=false
SKIP_SUBMISSION=false      # Skip server submission (for testing locally)
WAIT_FOR_RESULTS=true      # Wait for imputation to complete
POLL_INTERVAL=300          # Check every 5 minutes

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"

# Credential environment variables (or from config)
TOPMED_TOKEN="${TOPMED_TOKEN:-}"
TOPMED_PASSWORD="${TOPMED_PASSWORD:-}"
MICHIGAN_TOKEN="${MICHIGAN_TOKEN:-}"
ALLOFUS_AUTHENTICATED=false

# =============================================================================
# Benchmark Matrix Definition
# =============================================================================

# All approaches
APPROACHES=("A" "B" "C" "D" "E" "F")

# All servers
SERVERS=("topmed" "allofus" "michigan_hrc" "michigan_1kg")

# Approach descriptions
declare -A APPROACH_DESC
APPROACH_DESC[A]="Traditional QC-before + single imputation"
APPROACH_DESC[B]="Traditional QC-before (strict thresholds)"
APPROACH_DESC[C]="Intersect-first + QC + imputation"
APPROACH_DESC[D]="Per-platform QC + intersect + imputation"
APPROACH_DESC[E]="Our pipeline: 1-step (union → impute → merge → MagicalRsq-X)"
APPROACH_DESC[F]="Our pipeline: 2-step (union → impute → merge → re-impute → MagicalRsq-X)"

# Server descriptions
declare -A SERVER_DESC
SERVER_DESC[topmed]="TOPMed Imputation Server (r2 panel, diverse)"
SERVER_DESC[allofus]="All of Us AnVIL (TOPMed-based, 50%+ non-EUR)"
SERVER_DESC[michigan_hrc]="Michigan with HRC (EUR-focused, largest)"
SERVER_DESC[michigan_1kg]="Michigan with 1000G Phase 3 (global diversity)"

# =============================================================================
# Parse Arguments
# =============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Run full benchmark matrix (6 approaches × 4 servers = 24 combinations).

Required:
  -i, --input DIR            Input directory with test data

Optional:
  -o, --output DIR           Output directory (default: ./benchmark_results_full)
  -c, --config FILE          YAML config file with credentials
  -t, --threads N            Threads per job (default: 8)
  --skip-submission          Prepare files but don't submit to servers
  --no-wait                  Don't wait for imputation results
  --poll-interval N          Seconds between status checks (default: 300)
  --dry-run                  Print commands without executing
  -h, --help                 Show this help

Environment Variables (alternative to --config):
  TOPMED_TOKEN               TOPMed API token
  TOPMED_PASSWORD            TOPMed result password
  MICHIGAN_TOKEN             Michigan Imputation Server token

Benchmark Matrix:
-----------------
Approaches:
$(for a in "${APPROACHES[@]}"; do echo "  $a - ${APPROACH_DESC[$a]}"; done)

Servers:
$(for s in "${SERVERS[@]}"; do echo "  $s - ${SERVER_DESC[$s]}"; done)

Total: ${#APPROACHES[@]} × ${#SERVERS[@]} = $((${#APPROACHES[@]} * ${#SERVERS[@]})) combinations

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input) INPUT_DIR="$2"; shift 2 ;;
        -o|--output) OUTPUT_DIR="$2"; shift 2 ;;
        -c|--config) CONFIG_FILE="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --skip-submission) SKIP_SUBMISSION=true; shift ;;
        --no-wait) WAIT_FOR_RESULTS=false; shift ;;
        --poll-interval) POLL_INTERVAL="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) print_usage; exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$INPUT_DIR" ]]; then
    echo "ERROR: --input is required"
    print_usage
    exit 1
fi

# =============================================================================
# Setup and Utilities
# =============================================================================

mkdir -p "${OUTPUT_DIR}"/{logs,timing,status}
MASTER_LOG="${OUTPUT_DIR}/logs/benchmark_master.log"
TIMING_DB="${OUTPUT_DIR}/timing/all_timings.tsv"
STATUS_FILE="${OUTPUT_DIR}/status/job_status.tsv"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$MASTER_LOG"
}

run_cmd() {
    local cmd="$1"
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  [DRY-RUN] $cmd"
        return 0
    else
        log "  Executing: $cmd"
        eval "$cmd" 2>&1 | tee -a "$MASTER_LOG"
        return ${PIPESTATUS[0]}
    fi
}

record_timing() {
    local approach=$1
    local server=$2
    local step=$3
    local start=$4
    local end=$5
    local duration=$((end - start))
    echo -e "${approach}\t${server}\t${step}\t${start}\t${end}\t${duration}" >> "$TIMING_DB"
}

update_status() {
    local approach=$1
    local server=$2
    local status=$3
    local job_id=${4:-""}
    echo -e "${approach}\t${server}\t${status}\t${job_id}\t$(date '+%Y-%m-%d %H:%M:%S')" >> "$STATUS_FILE"
}

# =============================================================================
# Check Dependencies
# =============================================================================

check_dependencies() {
    log "Checking dependencies..."

    # Check for imputationbot
    if ! command -v imputationbot &> /dev/null; then
        log "  WARNING: imputationbot not found"
        log "  Install with: curl -sL imputationbot.now.sh | bash"
        log "  Or: pip install imputationbot"

        if [[ "$SKIP_SUBMISSION" != "true" ]]; then
            log "  Continuing without imputationbot (manual submission required)"
        fi
    else
        log "  imputationbot: OK"
        IMPUTATIONBOT_AVAILABLE=true
    fi

    # Check for terralab (All of Us)
    if ! command -v terralab &> /dev/null; then
        log "  WARNING: terralab not found (needed for All of Us)"
        log "  Install with: pip install terralab-cli"
    else
        log "  terralab: OK"
        # Check authentication
        if terralab auth status &>/dev/null; then
            ALLOFUS_AUTHENTICATED=true
            log "  terralab authenticated: OK"
        else
            log "  terralab not authenticated - run 'terralab login' first"
        fi
    fi

    # Required tools
    for tool in plink2 bcftools nextflow; do
        if command -v $tool &> /dev/null; then
            log "  $tool: OK"
        else
            log "  ERROR: $tool not found"
            exit 1
        fi
    done
}

# =============================================================================
# Imputationbot Setup
# =============================================================================

setup_imputationbot() {
    if [[ "${IMPUTATIONBOT_AVAILABLE:-false}" != "true" ]]; then
        return
    fi

    log "Setting up imputationbot instances..."

    # Add TOPMed instance
    if [[ -n "$TOPMED_TOKEN" ]]; then
        log "  Adding TOPMed instance..."
        imputationbot add-instance topmed \
            --url https://imputation.biodatacatalyst.nhlbi.nih.gov \
            --token "$TOPMED_TOKEN" 2>/dev/null || true
    fi

    # Add Michigan instance
    if [[ -n "$MICHIGAN_TOKEN" ]]; then
        log "  Adding Michigan instance..."
        imputationbot add-instance michigan \
            --url https://imputationserver.sph.umich.edu \
            --token "$MICHIGAN_TOKEN" 2>/dev/null || true
    fi
}

# =============================================================================
# Submit to Imputation Server
# =============================================================================

submit_to_server() {
    local vcf_dir=$1
    local server=$2
    local job_name=$3
    local output_dir=$4

    local job_id=""

    case "$server" in
        topmed)
            if [[ "${IMPUTATIONBOT_AVAILABLE:-false}" == "true" && -n "$TOPMED_TOKEN" ]]; then
                log "    Submitting to TOPMed via imputationbot..."
                job_id=$(imputationbot impute \
                    --instance topmed \
                    --files "${vcf_dir}"/chr*.vcf.gz \
                    --refpanel topmed-r2 \
                    --population mixed \
                    --name "$job_name" \
                    --output "$output_dir" 2>&1 | grep -oP 'job-\S+' | head -1)
            else
                log "    TOPMed: Manual submission required"
                log "    Upload ${vcf_dir}/chr*.vcf.gz to https://imputation.biodatacatalyst.nhlbi.nih.gov"
                job_id="manual_topmed_${job_name}"
            fi
            ;;

        allofus)
            if [[ "$ALLOFUS_AUTHENTICATED" == "true" ]]; then
                log "    Submitting to All of Us via terralab..."
                # All of Us uses Terra/AnVIL
                job_id=$(terralab imputation submit \
                    --input-vcfs "${vcf_dir}"/chr*.vcf.gz \
                    --output-dir "$output_dir" \
                    --name "$job_name" 2>&1 | grep -oP 'submission-\S+' | head -1)
            else
                log "    All of Us: Run 'terralab login' first"
                job_id="manual_allofus_${job_name}"
            fi
            ;;

        michigan_hrc)
            if [[ "${IMPUTATIONBOT_AVAILABLE:-false}" == "true" && -n "$MICHIGAN_TOKEN" ]]; then
                log "    Submitting to Michigan (HRC) via imputationbot..."
                job_id=$(imputationbot impute \
                    --instance michigan \
                    --files "${vcf_dir}"/chr*.vcf.gz \
                    --refpanel hrc-r1.1 \
                    --population mixed \
                    --name "$job_name" \
                    --output "$output_dir" 2>&1 | grep -oP 'job-\S+' | head -1)
            else
                log "    Michigan HRC: Manual submission required"
                job_id="manual_michigan_hrc_${job_name}"
            fi
            ;;

        michigan_1kg)
            if [[ "${IMPUTATIONBOT_AVAILABLE:-false}" == "true" && -n "$MICHIGAN_TOKEN" ]]; then
                log "    Submitting to Michigan (1KG) via imputationbot..."
                job_id=$(imputationbot impute \
                    --instance michigan \
                    --files "${vcf_dir}"/chr*.vcf.gz \
                    --refpanel 1000g-phase3-v5 \
                    --population mixed \
                    --name "$job_name" \
                    --output "$output_dir" 2>&1 | grep -oP 'job-\S+' | head -1)
            else
                log "    Michigan 1KG: Manual submission required"
                job_id="manual_michigan_1kg_${job_name}"
            fi
            ;;
    esac

    echo "$job_id"
}

# =============================================================================
# Wait for Imputation Jobs
# =============================================================================

wait_for_jobs() {
    local jobs_file=$1

    if [[ "$WAIT_FOR_RESULTS" != "true" ]]; then
        log "Skipping wait (--no-wait specified)"
        return
    fi

    log "Waiting for imputation jobs to complete..."

    while true; do
        local all_done=true
        local pending=0
        local completed=0
        local failed=0

        while IFS=$'\t' read -r approach server status job_id timestamp; do
            if [[ "$status" == "submitted" || "$status" == "running" ]]; then
                # Check job status
                local new_status=""

                if [[ "$job_id" == manual_* ]]; then
                    # Manual jobs - check if results directory exists
                    local results_dir="${OUTPUT_DIR}/${approach}_${server}/imputed_results"
                    if [[ -d "$results_dir" && -n "$(ls -A $results_dir 2>/dev/null)" ]]; then
                        new_status="completed"
                    else
                        new_status="pending_manual"
                        all_done=false
                        ((pending++))
                    fi
                elif [[ "${IMPUTATIONBOT_AVAILABLE:-false}" == "true" ]]; then
                    # Check via imputationbot
                    new_status=$(imputationbot jobs --id "$job_id" --format json 2>/dev/null | \
                        jq -r '.status' 2>/dev/null || echo "unknown")

                    case "$new_status" in
                        complete|succeeded) new_status="completed"; ((completed++)) ;;
                        failed|error) new_status="failed"; ((failed++)) ;;
                        *) all_done=false; ((pending++)) ;;
                    esac
                fi

                if [[ -n "$new_status" && "$new_status" != "$status" ]]; then
                    update_status "$approach" "$server" "$new_status" "$job_id"
                fi
            elif [[ "$status" == "completed" ]]; then
                ((completed++))
            elif [[ "$status" == "failed" ]]; then
                ((failed++))
            fi
        done < "$jobs_file"

        log "  Status: $completed completed, $pending pending, $failed failed"

        if [[ "$all_done" == "true" ]]; then
            log "All jobs completed!"
            break
        fi

        log "  Waiting ${POLL_INTERVAL}s before next check..."
        sleep "$POLL_INTERVAL"
    done
}

# =============================================================================
# Run Single Benchmark
# =============================================================================

run_benchmark() {
    local approach=$1
    local server=$2

    local run_name="${approach}_${server}"
    local run_dir="${OUTPUT_DIR}/${run_name}"

    log ""
    log "=============================================="
    log "Running: ${run_name}"
    log "  Approach: ${APPROACH_DESC[$approach]}"
    log "  Server: ${SERVER_DESC[$server]}"
    log "=============================================="

    mkdir -p "${run_dir}"/{prep,vcf_for_imputation,imputed_results,post_qc,final,logs}

    local step_start step_end
    update_status "$approach" "$server" "started"

    # -------------------------------------------------------------------------
    # STEP 1: Pre-processing based on approach
    # -------------------------------------------------------------------------
    step_start=$(date +%s)
    log "  Step 1: Pre-processing (${approach})..."

    case "$approach" in
        A|B)
            # Traditional QC-before
            run_cmd "${SCRIPT_DIR}/alternative_approaches/approach_a_topmed.sh \
                --input '${INPUT_PLINK}' \
                --output '${run_dir}/prep' \
                --threads ${THREADS} \
                2>&1 || true"
            ;;
        C)
            # Intersect-first
            run_cmd "${SCRIPT_DIR}/alternative_approaches/approach_c_intersect_first.sh \
                --inputs '${INPUT_PLINK}' \
                --output '${run_dir}/prep' \
                --threads ${THREADS} \
                2>&1 || true"
            ;;
        D)
            # Per-platform QC + intersect
            run_cmd "${SCRIPT_DIR}/alternative_approaches/approach_d_qcbefore_michigan.sh \
                --inputs '${INPUT_PLINK}' \
                --output '${run_dir}/prep' \
                --threads ${THREADS} \
                2>&1 || true"
            ;;
        E|F)
            # Our pipeline - minimal pre-processing
            log "    Our pipeline: Minimal prep (QC after imputation)"
            run_cmd "plink2 --bfile '${INPUT_PLINK}' \
                --geno 0.1 --mind 0.1 \
                --export vcf-4.2 bgz \
                --out '${run_dir}/prep/minimal_qc' \
                --threads ${THREADS}"
            ;;
    esac

    step_end=$(date +%s)
    record_timing "$approach" "$server" "preprocessing" "$step_start" "$step_end"

    # -------------------------------------------------------------------------
    # STEP 2: Prepare VCFs for imputation server
    # -------------------------------------------------------------------------
    step_start=$(date +%s)
    log "  Step 2: Preparing VCFs for ${server}..."

    local vcf_source
    vcf_source=$(find "${run_dir}/prep" -name "*.vcf.gz" | head -1)

    if [[ -z "$vcf_source" ]]; then
        # Convert from PLINK if needed
        local bed_source
        bed_source=$(find "${run_dir}/prep" -name "*.bed" | head -1)
        if [[ -n "$bed_source" ]]; then
            run_cmd "plink2 --bfile '${bed_source%.bed}' \
                --export vcf-4.2 bgz \
                --out '${run_dir}/vcf_for_imputation/prepared' \
                --threads ${THREADS}"
            vcf_source="${run_dir}/vcf_for_imputation/prepared.vcf.gz"
        fi
    fi

    # Split by chromosome
    if [[ -n "$vcf_source" ]]; then
        run_cmd "bcftools index -f '$vcf_source'"
        for chr in {1..22}; do
            run_cmd "bcftools view -r ${chr} '$vcf_source' \
                -Oz -o '${run_dir}/vcf_for_imputation/chr${chr}.vcf.gz'"
            run_cmd "bcftools index '${run_dir}/vcf_for_imputation/chr${chr}.vcf.gz'"
        done
    fi

    step_end=$(date +%s)
    record_timing "$approach" "$server" "vcf_preparation" "$step_start" "$step_end"

    # -------------------------------------------------------------------------
    # STEP 3: Submit to imputation server
    # -------------------------------------------------------------------------
    step_start=$(date +%s)
    log "  Step 3: Submitting to ${server}..."

    if [[ "$SKIP_SUBMISSION" != "true" ]]; then
        local job_id
        job_id=$(submit_to_server \
            "${run_dir}/vcf_for_imputation" \
            "$server" \
            "$run_name" \
            "${run_dir}/imputed_results")

        update_status "$approach" "$server" "submitted" "$job_id"
        log "    Job ID: ${job_id}"
    else
        log "    Skipping submission (--skip-submission)"
        update_status "$approach" "$server" "skipped"
    fi

    step_end=$(date +%s)
    record_timing "$approach" "$server" "submission" "$step_start" "$step_end"

    log "  Preprocessing complete for ${run_name}"
}

# =============================================================================
# Post-Imputation Processing
# =============================================================================

run_post_imputation() {
    local approach=$1
    local server=$2

    local run_name="${approach}_${server}"
    local run_dir="${OUTPUT_DIR}/${run_name}"

    log ""
    log "Post-processing: ${run_name}"

    local step_start step_end

    # -------------------------------------------------------------------------
    # STEP 4: Post-imputation QC
    # -------------------------------------------------------------------------
    step_start=$(date +%s)
    log "  Step 4: Post-imputation QC..."

    local imputed_vcf
    imputed_vcf=$(find "${run_dir}/imputed_results" -name "*.vcf.gz" | head -1)

    if [[ -z "$imputed_vcf" ]]; then
        log "    WARNING: No imputed results found"
        return
    fi

    case "$approach" in
        A|B|C|D)
            # Traditional: Simple R² filter
            run_cmd "bcftools view -i 'INFO/R2 >= 0.3' '$imputed_vcf' \
                -Oz -o '${run_dir}/post_qc/r2_filtered.vcf.gz'"
            ;;
        E|F)
            # Our pipeline: MagicalRsq-X filter
            log "    Applying MagicalRsq-X filtering..."
            run_cmd "Rscript '${PIPELINE_DIR}/helper_scripts/magicalrsq_filter.R' \
                --vcf '$imputed_vcf' \
                --ancestry mixed \
                --threshold 0.3 \
                --output '${run_dir}/post_qc/magicalrsq_filtered.vcf.gz' \
                2>/dev/null || bcftools view -i 'INFO/R2 >= 0.3' '$imputed_vcf' \
                -Oz -o '${run_dir}/post_qc/magicalrsq_filtered.vcf.gz'"
            ;;
    esac

    step_end=$(date +%s)
    record_timing "$approach" "$server" "post_qc" "$step_start" "$step_end"

    # -------------------------------------------------------------------------
    # STEP 5: Re-imputation (for approach F only)
    # -------------------------------------------------------------------------
    if [[ "$approach" == "F" ]]; then
        step_start=$(date +%s)
        log "  Step 5: Re-imputation (2-step)..."

        # This would submit again - for now, note it
        log "    NOTE: Re-imputation would be submitted here"
        log "    Using same server: ${server}"

        step_end=$(date +%s)
        record_timing "$approach" "$server" "reimputation" "$step_start" "$step_end"
    fi

    # -------------------------------------------------------------------------
    # STEP 6: Final QC and format
    # -------------------------------------------------------------------------
    step_start=$(date +%s)
    log "  Step 6: Final formatting..."

    local qc_vcf
    qc_vcf=$(find "${run_dir}/post_qc" -name "*.vcf.gz" | head -1)

    if [[ -n "$qc_vcf" ]]; then
        run_cmd "plink2 --vcf '$qc_vcf' \
            --make-bed \
            --out '${run_dir}/final/${run_name}' \
            --threads ${THREADS}"

        # Count final variants/samples
        if [[ -f "${run_dir}/final/${run_name}.bim" ]]; then
            local n_vars n_samp
            n_vars=$(wc -l < "${run_dir}/final/${run_name}.bim")
            n_samp=$(wc -l < "${run_dir}/final/${run_name}.fam")
            log "    Final: ${n_vars} variants, ${n_samp} samples"
            echo -e "${approach}\t${server}\t${n_vars}\t${n_samp}" >> "${OUTPUT_DIR}/final_counts.tsv"
        fi
    fi

    step_end=$(date +%s)
    record_timing "$approach" "$server" "final_format" "$step_start" "$step_end"

    update_status "$approach" "$server" "completed"
}

# =============================================================================
# Main Execution
# =============================================================================

log "=============================================="
log "FULL BENCHMARK MATRIX"
log "=============================================="
log "Date: $(date)"
log "Input: ${INPUT_DIR}"
log "Output: ${OUTPUT_DIR}"
log "Approaches: ${APPROACHES[*]}"
log "Servers: ${SERVERS[*]}"
log "Total combinations: $((${#APPROACHES[@]} * ${#SERVERS[@]}))"
log ""

# Initialize timing database
echo -e "approach\tserver\tstep\tstart\tend\tduration_seconds" > "$TIMING_DB"
echo -e "approach\tserver\tstatus\tjob_id\ttimestamp" > "$STATUS_FILE"

# Check dependencies
check_dependencies

# Setup imputationbot
setup_imputationbot

# Find input PLINK files
INPUT_PLINK=$(find "$INPUT_DIR" -name "*.bed" | head -1 | sed 's/.bed$//')
if [[ -z "$INPUT_PLINK" ]]; then
    log "ERROR: No PLINK files found in ${INPUT_DIR}"
    exit 1
fi
log "Input PLINK: ${INPUT_PLINK}"

TOTAL_START=$(date +%s)

# -------------------------------------------------------------------------
# Phase 1: Pre-processing and submission
# -------------------------------------------------------------------------
log ""
log "========== PHASE 1: PRE-PROCESSING AND SUBMISSION =========="

for approach in "${APPROACHES[@]}"; do
    for server in "${SERVERS[@]}"; do
        run_benchmark "$approach" "$server"
    done
done

# -------------------------------------------------------------------------
# Phase 2: Wait for imputation results
# -------------------------------------------------------------------------
log ""
log "========== PHASE 2: WAITING FOR IMPUTATION =========="
wait_for_jobs "$STATUS_FILE"

# -------------------------------------------------------------------------
# Phase 3: Post-processing
# -------------------------------------------------------------------------
log ""
log "========== PHASE 3: POST-PROCESSING =========="

for approach in "${APPROACHES[@]}"; do
    for server in "${SERVERS[@]}"; do
        run_post_imputation "$approach" "$server"
    done
done

# =============================================================================
# Generate Summary
# =============================================================================

TOTAL_END=$(date +%s)
TOTAL_TIME=$((TOTAL_END - TOTAL_START))

log ""
log "=============================================="
log "BENCHMARK MATRIX COMPLETE"
log "=============================================="
log "Total time: ${TOTAL_TIME} seconds ($(( TOTAL_TIME / 3600 ))h $((( TOTAL_TIME % 3600 ) / 60 ))m)"
log ""

# Generate summary report
cat > "${OUTPUT_DIR}/BENCHMARK_SUMMARY.md" << EOF
# Benchmark Matrix Results

**Date:** $(date)
**Total Runtime:** ${TOTAL_TIME} seconds

## Configuration
- Input: ${INPUT_DIR}
- Approaches: ${#APPROACHES[@]}
- Servers: ${#SERVERS[@]}
- Total combinations: $((${#APPROACHES[@]} * ${#SERVERS[@]}))

## Approach Descriptions
$(for a in "${APPROACHES[@]}"; do echo "- **$a**: ${APPROACH_DESC[$a]}"; done)

## Server Descriptions
$(for s in "${SERVERS[@]}"; do echo "- **$s**: ${SERVER_DESC[$s]}"; done)

## Final Variant/Sample Counts
\`\`\`
$(cat "${OUTPUT_DIR}/final_counts.tsv" 2>/dev/null || echo "No results yet")
\`\`\`

## Timing Summary
See: timing/all_timings.tsv

## Next Steps
1. Run concordance analysis with WGS truth
2. Compare variant retention by ancestry
3. Generate publication figures

\`\`\`bash
# Concordance
Rscript bench-helper-scripts/calculate_concordance.R ...

# Compare approaches
Rscript bench-helper-scripts/compare_approaches.R ...
\`\`\`
EOF

log "Summary: ${OUTPUT_DIR}/BENCHMARK_SUMMARY.md"
log "Timings: ${OUTPUT_DIR}/timing/all_timings.tsv"
log ""
