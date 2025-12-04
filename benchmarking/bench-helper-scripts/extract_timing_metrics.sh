#!/bin/bash
################################################################################
# extract_timing_metrics.sh
#
# Extracts timing metrics from benchmark runs for comparison.
# Parses timing.log files from each approach and generates summary.
#
# Usage:
#   ./extract_timing_metrics.sh \
#       --results-dir /path/to/benchmark_results \
#       --output timing_comparison.txt
#
################################################################################

set -euo pipefail

RESULTS_DIR=""
OUTPUT_FILE="timing_comparison.txt"

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Extracts timing metrics from benchmark runs.

Required:
  -r, --results-dir DIR    Directory containing benchmark results

Optional:
  -o, --output FILE        Output file (default: timing_comparison.txt)
  -h, --help               Show this help

Expected directory structure:
  results_dir/
  ├── approach_a_topmed/logs/timing.log
  ├── approach_b_michigan/logs/timing.log
  ├── ours_1step_topmed/logs/timing.log
  └── ...

Output format:
  Approach, Step, Duration (seconds), Duration (human-readable)

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -r|--results-dir)
            RESULTS_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_FILE="$2"
            shift 2
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

if [[ -z "$RESULTS_DIR" ]]; then
    echo "ERROR: --results-dir is required"
    print_usage
    exit 1
fi

# =============================================================================
# Functions
# =============================================================================

seconds_to_human() {
    local seconds=$1
    local days=$((seconds / 86400))
    local hours=$(((seconds % 86400) / 3600))
    local minutes=$(((seconds % 3600) / 60))
    local secs=$((seconds % 60))

    if [[ $days -gt 0 ]]; then
        echo "${days}d ${hours}h ${minutes}m"
    elif [[ $hours -gt 0 ]]; then
        echo "${hours}h ${minutes}m ${secs}s"
    elif [[ $minutes -gt 0 ]]; then
        echo "${minutes}m ${secs}s"
    else
        echo "${secs}s"
    fi
}

parse_timing_log() {
    local log_file=$1
    local approach=$2

    if [[ ! -f "$log_file" ]]; then
        echo "  WARNING: No timing log for ${approach}"
        return
    fi

    # Parse START/END pairs
    declare -A start_times
    declare -A end_times

    while IFS=' ' read -r step action timestamp; do
        if [[ "$action" == "START" ]]; then
            start_times[$step]=$timestamp
        elif [[ "$action" == "END" ]]; then
            end_times[$step]=$timestamp
        fi
    done < "$log_file"

    # Calculate durations
    for step in "${!start_times[@]}"; do
        if [[ -n "${end_times[$step]:-}" ]]; then
            duration=$((${end_times[$step]} - ${start_times[$step]}))
            human=$(seconds_to_human $duration)
            echo "${approach},${step},${duration},${human}"
        fi
    done
}

# =============================================================================
# Main
# =============================================================================

echo "=============================================="
echo "Extracting Timing Metrics"
echo "=============================================="
echo "Results directory: ${RESULTS_DIR}"
echo "Output file: ${OUTPUT_FILE}"
echo ""

# Header
echo "Approach,Step,Duration_Seconds,Duration_Human" > "$OUTPUT_FILE"

# Find all timing logs
total_times=""

for approach_dir in "${RESULTS_DIR}"/*/; do
    approach=$(basename "$approach_dir")

    # Check for timing log in various locations
    timing_log=""
    for possible in "logs/timing.log" "timing.log" ".nextflow.log"; do
        if [[ -f "${approach_dir}${possible}" ]]; then
            timing_log="${approach_dir}${possible}"
            break
        fi
    done

    if [[ -n "$timing_log" ]]; then
        echo "Processing: ${approach}"
        parse_timing_log "$timing_log" "$approach" >> "$OUTPUT_FILE"

        # Calculate total time
        total_start=$(grep "START" "$timing_log" 2>/dev/null | head -1 | awk '{print $3}' || echo "")
        total_end=$(grep "END" "$timing_log" 2>/dev/null | tail -1 | awk '{print $3}' || echo "")

        if [[ -n "$total_start" && -n "$total_end" ]]; then
            total_duration=$((total_end - total_start))
            total_human=$(seconds_to_human $total_duration)
            echo "${approach},TOTAL,${total_duration},${total_human}" >> "$OUTPUT_FILE"
            total_times="${total_times}${approach}: ${total_human}\n"
        fi
    else
        echo "Skipping: ${approach} (no timing log found)"
    fi
done

# =============================================================================
# Summary
# =============================================================================

echo ""
echo "=============================================="
echo "TIMING SUMMARY"
echo "=============================================="

# Generate summary table
echo ""
echo "Total Pipeline Run Times:"
echo "--------------------------"
printf "$total_times" | sort -t: -k2 -n

echo ""
echo "Per-Step Breakdown:"
echo "-------------------"

# Group by step across approaches
for step in STEP1 STEP2 STEP3 STEP4 STEP5 QC IMPUTATION MERGE; do
    step_data=$(grep ",${step}" "$OUTPUT_FILE" 2>/dev/null || true)
    if [[ -n "$step_data" ]]; then
        echo ""
        echo "${step}:"
        echo "$step_data" | while IFS=',' read -r approach step_name duration human; do
            printf "  %-30s %s\n" "$approach" "$human"
        done
    fi
done

echo ""
echo "Full results in: ${OUTPUT_FILE}"
echo ""

# =============================================================================
# Generate comparison plot data
# =============================================================================

PLOT_DATA="${OUTPUT_FILE%.txt}_plot.txt"

echo "# Data for plotting (approach, total_seconds)" > "$PLOT_DATA"
grep ",TOTAL," "$OUTPUT_FILE" | awk -F',' '{print $1, $3}' >> "$PLOT_DATA"

echo "Plot data saved to: ${PLOT_DATA}"
echo ""
echo "To generate plot in R:"
echo "  data <- read.table('${PLOT_DATA}', header=FALSE)"
echo "  barplot(data\$V2, names.arg=data\$V1, main='Pipeline Runtime Comparison')"
echo ""
