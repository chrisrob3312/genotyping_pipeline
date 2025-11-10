#!/bin/bash
#
# validate_samplesheet.sh
#
# Validates sample sheet format and file existence before running pipeline
#
# Usage: bash validate_samplesheet.sh sample_sheet.csv
#

set -euo pipefail

# Check if sample sheet provided
if [[ $# -eq 0 ]]; then
    echo "ERROR: No sample sheet provided!"
    echo "Usage: bash validate_samplesheet.sh sample_sheet.csv"
    exit 1
fi

SAMPLE_SHEET="$1"

# Check if file exists
if [[ ! -f "${SAMPLE_SHEET}" ]]; then
    echo "ERROR: Sample sheet not found: ${SAMPLE_SHEET}"
    exit 1
fi

echo "========================================"
echo "Validating Sample Sheet"
echo "========================================"
echo "File: ${SAMPLE_SHEET}"
echo ""

# Check if CSV has header
first_line=$(head -1 "${SAMPLE_SHEET}")
expected_header="platform_id,batch_id,input_path,file_type,build,file_structure"

if [[ "${first_line}" != "${expected_header}" ]]; then
    echo "ERROR: Invalid header!"
    echo "Expected: ${expected_header}"
    echo "Got:      ${first_line}"
    exit 1
fi

echo "✓ Header valid"
echo ""

# Initialize counters
n_entries=0
n_errors=0
n_warnings=0

# Valid options
valid_file_types=("plink" "vcf")
valid_builds=("hg19" "hg38")
valid_structures=("individual_samples" "individual_chr_split" "merged_batch" "merged_chr_split")

echo "Checking entries..."
echo ""

# Read and validate each line (skip header)
tail -n +2 "${SAMPLE_SHEET}" | while IFS=',' read -r platform batch path type build structure; do
    
    n_entries=$((n_entries + 1))
    
    echo "Entry ${n_entries}: ${platform}/${batch}"
    echo "  Path: ${path}"
    echo "  Type: ${type} | Build: ${build} | Structure: ${structure}"
    
    # Check platform_id not empty
    if [[ -z "${platform}" ]]; then
        echo "  ✗ ERROR: platform_id is empty"
        n_errors=$((n_errors + 1))
    fi
    
    # Check batch_id not empty
    if [[ -z "${batch}" ]]; then
        echo "  ✗ ERROR: batch_id is empty"
        n_errors=$((n_errors + 1))
    fi
    
    # Check input_path not empty
    if [[ -z "${path}" ]]; then
        echo "  ✗ ERROR: input_path is empty"
        n_errors=$((n_errors + 1))
    fi
    
    # Validate file_type
    if [[ ! " ${valid_file_types[@]} " =~ " ${type} " ]]; then
        echo "  ✗ ERROR: Invalid file_type '${type}'. Must be: plink, vcf"
        n_errors=$((n_errors + 1))
    fi
    
    # Validate build
    if [[ ! " ${valid_builds[@]} " =~ " ${build} " ]]; then
        echo "  ✗ ERROR: Invalid build '${build}'. Must be: hg19, hg38"
        n_errors=$((n_errors + 1))
    fi
    
    # Validate file_structure
    if [[ ! " ${valid_structures[@]} " =~ " ${structure} " ]]; then
        echo "  ✗ ERROR: Invalid file_structure '${structure}'"
        echo "          Must be: individual_samples, individual_chr_split, merged_batch, merged_chr_split"
        n_errors=$((n_errors + 1))
    fi
    
    # Check if path exists
    if [[ ! -z "${path}" ]]; then
        if [[ "${structure}" == "individual_samples" ]] || [[ "${structure}" == "individual_chr_split" ]]; then
            # Should be a directory
            if [[ ! -d "${path}" ]]; then
                echo "  ✗ ERROR: Directory not found: ${path}"
                n_errors=$((n_errors + 1))
            else
                echo "  ✓ Directory exists"
                
                # Check for files
                if [[ "${type}" == "plink" ]]; then
                    n_bed=$(find "${path}" -name "*.bed" -type f | wc -l)
                    if [[ ${n_bed} -eq 0 ]]; then
                        echo "  ⚠ WARNING: No .bed files found in directory"
                        n_warnings=$((n_warnings + 1))
                    else
                        echo "  ✓ Found ${n_bed} .bed files"
                    fi
                else
                    n_vcf=$(find "${path}" -name "*.vcf.gz" -type f | wc -l)
                    if [[ ${n_vcf} -eq 0 ]]; then
                        echo "  ⚠ WARNING: No .vcf.gz files found in directory"
                        n_warnings=$((n_warnings + 1))
                    else
                        echo "  ✓ Found ${n_vcf} .vcf.gz files"
                    fi
                fi
            fi
        else
            # Should be a file prefix
            if [[ "${type}" == "plink" ]]; then
                test_file="${path}.bed"
                if [[ ! -f "${test_file}" ]]; then
                    echo "  ✗ ERROR: File not found: ${test_file}"
                    n_errors=$((n_errors + 1))
                else
                    echo "  ✓ Files exist"
                    
                    # Check all three files
                    for ext in bed bim fam; do
                        if [[ ! -f "${path}.${ext}" ]]; then
                            echo "  ✗ ERROR: Missing ${path}.${ext}"
                            n_errors=$((n_errors + 1))
                        fi
                    done
                fi
            else
                test_file="${path}.vcf.gz"
                if [[ ! -f "${test_file}" ]]; then
                    echo "  ✗ ERROR: File not found: ${test_file}"
                    n_errors=$((n_errors + 1))
                else
                    echo "  ✓ File exists"
                    
                    # Check if indexed
                    if [[ ! -f "${test_file}.tbi" ]] && [[ ! -f "${test_file}.csi" ]]; then
                        echo "  ⚠ WARNING: VCF not indexed (.tbi or .csi missing)"
                        n_warnings=$((n_warnings + 1))
                    fi
                fi
            fi
        fi
    fi
    
    echo ""
    
done

echo "========================================"
echo "Validation Summary"
echo "========================================"
echo "Entries checked: ${n_entries}"
echo "Errors: ${n_errors}"
echo "Warnings: ${n_warnings}"
echo ""

if [[ ${n_errors} -gt 0 ]]; then
    echo "✗ VALIDATION FAILED"
    echo "Please fix errors before running pipeline"
    exit 1
elif [[ ${n_warnings} -gt 0 ]]; then
    echo "⚠ VALIDATION PASSED WITH WARNINGS"
    echo "Review warnings - pipeline may still run but might encounter issues"
    exit 0
else
    echo "✓ VALIDATION PASSED"
    echo "Sample sheet is ready for pipeline execution"
    exit 0
fi
