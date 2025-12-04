#!/bin/bash
################################################################################
# CONTAINER VERIFICATION SCRIPT
################################################################################
#
# LOCATION: scripts/verify_containers.sh
# PURPOSE: Comprehensive testing of all Apptainer containers after build
# CALLED BY: Module0_Apptainer_Build.nf (validateContainers process)
#
# WHAT THIS TESTS:
# - Container executability
# - All required software tools are present
# - Software versions are correct
# - Critical plugins/modules are available
# - Dependencies are properly linked
#
# USAGE:
#   bash scripts/verify_containers.sh /path/to/containers/
#
# OUTPUT:
#   - STDOUT: Color-coded test results
#   - validation_report.txt: Summary of all tests
#   - validation_detailed.log: Detailed test output
#   - container_inventory.csv: Spreadsheet of all tools
#
################################################################################

set -uo pipefail  # Don't use -e because we want to continue testing after failures

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Check if container directory provided
if [ "$#" -ne 1 ]; then
    echo -e "${RED}ERROR: Missing container directory${NC}"
    echo "Usage: $0 /path/to/containers/"
    exit 1
fi

CONTAINER_DIR="$1"

if [ ! -d "$CONTAINER_DIR" ]; then
    echo -e "${RED}ERROR: Container directory does not exist: $CONTAINER_DIR${NC}"
    exit 1
fi

# Create validation output directory
VALIDATION_DIR="${CONTAINER_DIR}/validation_logs"
mkdir -p "$VALIDATION_DIR"

REPORT="${VALIDATION_DIR}/validation_report.txt"
DETAILED_LOG="${VALIDATION_DIR}/validation_detailed.log"
INVENTORY="${VALIDATION_DIR}/container_inventory.csv"

# Initialize report files
cat > "$REPORT" << 'EOF'
================================================================================
CONTAINER VALIDATION REPORT
================================================================================
Generated: $(date)

EOF

echo "Timestamp,Container,Test,Result,Details" > "$INVENTORY"

# Counters
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

################################################################################
# HELPER FUNCTIONS
################################################################################

# Test a command in a container
test_container_cmd() {
    local container="$1"
    local test_name="$2"
    local cmd="$3"
    local expected_pattern="$4"  # Optional: regex pattern to match in output
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    echo "" >> "$DETAILED_LOG"
    echo "=== Testing: $test_name ===" >> "$DETAILED_LOG"
    echo "Container: $container" >> "$DETAILED_LOG"
    echo "Command: $cmd" >> "$DETAILED_LOG"
    
    # Run the test
    if output=$(apptainer exec "$container" bash -c "$cmd" 2>&1); then
        echo "Output: $output" >> "$DETAILED_LOG"
        
        # Check for expected pattern if provided
        if [ -n "$expected_pattern" ]; then
            if echo "$output" | grep -qE "$expected_pattern"; then
                echo -e "  ${GREEN}✓${NC} $test_name"
                echo "Result: PASS (matched pattern: $expected_pattern)" >> "$DETAILED_LOG"
                echo "$(date '+%Y-%m-%d %H:%M:%S'),$container,$test_name,PASS,$output" >> "$INVENTORY"
                PASSED_TESTS=$((PASSED_TESTS + 1))
                return 0
            else
                echo -e "  ${RED}✗${NC} $test_name (pattern not matched)"
                echo "Result: FAIL (expected pattern not found: $expected_pattern)" >> "$DETAILED_LOG"
                echo "$(date '+%Y-%m-%d %H:%M:%S'),$container,$test_name,FAIL,Pattern not matched" >> "$INVENTORY"
                FAILED_TESTS=$((FAILED_TESTS + 1))
                return 1
            fi
        else
            echo -e "  ${GREEN}✓${NC} $test_name"
            echo "Result: PASS" >> "$DETAILED_LOG"
            echo "$(date '+%Y-%m-%d %H:%M:%S'),$container,$test_name,PASS,$output" >> "$INVENTORY"
            PASSED_TESTS=$((PASSED_TESTS + 1))
            return 0
        fi
    else
        echo -e "  ${RED}✗${NC} $test_name (command failed)"
        echo "Output: $output" >> "$DETAILED_LOG"
        echo "Result: FAIL" >> "$DETAILED_LOG"
        echo "$(date '+%Y-%m-%d %H:%M:%S'),$container,$test_name,FAIL,Command failed" >> "$INVENTORY"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
}

################################################################################
# BEGIN TESTING
################################################################################

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  CONTAINER VALIDATION${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

################################################################################
# CONTAINER 1: PLINK 1.9
################################################################################
echo -e "${BLUE}[1/7]${NC} Testing plink_1.9.sif..."
PLINK19="${CONTAINER_DIR}/plink_1.9.sif"

if [ ! -f "$PLINK19" ]; then
    echo -e "  ${RED}✗${NC} Container not found: plink_1.9.sif"
    echo "FAIL: plink_1.9.sif not found" >> "$REPORT"
else
    test_container_cmd "$PLINK19" "PLINK 1.9 version" "plink --version" "PLINK v1\.9"
    test_container_cmd "$PLINK19" "PLINK 1.9 basic commands" "plink --help | head -5" "PLINK.*Usage"
fi

################################################################################
# CONTAINER 2: PLINK 2.0
################################################################################
echo -e "\n${BLUE}[2/7]${NC} Testing plink_2.0.sif..."
PLINK20="${CONTAINER_DIR}/plink_2.0.sif"

if [ ! -f "$PLINK20" ]; then
    echo -e "  ${RED}✗${NC} Container not found: plink_2.0.sif"
    echo "FAIL: plink_2.0.sif not found" >> "$REPORT"
else
    test_container_cmd "$PLINK20" "PLINK 2.0 version" "plink2 --version" "PLINK v2\."
    test_container_cmd "$PLINK20" "PLINK 2.0 basic commands" "plink2 --help | head -5" "PLINK.*Usage"
fi

################################################################################
# CONTAINER 3: BCFTOOLS CORE (CRITICAL - has +fixref plugin)
################################################################################
echo -e "\n${BLUE}[3/7]${NC} Testing bcftools.sif..."
BCFTOOLS="${CONTAINER_DIR}/bcftools.sif"

if [ ! -f "$BCFTOOLS" ]; then
    echo -e "  ${RED}✗${NC} Container not found: bcftools.sif"
    echo "FAIL: bcftools.sif not found" >> "$REPORT"
else
    test_container_cmd "$BCFTOOLS" "bcftools version" "bcftools --version" "bcftools 1\."
    test_container_cmd "$BCFTOOLS" "bcftools +fixref plugin" "bcftools plugin -l" "fixref.*Fix REF"
    test_container_cmd "$BCFTOOLS" "vcftools" "vcftools --version" "VCFtools"
    test_container_cmd "$BCFTOOLS" "htslib tabix" "tabix --version" "tabix.*htslib"
    test_container_cmd "$BCFTOOLS" "bgzip" "bgzip --version" "bgzip.*htslib"
fi

################################################################################
# CONTAINER 4: PERL + CROSSMAP
################################################################################
echo -e "\n${BLUE}[4/7]${NC} Testing perl_crossmap.sif..."
PERL_CROSSMAP="${CONTAINER_DIR}/perl_crossmap.sif"

if [ ! -f "$PERL_CROSSMAP" ]; then
    echo -e "  ${RED}✗${NC} Container not found: perl_crossmap.sif"
    echo "FAIL: perl_crossmap.sif not found" >> "$REPORT"
else
    test_container_cmd "$PERL_CROSSMAP" "Perl version" "perl --version" "perl.*v5\."
    test_container_cmd "$PERL_CROSSMAP" "CrossMap" "CrossMap -h" "CrossMap.*Usage"
    test_container_cmd "$PERL_CROSSMAP" "Python for CrossMap" "python3 --version" "Python 3\."
    test_container_cmd "$PERL_CROSSMAP" "PLINK in CrossMap container" "plink --version" "PLINK"
    test_container_cmd "$PERL_CROSSMAP" "bcftools in CrossMap container" "bcftools --version" "bcftools"
fi

################################################################################
# CONTAINER 5: PYTHON API (imputationbot + terralab)
################################################################################
echo -e "\n${BLUE}[5/7]${NC} Testing python_api.sif..."
PYTHON_API="${CONTAINER_DIR}/python_api.sif"

if [ ! -f "$PYTHON_API" ]; then
    echo -e "  ${RED}✗${NC} Container not found: python_api.sif"
    echo "FAIL: python_api.sif not found" >> "$REPORT"
else
    test_container_cmd "$PYTHON_API" "Python version" "python --version" "Python 3\."
    test_container_cmd "$PYTHON_API" "imputationbot" "imputationbot --version" "imputationbot"
    test_container_cmd "$PYTHON_API" "terralab CLI" "terralab --version" "terralab"
    test_container_cmd "$PYTHON_API" "requests library" "python -c 'import requests; print(requests.__version__)'" "[0-9]+\."
    test_container_cmd "$PYTHON_API" "pandas library" "python -c 'import pandas; print(pandas.__version__)'" "[0-9]+\."
fi

################################################################################
# CONTAINER 6: R GENETICS (GENESIS + MagicalRsq + tidyverse)
################################################################################
echo -e "\n${BLUE}[6/7]${NC} Testing r_genetics.sif..."
R_GENETICS="${CONTAINER_DIR}/r_genetics.sif"

if [ ! -f "$R_GENETICS" ]; then
    echo -e "  ${RED}✗${NC} Container not found: r_genetics.sif"
    echo "FAIL: r_genetics.sif not found" >> "$REPORT"
else
    test_container_cmd "$R_GENETICS" "R version" "R --version" "R version [4-9]\."
    test_container_cmd "$R_GENETICS" "GENESIS package" "Rscript -e 'library(GENESIS); packageVersion(\"GENESIS\")'" "GENESIS"
    test_container_cmd "$R_GENETICS" "GWASTools package" "Rscript -e 'library(GWASTools); packageVersion(\"GWASTools\")'" "GWASTools"
    test_container_cmd "$R_GENETICS" "tidyverse" "Rscript -e 'library(tidyverse); packageVersion(\"tidyverse\")'" "tidyverse"
    test_container_cmd "$R_GENETICS" "data.table" "Rscript -e 'library(data.table); packageVersion(\"data.table\")'" "data.table"
    test_container_cmd "$R_GENETICS" "ggplot2" "Rscript -e 'library(ggplot2); packageVersion(\"ggplot2\")'" "ggplot2"
fi

################################################################################
# CONTAINER 7: ANCESTRY SUITE (ADMIXTURE + RFMix + FLARE + G-NOMIX)
################################################################################
echo -e "\n${BLUE}[7/7]${NC} Testing ancestry_suite.sif..."
ANCESTRY="${CONTAINER_DIR}/ancestry_suite.sif"

if [ ! -f "$ANCESTRY" ]; then
    echo -e "  ${RED}✗${NC} Container not found: ancestry_suite.sif"
    echo "FAIL: ancestry_suite.sif not found" >> "$REPORT"
else
    test_container_cmd "$ANCESTRY" "ADMIXTURE" "admixture --version" "ADMIXTURE.*Version"
    test_container_cmd "$ANCESTRY" "RFMix v1" "which rfmix_v1" "rfmix_v1"
    test_container_cmd "$ANCESTRY" "RFMix v2" "rfmix --version" "RFMix.*v2"
    test_container_cmd "$ANCESTRY" "Python for FLARE" "python --version" "Python 3\."
    test_container_cmd "$ANCESTRY" "bcftools for ancestry" "bcftools --version" "bcftools"
fi

################################################################################
# GENERATE FINAL REPORT
################################################################################

echo "" >> "$REPORT"
echo "========================================" >> "$REPORT"
echo "TEST SUMMARY" >> "$REPORT"
echo "========================================" >> "$REPORT"
echo "Total tests run: $TOTAL_TESTS" >> "$REPORT"
echo "Passed: $PASSED_TESTS" >> "$REPORT"
echo "Failed: $FAILED_TESTS" >> "$REPORT"
echo "" >> "$REPORT"

if [ $FAILED_TESTS -eq 0 ]; then
    echo "OVERALL STATUS: ✓ ALL CONTAINERS VALIDATED SUCCESSFULLY!" >> "$REPORT"
    echo ""
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}✓ ALL CONTAINERS VALIDATED SUCCESSFULLY!${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo ""
    echo "Detailed results:"
    echo "  Report:    $REPORT"
    echo "  Log:       $DETAILED_LOG"
    echo "  Inventory: $INVENTORY"
    exit 0
else
    echo "OVERALL STATUS: ✗ VALIDATION FAILED - $FAILED_TESTS test(s) failed" >> "$REPORT"
    echo ""
    echo -e "${RED}========================================${NC}"
    echo -e "${RED}✗ VALIDATION FAILED${NC}"
    echo -e "${RED}$FAILED_TESTS test(s) failed${NC}"
    echo -e "${RED}========================================${NC}"
    echo ""
    echo "Check detailed logs:"
    echo "  Report:    $REPORT"
    echo "  Log:       $DETAILED_LOG"
    exit 1
fi
