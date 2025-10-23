#!/bin/bash

################################################################################
# Apptainer Container Build Script
# Builds all containers needed for the genotyping imputation pipeline
################################################################################

set -e  # Exit on error

# Configuration
BUILD_DIR="./apptainer_containers"
OUTPUT_DIR="./container_images"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}======================================================================"
echo -e "Building Apptainer Containers for Genotyping Pipeline"
echo -e "======================================================================${NC}"
echo ""

# Create directories
mkdir -p "${BUILD_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Function to build a container
build_container() {
    local def_file=$1
    local container_name=$2
    
    echo -e "${YELLOW}Building ${container_name}...${NC}"
    
    if [ ! -f "${def_file}" ]; then
        echo -e "${RED}ERROR: Definition file ${def_file} not found!${NC}"
        return 1
    fi
    
    apptainer build --force \
        "${OUTPUT_DIR}/${container_name}.sif" \
        "${def_file}"
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ ${container_name} built successfully!${NC}"
        echo -e "   Location: ${OUTPUT_DIR}/${container_name}.sif"
        echo -e "   Size: $(du -h ${OUTPUT_DIR}/${container_name}.sif | cut -f1)"
        echo ""
    else
        echo -e "${RED}✗ Failed to build ${container_name}${NC}"
        return 1
    fi
}

# Build all containers
echo -e "${BLUE}Starting container builds...${NC}"
echo ""

# 1. PLINK tools
build_container "plink.def" "plink"

# 2. R genetics
build_container "r_genetics.def" "r_genetics"

# 3. Perl and VCF tools
build_container "perl_vcftools.def" "perl_vcftools"

# 4. Ancestry tools
build_container "ancestry_tools.def" "ancestry_tools"

# 5. Python tools
build_container "python_tools.def" "python_tools"

echo -e "${GREEN}======================================================================"
echo -e "All containers built successfully!"
echo -e "======================================================================${NC}"
echo ""

# Summary
echo -e "${BLUE}Container Summary:${NC}"
echo ""
ls -lh "${OUTPUT_DIR}"/*.sif | awk '{print "  " $9 " - " $5}'
echo ""

# Create manifest file
echo -e "${BLUE}Creating container manifest...${NC}"
cat > "${OUTPUT_DIR}/container_manifest.txt" << EOF
Apptainer Container Manifest
Generated: $(date)
==================================================

Container Details:
--------------------------------------------------

1. plink.sif
   Tools: PLINK 1.9, PLINK 2.0
   Purpose: Genotyping QC and analysis
   Used in: Modules 1, 3, 4, 6

2. r_genetics.sif
   Tools: R 4.3.1, GENESIS, MagicalRsq, ggplot2, data.table
   Purpose: Statistical genetics analysis
   Used in: Modules 3, 6, 7

3. perl_vcftools.sif
   Tools: Perl, bcftools, vcftools, CrossMap, samtools
   Purpose: VCF manipulation and file conversion
   Used in: Modules 1, 3, 4, 5, 6

4. ancestry_tools.sif
   Tools: ADMIXTURE, RFMix v1/v2, FLARE, g-nomix, Graf-anc, SHAPEIT4, Eagle
   Purpose: Global and local ancestry inference
   Used in: Module 7

5. python_tools.sif
   Tools: Python 3.10, pandas, pysam, API libraries
   Purpose: API interactions and utilities
   Used in: Modules 2, 5 (TOPMed/AnVIL API)

==================================================
Build Information:
  Build Date: $(date)
  Build Host: $(hostname)
  Apptainer Version: $(apptainer --version)

Container Sizes:
$(ls -lh "${OUTPUT_DIR}"/*.sif | awk '{print "  " $9 " - " $5}')

==================================================
Usage Instructions:

To use these containers in Nextflow:
  1. Copy all .sif files to your HPC shared storage or cloud bucket
  2. Update nextflow.config with the correct paths
  3. Run pipeline with: nextflow run main.nf -profile apptainer

To test a container:
  apptainer exec <container_name>.sif <command>

Examples:
  apptainer exec plink.sif plink --version
  apptainer exec r_genetics.sif Rscript my_script.R
  apptainer exec ancestry_tools.sif admixture --version

==================================================
EOF

echo -e "${GREEN}Manifest created: ${OUTPUT_DIR}/container_manifest.txt${NC}"
echo ""

# Optional: Test containers
echo -e "${YELLOW}Do you want to test the containers? (y/n)${NC}"
read -r response

if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
    echo -e "${BLUE}Testing containers...${NC}"
    echo ""
    
    echo "Testing plink.sif..."
    apptainer exec "${OUTPUT_DIR}/plink.sif" plink --version
    
    echo ""
    echo "Testing r_genetics.sif..."
    apptainer exec "${OUTPUT_DIR}/r_genetics.sif" R --version
    
    echo ""
    echo "Testing perl_vcftools.sif..."
    apptainer exec "${OUTPUT_DIR}/perl_vcftools.sif" bcftools --version
    
    echo ""
    echo "Testing ancestry_tools.sif..."
    apptainer exec "${OUTPUT_DIR}/ancestry_tools.sif" admixture --version
    
    echo ""
    echo "Testing python_tools.sif..."
    apptainer exec "${OUTPUT_DIR}/python_tools.sif" python3 --version
    
    echo -e "${GREEN}All container tests passed!${NC}"
fi

echo ""
echo -e "${GREEN}======================================================================"
echo -e "Container Build Complete!"
echo -e "======================================================================${NC}"
echo ""
echo -e "Next steps:"
echo -e "  1. Review ${OUTPUT_DIR}/container_manifest.txt"
echo -e "  2. Copy .sif files to your deployment location:"
echo -e "     - HPC: cp ${OUTPUT_DIR}/*.sif /shared/containers/"
echo -e "     - Cloud: gsutil cp ${OUTPUT_DIR}/*.sif gs://your-bucket/containers/"
echo -e "  3. Update nextflow.config with container paths"
echo -e "  4. Run your pipeline!"
echo ""
