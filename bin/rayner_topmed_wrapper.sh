#!/bin/bash
# Modified Will Rayner wrapper for TOPMed alignment
# Performs ONLY strand/allele checks, NO frequency filtering
# Version: 2.0 - Updated for pipeline v2

set -euo pipefail

BFILE=$1
RAYNER_SCRIPT=$2
TOPMED_REF=$3
OUTPUT_PREFIX=$4

echo "================================================"
echo "Will Rayner TOPMed Alignment Check"
echo "Mode: Strand/Allele only (NO QC filtering)"
echo "================================================"

# Step 1: Calculate frequencies (required by Rayner script)
plink --bfile ${BFILE} \
    --freq \
    --out ${OUTPUT_PREFIX}_freq

# Step 2: Run Rayner script
perl ${RAYNER_SCRIPT} \
    -b ${BFILE}.bim \
    -f ${OUTPUT_PREFIX}_freq.frq \
    -r ${TOPMED_REF} \
    -h \
    -o ${OUTPUT_PREFIX}

# Step 3: Filter the Run-plink.sh to keep ONLY alignment fixes
# Remove lines that contain QC filtering commands
if [[ -f "Run-plink.sh" ]]; then
    echo "Filtering Run-plink.sh to remove QC steps..."
    
    # Keep only: --flip, --a1-allele, --exclude (for impossible fixes)
    # Remove: --maf, --geno, --hwe, --mind (these are QC, not alignment)
    
    grep -v "\\--maf" Run-plink.sh | \
    grep -v "\\--geno" | \
    grep -v "\\--hwe" | \
    grep -v "\\--mind" > ${OUTPUT_PREFIX}_alignment_only.sh
    
    # Make executable
    chmod +x ${OUTPUT_PREFIX}_alignment_only.sh
    
    echo "Created: ${OUTPUT_PREFIX}_alignment_only.sh"
    echo "This script contains ONLY strand/allele alignment fixes"
    
    # Execute alignment fixes
    bash ${OUTPUT_PREFIX}_alignment_only.sh
    
else
    echo "WARNING: Run-plink.sh not created by Rayner script"
    exit 1
fi

# Step 4: Generate summary report
cat > ${OUTPUT_PREFIX}_rayner_summary.txt <<EOF
================================================
Will Rayner TOPMed Alignment Summary
================================================
Input: ${BFILE}
Reference: TOPMed

ACTIONS TAKEN:
$(grep -c "flip" ${OUTPUT_PREFIX}_alignment_only.sh || echo "0") strand flips
$(grep -c "a1-allele" ${OUTPUT_PREFIX}_alignment_only.sh || echo "0") allele assignments
$(grep -c "exclude" ${OUTPUT_PREFIX}_alignment_only.sh || echo "0") variants excluded (unfixable)

QC FILTERS: NONE (performed in Module 1 separately)

Output: ${OUTPUT_PREFIX}-updated.{bed,bim,fam}
================================================
EOF

cat ${OUTPUT_PREFIX}_rayner_summary.txt

echo "✓ Rayner alignment complete"
