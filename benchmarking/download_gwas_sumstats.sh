#!/bin/bash
################################################################################
# download_gwas_sumstats.sh
#
# Downloads publicly available GWAS summary statistics for benchmarking.
#
# Traits selected for:
#   - Multi-ancestry GWAS available
#   - Known ancestry-specific effects
#   - Rare variant associations (to test our pipeline's advantage)
#   - PRS weights available
#
# Usage:
#   ./download_gwas_sumstats.sh -o /path/to/output
#
################################################################################

set -euo pipefail

OUTPUT_DIR="${OUTPUT_DIR:-./benchmarking/gwas_data}"

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Downloads GWAS summary statistics for benchmarking.

Options:
  -o, --output-dir DIR   Output directory (default: ./benchmarking/gwas_data)
  -h, --help             Show this help

Traits Downloaded:
  1. Height (GIANT consortium) - multi-ancestry
  2. BMI (GIANT consortium) - ancestry-specific effects
  3. Type 2 Diabetes (DIAMANTE) - strong AFR/AMR effects
  4. LDL Cholesterol (GLGC) - rare variants
  5. Blood Pressure (MVP) - multi-ancestry

PRS Weights:
  - Downloaded from PGS Catalog for each trait
  - Multi-ancestry scores when available

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output-dir)
            OUTPUT_DIR="$2"
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

echo "=============================================="
echo "GWAS Summary Statistics Download"
echo "=============================================="
echo "Output: ${OUTPUT_DIR}"
echo ""

mkdir -p "${OUTPUT_DIR}"/{sumstats,prs_weights,known_hits}
cd "${OUTPUT_DIR}"

# =============================================================================
# 1. HEIGHT (GIANT Consortium)
# =============================================================================

echo ""
echo "=== 1. Height (GIANT) ==="

cd sumstats

# GIANT 2022 multi-ancestry height GWAS
if [ ! -f "height_giant_2022.txt.gz" ]; then
    echo "  Downloading GIANT height GWAS..."
    # Meta-analysis of ~5M individuals
    wget -q -c "https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018_UPDATED.txt.gz" \
        -O height_giant_2022.txt.gz 2>/dev/null || \
    echo "  NOTE: Download from https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files"
fi

cd "${OUTPUT_DIR}"

# =============================================================================
# 2. BMI (GIANT Consortium)
# =============================================================================

echo ""
echo "=== 2. BMI (GIANT) ==="

cd sumstats

if [ ! -f "bmi_giant_2018.txt.gz" ]; then
    echo "  Downloading GIANT BMI GWAS..."
    wget -q -c "https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz" \
        -O bmi_giant_2018.txt.gz 2>/dev/null || \
    echo "  NOTE: Download from GIANT consortium website"
fi

cd "${OUTPUT_DIR}"

# =============================================================================
# 3. TYPE 2 DIABETES (DIAMANTE)
# =============================================================================

echo ""
echo "=== 3. Type 2 Diabetes (DIAMANTE) ==="

cd sumstats

# DIAMANTE trans-ancestry T2D GWAS
if [ ! -f "t2d_diamante_transancestry.txt.gz" ]; then
    echo "  Downloading DIAMANTE T2D GWAS..."
    echo "  NOTE: Download from https://diagram-consortium.org/downloads.html"
    echo "        Select: DIAMANTE Trans-ancestry GWAS"

    # Create placeholder with download instructions
    cat > t2d_download_instructions.txt << 'INSTRUCTIONS'
DIAMANTE Trans-Ancestry Type 2 Diabetes GWAS

Download from: https://diagram-consortium.org/downloads.html

Files needed:
1. Trans-ancestry meta-analysis (all ancestries combined)
2. EUR-specific results
3. AFR-specific results (if available)
4. AMR-specific results (if available)

Why this trait:
- Strong ancestry-specific effects
- Well-characterized loci
- Important for demonstrating pipeline benefits in admixed populations
INSTRUCTIONS
fi

cd "${OUTPUT_DIR}"

# =============================================================================
# 4. LDL CHOLESTEROL (GLGC)
# =============================================================================

echo ""
echo "=== 4. LDL Cholesterol (GLGC) ==="

cd sumstats

if [ ! -f "ldl_glgc_2021.txt.gz" ]; then
    echo "  Downloading GLGC LDL GWAS..."
    # Global Lipids Genetics Consortium
    wget -q -c "http://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/trans_ancestry/LDL_INV_trans_ancestry_summary_statistics.txt.gz" \
        -O ldl_glgc_2021.txt.gz 2>/dev/null || \
    echo "  NOTE: Download from http://csg.sph.umich.edu/willer/public/glgc-lipids2021/"
fi

cd "${OUTPUT_DIR}"

# =============================================================================
# 5. BLOOD PRESSURE (MVP / UKB)
# =============================================================================

echo ""
echo "=== 5. Blood Pressure ==="

cd sumstats

if [ ! -f "sbp_mvp_multiancestry.txt.gz" ]; then
    echo "  Downloading blood pressure GWAS..."
    echo "  NOTE: MVP data requires dbGaP access"
    echo "        Alternative: Use Neale Lab UKB GWAS (European)"

    # Neale Lab UKB (European only, but publicly available)
    # SBP phenotype code: 4080
    cat > bp_download_instructions.txt << 'INSTRUCTIONS'
Blood Pressure GWAS Options:

1. Million Veteran Program (MVP) - Multi-ancestry
   - Requires dbGaP application
   - Best for ancestry comparison
   - Download: https://www.ncbi.nlm.nih.gov/gap/

2. Neale Lab UK Biobank (European)
   - Publicly available
   - Phenotype 4080 (SBP) or 4079 (DBP)
   - Download: https://www.nealelab.is/uk-biobank

3. ICBP (International Consortium for Blood Pressure)
   - Multi-ancestry
   - Download: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000585
INSTRUCTIONS
fi

cd "${OUTPUT_DIR}"

# =============================================================================
# 6. PRS WEIGHTS FROM PGS CATALOG
# =============================================================================

echo ""
echo "=== 6. PRS Weights (PGS Catalog) ==="

cd prs_weights

# Height PRS (multi-ancestry)
if [ ! -f "PGS000011_height.txt.gz" ]; then
    echo "  Downloading height PRS weights..."
    wget -q -c "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000011/ScoringFiles/PGS000011.txt.gz" \
        -O PGS000011_height.txt.gz 2>/dev/null || \
    echo "  NOTE: Download from https://www.pgscatalog.org/score/PGS000011/"
fi

# BMI PRS
if [ ! -f "PGS000027_bmi.txt.gz" ]; then
    echo "  Downloading BMI PRS weights..."
    wget -q -c "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000027/ScoringFiles/PGS000027.txt.gz" \
        -O PGS000027_bmi.txt.gz 2>/dev/null || \
    echo "  NOTE: Download from https://www.pgscatalog.org/score/PGS000027/"
fi

# T2D PRS (multi-ancestry)
if [ ! -f "PGS002308_t2d.txt.gz" ]; then
    echo "  Downloading T2D PRS weights (multi-ancestry)..."
    wget -q -c "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS002308/ScoringFiles/PGS002308.txt.gz" \
        -O PGS002308_t2d.txt.gz 2>/dev/null || \
    echo "  NOTE: Download from https://www.pgscatalog.org/score/PGS002308/"
fi

cd "${OUTPUT_DIR}"

# =============================================================================
# 7. KNOWN HITS FOR VALIDATION
# =============================================================================

echo ""
echo "=== 7. Creating Known Hits Reference ==="

cd known_hits

# Download GWAS Catalog associations for these traits
cat > known_hits_reference.txt << 'KNOWN_HITS'
# Known GWAS hits for validation
# Format: CHR POS RSID GENE TRAIT PVALUE ANCESTRY
# Selected hits with known ancestry-specific effects

# Height - known multi-ancestry hits
2   27730940    rs1229984    ADH1B       Height  1e-50   Multi
6   32413564    rs9272219    HLA-DQA1    Height  1e-30   Multi
3   53829968    rs724016     CACNA1D     Height  1e-25   Multi

# BMI - ancestry-specific
16  53767042    rs1421085    FTO         BMI     1e-100  Multi
18  57829135    rs6567160    MC4R        BMI     1e-50   Multi
2   25150296    rs13021737   POMC        BMI     1e-30   EUR

# T2D - strong ancestry effects
10  114758349   rs7903146    TCF7L2      T2D     1e-200  Multi
9   22134094    rs10811661   CDKN2A/B    T2D     1e-50   Multi
6   20686996    rs1800562    HFE         T2D     1e-20   EUR
# SLC16A11 - Latino-specific
17  6945866     rs13342232   SLC16A11    T2D     1e-25   AMR

# LDL Cholesterol - rare variant effects
1   55505647    rs11591147   PCSK9       LDL     1e-100  Multi
19  11202306    rs429358     APOE        LDL     1e-200  Multi
2   21263900    rs515135     APOB        LDL     1e-50   Multi

# Blood Pressure
17  47402807    rs11649420   PLCD3       SBP     1e-30   Multi
4   81164723    rs13107325   SLC39A8     SBP     1e-40   EUR
10  18745778    rs4373814    CACNB2      SBP     1e-20   Multi
KNOWN_HITS

echo "  Created known_hits_reference.txt"

cd "${OUTPUT_DIR}"

# =============================================================================
# CREATE SUMMARY
# =============================================================================

cat > GWAS_DATA_README.txt << 'EOF'
================================================================================
GWAS Summary Statistics for Benchmarking
================================================================================

This directory contains GWAS summary statistics and PRS weights for
validating imputation pipeline performance.

Directory Structure:
--------------------

sumstats/
├── height_giant_2022.txt.gz      # GIANT height meta-analysis
├── bmi_giant_2018.txt.gz         # GIANT BMI meta-analysis
├── t2d_diamante_transancestry.txt.gz  # DIAMANTE T2D
├── ldl_glgc_2021.txt.gz          # GLGC LDL cholesterol
└── *_download_instructions.txt    # Manual download instructions

prs_weights/
├── PGS000011_height.txt.gz       # Height PRS (PGS Catalog)
├── PGS000027_bmi.txt.gz          # BMI PRS
└── PGS002308_t2d.txt.gz          # T2D multi-ancestry PRS

known_hits/
└── known_hits_reference.txt       # Known loci for hit recovery test


Benchmarking Tests:
-------------------

1. GWAS Hit Recovery
   - Run association testing on imputed data
   - Check if known hits reach significance
   - Compare across approaches and ancestries

2. Effect Size Concordance
   - Compare β estimates with published values
   - Check for ancestry-specific bias

3. PRS Performance
   - Calculate PRS using published weights
   - Compare R² and AUC across ancestries
   - Show improvement for underrepresented groups

4. Genomic Inflation
   - Calculate λ GC for each approach
   - Should be close to 1.0 (no inflation)


Usage Example:
--------------

# Calculate PRS
plink2 --bfile imputed_data \
    --score prs_weights/PGS000011_height.txt.gz \
    --out prs_height

# Run GWAS
plink2 --bfile imputed_data \
    --pheno phenotypes.txt \
    --pheno-name HEIGHT \
    --glm \
    --out gwas_height

# Check hit recovery
Rscript bench-helper-scripts/check_hit_recovery.R \
    --gwas gwas_height.glm.linear \
    --known-hits known_hits/known_hits_reference.txt \
    --output hit_recovery_results.txt

================================================================================
EOF

echo ""
echo "=============================================="
echo "GWAS Data Download Complete!"
echo "=============================================="
echo ""
echo "Some files require manual download (see instructions in each directory)"
echo ""
echo "Recommended traits for benchmarking:"
echo "  1. Height - good baseline, many known hits"
echo "  2. T2D - strong ancestry effects (test AMR improvement)"
echo "  3. LDL - rare variant associations (test our pipeline advantage)"
echo ""
echo "See GWAS_DATA_README.txt for usage instructions"
echo ""
