# Test Data and Benchmark Datasets

This directory contains test data for validating the genotyping pipeline and instructions for obtaining benchmark datasets.

## Table of Contents

1. [Quick Test Data](#quick-test-data) - Small files for CI/CD testing
2. [Module 1 Input Testing](#module-1-input-testing) - Test all input formats
3. [Benchmark Datasets](#benchmark-datasets) - WGS truth + array data for validation
4. [LAI Reference Data](#lai-reference-data) - Reference panels for ancestry inference

---

## Quick Test Data

### 1000 Genomes Subset (Recommended for Testing)

Download a small subset of 1000 Genomes data for quick testing:

```bash
# Create test data directory
cd test_data

# Download chromosome 22 (smallest, good for testing)
# 1000 Genomes Phase 3 - ~2500 samples
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi

# Download sample information (population labels)
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
```

### Create Smaller Test Subsets

```bash
# Extract 100 samples from diverse populations for quick tests
# AFR: YRI, LWK; EUR: CEU, GBR; EAS: CHB, JPT; SAS: GIH; AMR: MXL

# Get sample IDs for each population
awk '$2=="YRI" || $2=="CEU" || $2=="CHB" || $2=="GIH" {print $1}' \
    integrated_call_samples_v3.20130502.ALL.panel | head -100 > test_samples.txt

# Subset VCF to these samples and first 10,000 variants
bcftools view -S test_samples.txt \
    ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | \
    bcftools view -r 22:16000000-18000000 -Oz -o test_chr22_subset.vcf.gz

bcftools index test_chr22_subset.vcf.gz

# Convert to PLINK format
plink2 --vcf test_chr22_subset.vcf.gz \
    --make-bed \
    --out test_chr22_subset
```

---

## Module 1 Input Testing

Module 1 (Pre-Imputation QC) accepts multiple input formats. Create test files for each:

### Input Format Matrix

| Format | File Structure | Build | Test File |
|--------|---------------|-------|-----------|
| PLINK binary | Single set | hg19 | `test_hg19_single/` |
| PLINK binary | Single set | hg38 | `test_hg38_single/` |
| PLINK binary | Per-chromosome | hg19 | `test_hg19_perchr/` |
| PLINK binary | Per-chromosome | hg38 | `test_hg38_perchr/` |
| VCF | Single file | hg19 | `test_hg19.vcf.gz` |
| VCF | Single file | hg38 | `test_hg38.vcf.gz` |
| VCF | Per-chromosome | hg19 | `test_hg19_perchr_vcf/` |

### Create All Test Input Formats

```bash
#!/bin/bash
# create_test_inputs.sh - Generate test files for all Module 1 input formats

set -euo pipefail

INPUT_VCF="test_chr22_subset.vcf.gz"
PREFIX="pipeline_test"

echo "=== Creating test input formats for Module 1 ==="

# ------------------------------------
# 1. PLINK Binary - Single Set (hg19)
# ------------------------------------
echo "Creating PLINK single-set hg19..."
mkdir -p test_hg19_single

plink2 --vcf ${INPUT_VCF} \
    --make-bed \
    --out test_hg19_single/${PREFIX}

# ------------------------------------
# 2. PLINK Binary - Per-Chromosome (hg19)
# ------------------------------------
echo "Creating PLINK per-chromosome hg19..."
mkdir -p test_hg19_perchr

# Since we only have chr22, just copy it
cp test_hg19_single/${PREFIX}.bed test_hg19_perchr/${PREFIX}_chr22.bed
cp test_hg19_single/${PREFIX}.bim test_hg19_perchr/${PREFIX}_chr22.bim
cp test_hg19_single/${PREFIX}.fam test_hg19_perchr/${PREFIX}_chr22.fam

# ------------------------------------
# 3. PLINK Binary - Single Set (hg38 via liftover simulation)
# ------------------------------------
echo "Creating PLINK single-set hg38 (simulated)..."
mkdir -p test_hg38_single

# For testing, just copy and rename (in production, use CrossMap)
cp test_hg19_single/${PREFIX}.bed test_hg38_single/${PREFIX}.bed
cp test_hg19_single/${PREFIX}.fam test_hg38_single/${PREFIX}.fam

# Modify chromosome positions slightly to simulate hg38
# (This is for testing only - real data would use CrossMap)
awk '{$4=$4+10000; print}' test_hg19_single/${PREFIX}.bim > test_hg38_single/${PREFIX}.bim

# ------------------------------------
# 4. VCF - Single File (hg19)
# ------------------------------------
echo "Creating VCF single-file hg19..."
cp ${INPUT_VCF} test_hg19.vcf.gz
cp ${INPUT_VCF}.csi test_hg19.vcf.gz.csi 2>/dev/null || \
    bcftools index test_hg19.vcf.gz

# ------------------------------------
# 5. VCF - Per-Chromosome (hg19)
# ------------------------------------
echo "Creating VCF per-chromosome hg19..."
mkdir -p test_hg19_perchr_vcf

cp ${INPUT_VCF} test_hg19_perchr_vcf/${PREFIX}_chr22.vcf.gz
bcftools index test_hg19_perchr_vcf/${PREFIX}_chr22.vcf.gz

echo ""
echo "=== Test input files created ==="
echo ""
echo "Directory structure:"
find test_* -type f | head -30
```

### Sample Sheet for Testing

Create `test_sample_sheet.csv`:

```csv
platform_id,batch_id,input_path,file_type,build,file_structure
Illumina_Omni,batch1,test_data/test_hg19_single/pipeline_test,plink,hg19,single
Illumina_Omni,batch2,test_data/test_hg19_perchr/pipeline_test,plink,hg19,per_chromosome
Illumina_Omni,batch3,test_data/test_hg38_single/pipeline_test,plink,hg38,single
Illumina_Omni,batch4,test_data/test_hg19.vcf.gz,vcf,hg19,single
Illumina_Omni,batch5,test_data/test_hg19_perchr_vcf/pipeline_test,vcf,hg19,per_chromosome
```

---

## Benchmark Datasets

### For Imputation Accuracy Benchmarking

These datasets have **paired WGS truth data and genotyping array data** for validation:

#### 1. 1000 Genomes + Array Data (Recommended - Public)

**Source:** [1000 Genomes FTP](https://www.internationalgenome.org/data)

| Data Type | Samples | Ancestries | Access |
|-----------|---------|------------|--------|
| WGS (30x) | 2,504 | AFR, EUR, EAS, SAS, AMR | Public FTP |
| Omni 2.5M array | 2,141 | Same | Public FTP |
| Affymetrix 6.0 | 1,248 | Same | Public FTP |

```bash
# Download WGS data (truth set)
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

# Download Omni 2.5M array data
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/broad_intensities/Omni25_genotypes_2141_samples.b37.v2.vcf.gz

# Download Affymetrix 6.0 array data
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/coriell_affy6_intensities/ALL.wgs.nhgri_coriell_affy_6.20131213.snps_indels.chip.genotypes.vcf.gz

# Sample information
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
```

**Benchmark Strategy:**
1. Start with array data (Omni 2.5M or Affy 6.0)
2. Run pipeline imputation
3. Compare imputed genotypes to WGS truth
4. Calculate concordance by MAF bin and ancestry

---

#### 2. GIAB (Genome in a Bottle) - High-Quality Truth

**Source:** [NIST GIAB](https://www.nist.gov/programs-projects/genome-bottle)

| Sample | Ancestry | WGS | Array Data |
|--------|----------|-----|------------|
| HG001/NA12878 | European (CEU) | 300x Illumina | Coriell (various arrays) |
| HG002/NA24385 | Ashkenazi Jewish | 300x Illumina | Coriell |
| HG005/NA24631 | Han Chinese | 300x Illumina | Coriell |

```bash
# Download GIAB truth VCF (high-confidence calls)
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed

# For array data, order from Coriell:
# https://www.coriell.org/0/Sections/Search/Sample_Detail.aspx?Ref=NA12878
```

---

#### 3. BioMe / Mount Sinai (dbGaP - Requires Access)

**dbGaP Accession:** [phs001644](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001644)

| Data Type | Samples | Ancestries |
|-----------|---------|------------|
| WGS | ~30,000 | AFR, EUR, Hispanic/Latino |
| Arrays | Various | Illumina Core, OmniExpress, MEGA, Omni 2.5M |

**Note:** Requires dbGaP application approval.

---

#### 4. HCHS/SOL - Hispanic/Latino Benchmark (dbGaP)

**dbGaP Accession:** phs000810

Best for benchmarking admixed Hispanic/Latino populations with:
- Array genotypes
- TOPMed WGS for truth comparison

---

### Ancestry-Stratified Benchmark Samples

For comprehensive benchmarking across ancestries, use these 1000 Genomes populations:

| Ancestry | Populations | N Samples | Notes |
|----------|-------------|-----------|-------|
| **African (AFR)** | YRI, LWK, GWD, MSL, ESN, ASW, ACB | ~660 | Includes African Americans |
| **European (EUR)** | CEU, TSI, FIN, GBR, IBS | ~500 | |
| **East Asian (EAS)** | CHB, JPT, CHS, CDX, KHV | ~500 | |
| **South Asian (SAS)** | GIH, PJL, BEB, STU, ITU | ~490 | |
| **American (AMR)** | MXL, PUR, CLM, PEL | ~350 | Admixed Latino |

---

## LAI Reference Data

### 1000 Genomes Reference Panel for LAI

```bash
# Download phased reference haplotypes (all chromosomes)
for chr in {1..22}; do
    wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
    wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
done

# Download population labels
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
```

### Genetic Maps (HapMap/1000G)

```bash
# Download genetic maps from SHAPEIT4/Eagle
wget -c https://github.com/odelaneau/shapeit4/raw/master/maps/genetic_maps.b37.tar.gz
tar -xzf genetic_maps.b37.tar.gz -C resources/genetic_maps/

# Or from Eagle
wget -c https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg19_withX.txt.gz
```

### Create LAI Reference Files

```bash
# After downloading 1000 Genomes data, create LAI reference files:

# 1. Create sample map (for RFMix v2, FLARE)
awk 'NR>1 {print $1"\t"$2}' integrated_call_samples_v3.20130502.ALL.panel \
    > resources/ancestry_references/reference_sample_map.txt

# 2. Subset to reference populations only (exclude admixed AMR)
grep -E "YRI|LWK|GWD|MSL|ESN|CEU|TSI|FIN|GBR|IBS|CHB|JPT|CHS|CDX|KHV|GIH|PJL|BEB|STU|ITU" \
    resources/ancestry_references/reference_sample_map.txt \
    > resources/ancestry_references/reference_sample_map_noAMR.txt

# 3. Convert to RFMix v1 format (per chromosome)
for chr in {1..22}; do
    python helper_scripts/convert_vcf_to_rfmix1.py \
        --vcf ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
        --sample-map resources/ancestry_references/reference_sample_map.txt \
        --genetic-map resources/genetic_maps/chr${chr}.map \
        --output-prefix resources/ancestry_references/rfmix1_reference \
        --chromosome ${chr} \
        --verbose
done
```

---

## Running Benchmark Tests

### Full Pipeline Benchmark Script

```bash
#!/bin/bash
# benchmark_pipeline.sh - Run full benchmark comparison

# Configuration
ARRAY_DATA="test_data/omni25_subset.vcf.gz"
WGS_TRUTH="test_data/wgs_truth.vcf.gz"
OUTPUT_DIR="benchmark_results"

# 1. Run pipeline on array data
nextflow run main.nf \
    --sample_sheet test_data/benchmark_sample_sheet.csv \
    --outdir ${OUTPUT_DIR}/pipeline_output \
    --imputation_service topmed \
    -profile slurm

# 2. Compare imputed to WGS truth
bcftools isec \
    -p ${OUTPUT_DIR}/comparison \
    ${OUTPUT_DIR}/pipeline_output/imputed.vcf.gz \
    ${WGS_TRUTH}

# 3. Calculate concordance by MAF bin
Rscript helper_scripts/calculate_concordance.R \
    --imputed ${OUTPUT_DIR}/pipeline_output/imputed.vcf.gz \
    --truth ${WGS_TRUTH} \
    --ancestry resources/ancestry_references/sample_ancestry.txt \
    --output ${OUTPUT_DIR}/benchmark_report.html
```

---

## Summary of Data Sources

| Purpose | Dataset | Access | Diversity |
|---------|---------|--------|-----------|
| **Quick Testing** | 1KG subset (chr22) | Public FTP | 26 populations |
| **Imputation Benchmark** | 1KG + Omni 2.5M | Public FTP | Global |
| **High-Accuracy Truth** | GIAB HG001-HG005 | Public FTP | EUR, ASJ, CHN |
| **Large-Scale Diverse** | BioMe/MLOF | dbGaP | AFR, EUR, Latino |
| **Hispanic/Latino Focus** | HCHS/SOL | dbGaP | Latino subgroups |
| **LAI Reference** | 1KG Phase 3 | Public FTP | 5 superpopulations |

---

## References

- [1000 Genomes Project](https://www.internationalgenome.org/)
- [GIAB - Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle)
- [TOPMed Imputation Server](https://imputation.biodatacatalyst.nhlbi.nih.gov/)
- [Comprehensive SNP Array Evaluation](https://www.nature.com/articles/s41598-022-22215-y)
- [TOPMed Imputation for Diverse Populations](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008500)
