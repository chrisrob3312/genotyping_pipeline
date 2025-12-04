# Genotyping Pipeline Testing and Benchmarking Guide

## Table of Contents
1. [Issues to Fix Before Testing](#1-issues-to-fix-before-testing)
2. [Module Testing Checklist](#2-module-testing-checklist)
3. [Benchmarking Workflow](#3-benchmarking-workflow)
4. [Comparative Analysis Against Existing Pipelines](#4-comparative-analysis-against-existing-pipelines)
5. [Publication Strategy](#5-publication-strategy)

---

## 1. Issues to Fix Before Testing

### Critical Issues (Must Fix)

- [ ] **Parameter mismatch: `monitor_interval` vs `monitor_interval_minutes`**
  - Location: `modules/Module5_ReImputation.nf` vs `config/nextflow.config`
  - Fix: Add `params.monitor_interval` to config (in seconds)

- [ ] **Missing `magicalrsqx_dir` parameter**
  - Location: `modules/Module6_PostMergeQC.nf` references it
  - Fix: Add to `config/nextflow.config`:
    ```groovy
    magicalrsqx_dir = 'resources/magicalrsq_models'
    ```

- [ ] **Incorrect `file_structure` values in test data documentation**
  - Location: `test_data/README.md` lines 152-157
  - Fix: Change `"single"` to `"individual_samples"` in CSV examples

- [ ] **Missing `primary_ancestry` default value**
  - Location: `config/nextflow.config`
  - Fix: Add with default:
    ```groovy
    primary_ancestry = 'AFR'  // Default for imputation server selection
    ```

- [ ] **Empty containers directory**
  - Location: `resources/containers/`
  - Fix: Build containers using Module 0 or document manual build process

- [ ] **Missing reference files**
  - FASTA reference genome
  - Chain files (hg19 <-> hg38)
  - TOPMed/HRC reference for Rayner check

### High Priority Issues

- [ ] **Parameter naming inconsistencies**
  - `hwe_pvalue` vs `hwe_threshold` - standardize naming
  - `magicalrsq_models_dir` vs `magicalrsqx_dir` - use consistent name
  - `admixture_k` vs `admixture_k_min/max` - clarify usage

- [ ] **Container name inconsistencies**
  - Config: `plink.sif` vs Documentation: `plink_1.9.sif`
  - Standardize across all files

- [ ] **Missing parameter definitions** (add to config with defaults):
  - `n_pcs` - number of PCs for population stratification
  - `maf_threshold` - MAF filtering threshold
  - `admixture_cv_folds` - cross-validation folds for ADMIXTURE

- [ ] **Placeholder values to update**
  - `your_account` in SLURM profiles
  - `your-project` in examples
  - `your-repo` in documentation

- [ ] **Add VCF parsing library to `convert_vcf_to_rfmix1.py`**
  - Add: `import pysam` or `from cyvcf2 import VCF`

### Medium Priority Issues

- [ ] Create container definition (.def) files in `resources/containers/`
- [ ] Add README files to:
  - `resources/genetic_maps/`
  - `resources/magicalrsq_models/`
- [ ] Remove or implement stub files:
  - `examples/run_*.sh` (empty stubs)
  - `helper_scripts/prepare_references.sh` (empty)

---

## 2. Module Testing Checklist

### Pre-Testing Setup

- [ ] Install Nextflow (v23.04+)
- [ ] Install Apptainer/Singularity
- [ ] Clone repository
- [ ] Download/build all containers
- [ ] Download reference files:
  - [ ] Reference genome FASTA (hg19 and hg38)
  - [ ] Chain files
  - [ ] BRAVO Freeze 10 data
  - [ ] Genetic maps (from TractorWorkflow)
  - [ ] GRAF-anc data files
- [ ] Create test dataset (use `helper_scripts/create_test_data.sh`)

### Test Data Matrix

Test all 8 input format combinations:

| # | Format | Build | Structure | Status |
|---|--------|-------|-----------|--------|
| 1 | PLINK  | hg19  | individual_samples | [ ] |
| 2 | PLINK  | hg19  | combined | [ ] |
| 3 | PLINK  | hg38  | individual_samples | [ ] |
| 4 | PLINK  | hg38  | combined | [ ] |
| 5 | VCF    | hg19  | individual_samples | [ ] |
| 6 | VCF    | hg19  | combined | [ ] |
| 7 | VCF    | hg38  | individual_samples | [ ] |
| 8 | VCF    | hg38  | combined | [ ] |

---

### Module 0: Container Build

**Purpose**: Build Apptainer/Singularity containers

- [ ] All container definitions exist
- [ ] Containers build successfully
- [ ] Containers contain expected tools:
  - [ ] PLINK 1.9 and 2.0
  - [ ] BCFtools
  - [ ] SAMtools
  - [ ] Eagle2
  - [ ] ADMIXTURE
  - [ ] GRAF-anc
  - [ ] Python with pandas, numpy
  - [ ] R with required packages

**Test Command**:
```bash
nextflow run main.nf -profile test,singularity -entry MODULE0_BUILD_CONTAINERS
```

---

### Module 1: Input Validation & Standardization

**Purpose**: Validate inputs, convert formats, coordinate liftover

**Tests**:
- [ ] **Input validation**
  - [ ] Detects invalid sample sheet
  - [ ] Detects missing input files
  - [ ] Validates genome build specification
  - [ ] Validates file format specification

- [ ] **Format conversion**
  - [ ] VCF → PLINK conversion works
  - [ ] PLINK → VCF conversion works
  - [ ] Preserves all variants
  - [ ] Preserves sample IDs

- [ ] **Coordinate liftover**
  - [ ] hg19 → hg38 liftover works
  - [ ] hg38 → hg19 liftover works (if needed)
  - [ ] Unmapped variants logged
  - [ ] Strand flips handled correctly

- [ ] **Chromosome splitting**
  - [ ] Creates per-chromosome files
  - [ ] Handles chrX/chrY correctly
  - [ ] Outputs proper naming convention

**Test Command**:
```bash
nextflow run main.nf -profile test,singularity -entry MODULE1_VALIDATE \
  --sample_sheet test_data/sample_sheet.csv
```

**Verification**:
```bash
# Check output structure
ls -la results/module1/
# Verify variant counts match
bcftools stats results/module1/standardized/*.vcf.gz
```

---

### Module 2: Imputation Server Preparation

**Purpose**: QC and format data for imputation servers

**Tests**:
- [ ] **Pre-imputation QC**
  - [ ] Rayner/McCarthy HRC check runs
  - [ ] Strand flip correction applied
  - [ ] Allele frequency comparison works
  - [ ] Multi-allelic variants handled

- [ ] **Per-server formatting**
  - [ ] Michigan format correct
  - [ ] TOPMed format correct
  - [ ] NYGC format correct
  - [ ] All of Us format correct

- [ ] **Output validation**
  - [ ] VCF files are valid
  - [ ] Chromosome naming correct for each server
  - [ ] Files ready for upload

**Test Command**:
```bash
nextflow run main.nf -profile test,singularity -entry MODULE2_PREP_IMPUTATION \
  --imputation_servers "michigan,topmed"
```

**Verification**:
```bash
# Validate VCF for Michigan server
bcftools view -h results/module2/michigan/chr22.vcf.gz | head
# Check variant counts
bcftools stats results/module2/*/chr22.vcf.gz
```

---

### Module 3: Imputation Server Submission

**Purpose**: Submit to and download from imputation servers

**Tests**:
- [ ] **Server connectivity**
  - [ ] Michigan API authentication
  - [ ] TOPMed API authentication
  - [ ] NYGC/All of Us (manual process documented)

- [ ] **Job submission**
  - [ ] Files upload successfully
  - [ ] Job parameters set correctly
  - [ ] Job ID captured

- [ ] **Status monitoring**
  - [ ] Polls at correct interval
  - [ ] Handles server errors gracefully
  - [ ] Logs progress

- [ ] **Result download**
  - [ ] Downloads all result files
  - [ ] Decryption works (for encrypted results)
  - [ ] Verifies file integrity

**Test Command**:
```bash
# Dry run first
nextflow run main.nf -profile test,singularity -entry MODULE3_SUBMIT \
  --dry_run true --imputation_servers "michigan"
```

---

### Module 4: Post-Imputation QC

**Purpose**: Quality control of imputed data

**Tests**:
- [ ] **Quality filtering**
  - [ ] Rsq/INFO filtering works
  - [ ] MAF filtering applied
  - [ ] Configurable thresholds

- [ ] **Metrics calculation**
  - [ ] Per-variant Rsq extracted
  - [ ] MAF calculated
  - [ ] Sample-level metrics

- [ ] **MagicalRsq-X integration** (optional)
  - [ ] Model files loaded
  - [ ] Calibrated Rsq calculated
  - [ ] Comparison plots generated

- [ ] **QC reports**
  - [ ] HTML report generated
  - [ ] Variant counts by quality
  - [ ] MAF distribution plots

**Test Command**:
```bash
nextflow run main.nf -profile test,singularity -entry MODULE4_POST_QC \
  --rsq_threshold 0.3 --maf_filter 0.01
```

---

### Module 5: Server Result Comparison & Re-imputation

**Purpose**: Compare results across servers, identify failures, re-impute

**Tests**:
- [ ] **Cross-server comparison**
  - [ ] Variant overlap calculated
  - [ ] Concordance metrics computed
  - [ ] Server-specific variants identified

- [ ] **Region identification**
  - [ ] Low-quality regions detected
  - [ ] Failed samples identified
  - [ ] Ancestry-specific issues flagged

- [ ] **Re-imputation workflow**
  - [ ] Problem regions extracted
  - [ ] Secondary server submission works
  - [ ] Results merged correctly

**Test Command**:
```bash
nextflow run main.nf -profile test,singularity -entry MODULE5_COMPARE \
  --servers_to_compare "michigan,topmed"
```

---

### Module 6: Merge & Final QC

**Purpose**: Merge server results, final quality control

**Tests**:
- [ ] **Data merging**
  - [ ] Best-server-per-region logic works
  - [ ] Consensus calling correct
  - [ ] No duplicate variants

- [ ] **Final QC**
  - [ ] HWE filtering
  - [ ] Missingness filtering
  - [ ] Population-specific QC

- [ ] **Output formats**
  - [ ] PLINK format correct
  - [ ] VCF format correct
  - [ ] Dosage preserved

**Test Command**:
```bash
nextflow run main.nf -profile test,singularity -entry MODULE6_MERGE
```

---

### Module 7: Ancestry Estimation

**Purpose**: Global and local ancestry inference

**Tests**:
- [ ] **GRAF-anc**
  - [ ] Runs successfully
  - [ ] Ancestry proportions calculated
  - [ ] Classification output correct

- [ ] **ADMIXTURE**
  - [ ] Multiple K values tested
  - [ ] Cross-validation works
  - [ ] Optimal K selected

- [ ] **Local Ancestry (if enabled)**
  - [ ] **RFMix v2** runs correctly
  - [ ] **FLARE** runs correctly (if selected)
  - [ ] **G-NOMIX** runs correctly (if selected)
  - [ ] **RFMix v1** runs correctly (if selected)
  - [ ] Only selected method(s) run

- [ ] **Visualization**
  - [ ] PCA plots generated
  - [ ] ADMIXTURE bar plots
  - [ ] LAI karyograms (if LAI enabled)

**Test Command**:
```bash
nextflow run main.nf -profile test,singularity -entry MODULE7_ANCESTRY \
  --run_lai true --lai_methods "rfmix2"
```

---

### Integration Testing

**Full Pipeline Run**:
```bash
# Small test dataset
nextflow run main.nf -profile test,singularity \
  --sample_sheet test_data/sample_sheet.csv \
  --outdir results/integration_test \
  -resume

# Check all outputs
ls -laR results/integration_test/
```

**End-to-End Validation**:
- [ ] Pipeline completes without errors
- [ ] All expected output files present
- [ ] Logs capture all steps
- [ ] Resource usage within limits
- [ ] Results reproducible with `-resume`

---

### Troubleshooting Checklist

**Common Issues**:

| Symptom | Likely Cause | Solution |
|---------|--------------|----------|
| Container not found | Missing .sif file | Build with Module 0 or download |
| Liftover fails | Missing chain file | Download from UCSC |
| Rayner check fails | Wrong reference format | Use BRAVO→Rayner converter |
| GRAF-anc error | Missing data files | Download AncSnpPopAFs.txt |
| LAI not running | Wrong `lai_methods` format | Use comma-separated, no spaces |
| Memory error | Process limit too low | Increase in nextflow.config |
| API timeout | Server busy | Increase `monitor_interval` |

**Debug Commands**:
```bash
# View detailed error
nextflow log <run_name> -f status,hash,name,exit,error

# Check specific task
cat work/<hash>/.command.log
cat work/<hash>/.command.err

# Resume from failure
nextflow run main.nf -resume
```

---

## 3. Benchmarking Workflow

### Overview

Compare your pipeline against:
1. **Michigan Imputation Server** (baseline)
2. **TOPMed Imputation Server** (multi-ancestry reference)
3. **nf-core/phaseimpute** (Nextflow competitor)
4. **h3abionet/chipimputation** (African-focused)
5. **eQTL-Catalogue/genimpute** (production pipeline)

### Benchmarking Data Sources

#### Public Datasets with Paired Array + WGS

| Dataset | Ancestry | N Samples | Array | WGS | Access |
|---------|----------|-----------|-------|-----|--------|
| 1000 Genomes | Global | ~2,500 | Omni 2.5M, Affy 6.0 | 30x WGS | Public |
| GIAB | Mixed | 7 | Multiple | 50x+ WGS | Public |
| H3Africa | African | ~5,000 | H3Africa array | WGS subset | Controlled |
| gnomAD | Global | ~76,000 | - | WGS | Public (aggregates) |
| UK Biobank | European | 500,000 | UK Biobank array | 50k WGS | Controlled |

**Recommended Test Set**:
```
1. 1000 Genomes Omni 2.5M + WGS
   - 5 ancestries × 100 samples each = 500 samples
   - Ancestries: EUR (CEU), AFR (YRI), EAS (CHB), SAS (GIH), AMR (PUR)

2. GIAB samples for high-accuracy validation
   - HG001-HG007 with array genotypes
```

### Metrics to Measure

#### 1. Imputation Accuracy

```
Primary Metrics:
├── Dosage r² (squared correlation with WGS truth)
│   ├── Overall
│   ├── By MAF bin: <0.5%, 0.5-1%, 1-5%, 5-10%, 10-50%, >50%
│   └── By ancestry: EUR, AFR, EAS, SAS, AMR
├── Concordance rate (genotype match %)
├── Non-reference concordance
└── False positive rate (for rare variants)

Secondary Metrics:
├── Software Rsq/INFO score calibration
├── MagicalRsq vs true r²
└── Imputation certainty (mean GP)
```

#### 2. Coverage Metrics

```
├── Total variants imputed
├── Variants passing QC (Rsq ≥ 0.3, 0.8)
├── Overlap with WGS truth set
├── Unique variants per server
└── Failed regions count
```

#### 3. Computational Performance

```
├── Wall-clock time
│   ├── Per module
│   ├── Per chromosome
│   └── Total pipeline
├── Peak memory usage
├── CPU hours
├── Storage requirements
└── Cost (if cloud-based)
```

#### 4. Convenience/Usability

```
├── Setup complexity (time to first run)
├── Documentation quality
├── Error handling
├── Resume capability
├── Multi-server support
└── Ancestry-aware features
```

### Benchmarking Protocol

#### Phase 1: Data Preparation (Week 1)

```bash
# 1. Download 1000 Genomes data
./helper_scripts/create_test_data.sh

# 2. Create array-like subsets
plink2 --vcf 1kg_wgs.vcf.gz \
  --extract omni25_snplist.txt \
  --make-pgen \
  --out 1kg_omni25_subset

# 3. Split by ancestry
for pop in CEU YRI CHB GIH PUR; do
  plink2 --pfile 1kg_omni25_subset \
    --keep ${pop}_samples.txt \
    --make-pgen \
    --out 1kg_${pop}_omni25
done
```

#### Phase 2: Run Pipelines (Week 2-3)

```bash
# Your pipeline
nextflow run main.nf \
  --sample_sheet benchmarking/sample_sheet.csv \
  --imputation_servers "michigan,topmed" \
  --run_lai true \
  --outdir results/your_pipeline \
  -with-report benchmarking/your_pipeline_report.html \
  -with-trace benchmarking/your_pipeline_trace.txt \
  -with-timeline benchmarking/your_pipeline_timeline.html

# nf-core/phaseimpute
nextflow run nf-core/phaseimpute \
  --input benchmarking/samplesheet_nfcore.csv \
  --outdir results/nfcore_phaseimpute \
  -profile singularity

# Manual submission to Michigan/TOPMed servers
# (for baseline comparison)
```

#### Phase 3: Calculate Metrics (Week 4)

```python
# benchmarking/calculate_metrics.py
import pandas as pd
import numpy as np
from scipy.stats import pearsonr

def calculate_dosage_r2(imputed_dosage, true_genotype):
    """Calculate squared Pearson correlation."""
    mask = ~np.isnan(imputed_dosage) & ~np.isnan(true_genotype)
    if mask.sum() < 10:
        return np.nan
    r, _ = pearsonr(imputed_dosage[mask], true_genotype[mask])
    return r ** 2

def calculate_concordance(imputed_gt, true_gt):
    """Calculate genotype concordance rate."""
    mask = (imputed_gt >= 0) & (true_gt >= 0)
    if mask.sum() == 0:
        return np.nan
    return (imputed_gt[mask] == true_gt[mask]).mean()

def stratify_by_maf(variants_df, maf_bins):
    """Stratify metrics by MAF bin."""
    results = {}
    for low, high in maf_bins:
        mask = (variants_df['maf'] >= low) & (variants_df['maf'] < high)
        subset = variants_df[mask]
        results[f'{low}-{high}'] = {
            'n_variants': len(subset),
            'mean_r2': subset['dosage_r2'].mean(),
            'mean_concordance': subset['concordance'].mean()
        }
    return results

def stratify_by_ancestry(samples_df, ancestry_col='population'):
    """Stratify metrics by ancestry group."""
    return samples_df.groupby(ancestry_col).agg({
        'mean_r2': 'mean',
        'n_well_imputed': 'sum',
        'total_variants': 'sum'
    })
```

#### Phase 4: Generate Visualizations (Week 5)

```python
# benchmarking/plot_results.py
import matplotlib.pyplot as plt
import seaborn as sns

def plot_r2_by_maf(results_dict, output_file):
    """Plot imputation r² across MAF bins for multiple pipelines."""
    fig, ax = plt.subplots(figsize=(10, 6))

    maf_bins = ['<0.5%', '0.5-1%', '1-5%', '5-10%', '10-50%', '>50%']
    x = np.arange(len(maf_bins))
    width = 0.15

    for i, (pipeline, data) in enumerate(results_dict.items()):
        ax.bar(x + i*width, data['r2_values'], width, label=pipeline)

    ax.set_xlabel('Minor Allele Frequency')
    ax.set_ylabel('Mean Dosage r²')
    ax.set_title('Imputation Accuracy by MAF')
    ax.set_xticks(x + width * len(results_dict) / 2)
    ax.set_xticklabels(maf_bins)
    ax.legend()
    ax.set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)

def plot_r2_by_ancestry(results_dict, output_file):
    """Plot imputation r² across ancestry groups."""
    # Similar structure to above
    pass

def plot_computational_performance(trace_files, output_file):
    """Compare runtime and memory across pipelines."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Runtime comparison
    # Memory comparison

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)

def plot_coverage_venn(pipeline_variants, output_file):
    """Venn diagram of variants imputed by each pipeline."""
    from matplotlib_venn import venn3
    # Implementation
    pass
```

### Expected Results Summary Table

| Metric | Your Pipeline | Michigan | TOPMed | nf-core |
|--------|---------------|----------|--------|---------|
| **Accuracy** |
| Overall r² (MAF>5%) | | | | |
| Overall r² (MAF 1-5%) | | | | |
| Overall r² (MAF<1%) | | | | |
| EUR r² | | | | |
| AFR r² | | | | |
| AMR r² | | | | |
| **Coverage** |
| Total variants | | | | |
| Well-imputed (r²≥0.8) | | | | |
| Unique variants | | | | |
| **Performance** |
| Total runtime | | | | |
| Peak memory | | | | |
| CPU hours | | | | |
| **Features** |
| Multi-server | ✓ | - | - | - |
| LAI integration | ✓ | - | - | - |
| Auto-ancestry QC | ✓ | - | - | - |
| Resume support | ✓ | - | - | ✓ |

---

## 4. Comparative Analysis Against Existing Pipelines

### Key Differentiators of Your Pipeline

| Feature | Your Pipeline | nf-core/phaseimpute | h3abionet | eQTL-Catalogue |
|---------|---------------|---------------------|-----------|----------------|
| Multi-server submission | ✓ | - | - | - |
| Server result comparison | ✓ | - | - | - |
| Best-per-region merging | ✓ | - | - | - |
| Integrated LAI | ✓ | - | - | - |
| MagicalRsq-X support | ✓ | - | - | - |
| Ancestry-stratified QC | ✓ | - | ✓ | - |
| Auto re-imputation | ✓ | - | - | - |
| All of Us support | ✓ | - | - | - |

### Unique Selling Points

1. **Multi-Server Orchestration**: Only pipeline that submits to multiple servers and intelligently merges results
2. **Trans-Ancestry Optimization**: Designed for diverse populations with ancestry-aware QC
3. **Integrated Local Ancestry**: RFMix v2, FLARE, G-NOMIX built-in
4. **Quality Calibration**: MagicalRsq-X integration for better accuracy estimation
5. **Comprehensive Module System**: 7 specialized modules vs monolithic workflows

### Comparison Points for Publication

```
Accuracy Advantages:
├── Better rare variant imputation through multi-server consensus
├── Improved accuracy for admixed populations via LAI
├── Calibrated quality scores with MagicalRsq-X
└── Ancestry-specific server selection

Usability Advantages:
├── Single workflow for all major imputation servers
├── Automatic format conversion and liftover
├── Built-in troubleshooting (re-imputation of failed regions)
└── Comprehensive ancestry analysis

Reproducibility Advantages:
├── Nextflow DSL2 with containerization
├── Version-controlled containers
├── Detailed provenance tracking
└── Resume capability
```

---

## 5. Publication Strategy

### Target Journals

#### Tier 1: High-Impact Methods Journals

| Journal | IF | Focus | Fit |
|---------|-----|-------|-----|
| **Nature Methods** | 48 | Novel methods | If showing significant accuracy improvement |
| **Genome Biology** | 17 | Genomics methods | Good fit for comprehensive pipeline |
| **Genome Research** | 9 | Computational genomics | Strong fit |

#### Tier 2: Bioinformatics/Methods Journals

| Journal | IF | Focus | Fit |
|---------|-----|-------|-----|
| **Bioinformatics** | 6.9 | Software/methods | Excellent fit for pipeline paper |
| **NAR Genomics & Bioinformatics** | 4.8 | Genomics software | Good fit |
| **GigaScience** | 7.7 | Data-intensive research | Good for benchmark-heavy paper |
| **BMC Bioinformatics** | 3.3 | Methods | Good fallback |

#### Tier 3: Domain-Specific Journals

| Journal | IF | Focus | Fit |
|---------|-----|-------|-----|
| **AJHG** | 11 | Human genetics | If ancestry focus emphasized |
| **Genetic Epidemiology** | 2.5 | Methods in genetics | Good for methods paper |
| **Frontiers in Genetics** | 3.7 | Open access genetics | Good option |

### Recommended Strategy

**Primary Target**: **Bioinformatics** (Applications Note or Original Paper)
- Strong track record for pipeline papers
- Appropriate scope
- Reasonable review timeline

**Backup**: **GigaScience** or **NAR Genomics & Bioinformatics**
- Open access
- Appreciate comprehensive benchmarking
- Support data sharing

### Paper Structure Outline

```
Title: "A Multi-Server Genotype Imputation Pipeline with Integrated
       Ancestry Analysis for Trans-Ancestry Genetic Studies"

Abstract (250 words)
- Problem: Imputation accuracy varies by ancestry, single-server approaches suboptimal
- Solution: Multi-server orchestration with intelligent merging
- Results: X% improvement in rare variant imputation, Y% better for non-EUR
- Availability: GitHub, containers, documentation

Introduction
- Imputation importance for GWAS
- Current limitations (single server, ancestry bias)
- Need for trans-ancestry approaches

Methods
- Pipeline architecture (7 modules)
- Multi-server submission and comparison
- Intelligent result merging
- LAI integration
- Benchmarking design

Results
- Accuracy comparison vs existing pipelines
- MAF-stratified performance
- Ancestry-stratified performance
- Computational benchmarks
- Case study: Trans-ancestry imputation

Discussion
- Advantages for diverse populations
- Limitations
- Future directions

Availability
- GitHub repository
- Container registry
- Documentation
- Test data
```

### Required for Publication

- [ ] **Code availability**: GitHub with DOI (Zenodo)
- [ ] **Container availability**: Docker Hub or Singularity Hub
- [ ] **Documentation**: ReadTheDocs or GitHub Pages
- [ ] **Test data**: Publicly accessible
- [ ] **Benchmarking data**: Reproducible scripts
- [ ] **Example outputs**: Supplementary materials

### Timeline

| Week | Task |
|------|------|
| 1-2 | Fix critical issues, prepare test data |
| 3-4 | Run benchmarking experiments |
| 5-6 | Calculate metrics, generate figures |
| 7-8 | Write manuscript draft |
| 9 | Internal review, revisions |
| 10 | Submit to Bioinformatics |

---

## Appendix A: Quick Reference Commands

```bash
# Full pipeline test
nextflow run main.nf -profile test,singularity --outdir test_results

# Single module test
nextflow run main.nf -profile test,singularity -entry MODULE1_VALIDATE

# Resume failed run
nextflow run main.nf -resume

# Generate execution report
nextflow run main.nf -with-report report.html -with-trace trace.txt

# Clean work directory
nextflow clean -f

# View run history
nextflow log
```

## Appendix B: File Checklist

```
Required before testing:
├── config/nextflow.config (parameters updated)
├── resources/
│   ├── containers/*.sif (built)
│   ├── references/
│   │   ├── hg38.fa
│   │   └── hg19ToHg38.over.chain.gz
│   ├── ancestry_references/
│   │   ├── lai_reference/ (populated)
│   │   └── grafanc_data/AncSnpPopAFs.txt
│   └── genetic_maps/ (populated)
└── test_data/
    ├── sample_sheet.csv
    └── input_files/
```
