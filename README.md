# Genotyping Imputation & Ancestry Pipeline

A comprehensive Nextflow DSL2 pipeline for genotyping array data quality control, imputation, platform merging, and ancestry estimation.

---

## Overview

This pipeline processes genotyping array data through seven integrated modules:

```
                                    INPUT
                                      │
                                      ▼
    ┌─────────────────────────────────────────────────────────────┐
    │  MODULE 1: Pre-Imputation QC                                │
    │  • Sample discovery and format conversion                   │
    │  • Reference alignment (hg19→hg38 liftover if needed)      │
    │  • Platform-level merging (UNION)                          │
    │  • Light QC (call rates, duplicates, monomorphic)          │
    │  • Create service-specific VCFs                            │
    └─────────────────────────────────────────────────────────────┘
                          │                    │
                    TOPMed VCFs          AnVIL VCF
                          │                    │
                          ▼                    ▼
    ┌─────────────────────────────────────────────────────────────┐
    │  MODULE 2: Imputation                                       │
    │  • TOPMed Imputation Server (imputationbot)                │
    │  • All of Us AnVIL Service (terralab-cli)                  │
    │  • Automatic job monitoring and download                   │
    └─────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
    ┌─────────────────────────────────────────────────────────────┐
    │  MODULE 3: Post-Imputation QC                               │
    │  • MagicalRsq-X imputation quality filtering               │
    │  • Sample and variant call rate filters                    │
    │  • Per-platform QC reports                                 │
    └─────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
    ┌─────────────────────────────────────────────────────────────┐
    │  MODULE 4: Platform Merging                                 │
    │  • Cross-platform variant harmonization                    │
    │  • INTERSECTION merge (shared variants only)               │
    │  • Multi-allelic resolution                                │
    └─────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
    ┌─────────────────────────────────────────────────────────────┐
    │  MODULE 5: Re-Imputation (Optional)                         │
    │  • Re-impute merged dataset                                │
    │  • Fill gaps from platform merging                         │
    │  • Same services as Module 2                               │
    └─────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
    ┌─────────────────────────────────────────────────────────────┐
    │  MODULE 6: Post-Merge QC (Final)                            │
    │  • MagicalRsq-X second-pass filtering                      │
    │  • HWE, MAF, heterozygosity filters                        │
    │  • Relatedness detection (GENESIS/PLINK)                   │
    │  • PCA for population structure                            │
    └─────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
    ┌─────────────────────────────────────────────────────────────┐
    │  MODULE 7: Ancestry Estimation                              │
    │  • GRAF-anc global ancestry (8 major groups + subgroups)   │
    │  • ADMIXTURE (K=2-12 with cross-validation)                │
    │  • Local Ancestry Inference (RFMix v2, optional others)    │
    │  • Comprehensive visualizations                            │
    └─────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
                                   OUTPUT
                      (Analysis-ready datasets + ancestry)
```

---

## Quick Start

### 1. Prerequisites

- **Nextflow** >= 23.04.0
- **Apptainer/Singularity** (for containerized execution)
- **Python 3.12+** (for terralab-cli if using All of Us)

### 2. Clone and Setup

```bash
git clone https://github.com/your-repo/genotyping_pipeline.git
cd genotyping_pipeline

# Build containers (Module 0)
nextflow run modules/Module0_Apptainer_Build.nf -profile apptainer
```

### 3. Prepare Reference Files

Download and place reference files in `resources/references/`:
- `hg19.fa` + `hg19.fa.fai` - Human reference genome GRCh37
- `hg38.fa` + `hg38.fa.fai` - Human reference genome GRCh38
- `hg19ToHg38.over.chain.gz` - Liftover chain file
- `PASS.Variants.TOPMed_freeze10_hg38.tab.gz` - TOPMed reference

See `resources/references/README.md` for download links.

### 4. Create Sample Sheet

Create a CSV file with your input data:

```csv
platform_id,batch_id,input_path,file_type,build,file_structure
GSAv1,batch_2020_01,/data/gsa_v1/batch_2020_01,plink,hg19,individual_samples
GSAv2,batch_2021_03,/data/gsa_v2/batch_2021_03,plink,hg38,merged_batch
```

See `documentation/input_file_set_up.md` for detailed examples.

### 5. Configure Imputation Services

#### TOPMed Imputation Server
1. Create account at https://imputation.biodatacatalyst.nhlbi.nih.gov
2. Get your API token from the web interface
3. Set in config or command line: `--topmed_api_token YOUR_TOKEN`

#### All of Us / AnVIL Service (Optional)
1. Install terralab-cli: `pip install terralab-cli`
2. Authenticate: `terralab login` (opens browser)
3. Enable in config: `run_anvil = true`

### 6. Run Pipeline

```bash
# Basic run with TOPMed only
nextflow run main.nf \
    --sample_sheet samples.csv \
    --topmed_api_token YOUR_TOKEN \
    --topmed_password YOUR_PASSWORD \
    -profile apptainer,slurm

# Run with both imputation services
nextflow run main.nf \
    --sample_sheet samples.csv \
    --topmed_api_token YOUR_TOKEN \
    --topmed_password YOUR_PASSWORD \
    --run_anvil true \
    -profile apptainer,slurm

# Skip optional modules
nextflow run main.nf \
    --sample_sheet samples.csv \
    --skip_modules "5,7" \
    -profile apptainer
```

---

## Configuration

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sample_sheet` | required | Input CSV file |
| `outdir` | `results` | Output directory |
| `run_topmed` | `true` | Enable TOPMed imputation |
| `run_anvil` | `false` | Enable All of Us imputation |
| `skip_modules` | `""` | Comma-separated modules to skip |

### QC Thresholds

| Parameter | Default | Description |
|-----------|---------|-------------|
| `magicalrsq_threshold` | `0.3` | MagicalRsq-X R² threshold |
| `sample_call_rate` | `0.95` | Minimum sample call rate |
| `variant_call_rate` | `0.95` | Minimum variant call rate |
| `hwe_pvalue` | `1e-6` | HWE filter p-value |
| `kinship_threshold` | `0.177` | 1st-degree relatedness threshold |

### Ancestry Analysis

| Parameter | Default | Description |
|-----------|---------|-------------|
| `admixture_k` | `"5,6,7"` | ADMIXTURE K values |
| `run_lai` | `false` | Enable local ancestry inference |
| `lai_methods` | `"rfmix2"` | LAI methods (rfmix2,rfmix1,flare,gnomix) |

---

## Execution Profiles

```bash
# Local execution with Apptainer
nextflow run main.nf -profile apptainer

# SLURM cluster
nextflow run main.nf -profile apptainer,slurm

# PBS cluster
nextflow run main.nf -profile apptainer,pbs

# Test mode (reduced resources)
nextflow run main.nf -profile apptainer,test
```

---

## Output Structure

```
results/
├── module1/                    # Pre-imputation QC outputs
│   ├── 01_discovery/
│   ├── 02_prep/
│   ├── 03_merge/
│   ├── 04_alignment/
│   ├── 05_liftover/
│   ├── 06_platform_merge/
│   ├── 07_validation/
│   ├── 08_light_qc/
│   └── 09_vcf_creation/
├── module2/                    # Imputation outputs
│   ├── 01_topmed/
│   └── 02_anvil/
├── module3/                    # Post-imputation QC
├── module4/                    # Platform merging
├── module5/                    # Re-imputation (optional)
├── module6/                    # Final QC
│   ├── maf_filtered/
│   ├── no_maf_filter/
│   └── pca/
├── module7/                    # Ancestry estimation
│   ├── 01_grafanc/
│   ├── 02_admixture/
│   └── 03_local_ancestry/
└── pipeline_info/              # Execution reports
    ├── execution_timeline.html
    ├── execution_report.html
    └── execution_trace.txt
```

---

## Module Descriptions

### Module 0: Container Build
Builds all required Apptainer/Singularity containers from definition files.

### Module 1: Pre-Imputation QC
- Discovers and loads input files (4 structure types supported)
- Merges samples to batches, batches to platforms (UNION merge)
- Aligns to reference genome (bcftools +fixref)
- Lifts hg19 to hg38 (CrossMap)
- Applies light QC (call rates, duplicates)
- Creates service-specific VCFs

### Module 2: Imputation
- Submits to TOPMed and/or All of Us imputation services
- Monitors job progress automatically
- Downloads and organizes results
- Handles API rate limits and retries

### Module 3: Post-Imputation QC
- Applies MagicalRsq-X imputation quality filtering
- Filters by sample and variant call rates
- Generates per-platform QC reports

### Module 4: Platform Merging
- Merges multiple platforms (INTERSECTION)
- Resolves multi-allelic variants
- Harmonizes variant IDs

### Module 5: Re-Imputation (Optional)
- Re-imputes merged dataset
- Uses same services as Module 2
- Fills gaps from intersection merging

### Module 6: Post-Merge QC
- Second-pass MagicalRsq-X filtering
- HWE, MAF, heterozygosity filters
- Relatedness detection (GENESIS PCRelate or PLINK)
- PCA for population structure
- Creates MAF-filtered and unfiltered outputs

### Module 7: Ancestry Estimation
- **GRAF-anc**: Global ancestry (8 major groups, 50+ subgroups)
- **ADMIXTURE**: Unsupervised clustering (K=2-12)
- **Local Ancestry** (optional): RFMix v2, RFMix v1, FLARE, G-NOMIX
- Comprehensive visualizations and reports

---

## Documentation

- `documentation/input_file_set_up.md` - Sample sheet format guide
- `documentation/Module_1_Workflow.md` - Detailed Module 1 workflow
- `documentation/pipeline_workflow.md` - Complete pipeline workflow

---

## Requirements

### Software (containerized)
- PLINK 1.9 / PLINK 2.0
- bcftools 1.18+
- CrossMap
- Python 3.11+
- R 4.x with genetics packages
- imputationbot (TOPMed CLI)
- terralab-cli (All of Us CLI)
- GRAF-anc
- ADMIXTURE
- RFMix v2

### Resources
- **Storage**: ~50GB per 1000 samples (varies by imputation density)
- **Memory**: 4-64GB depending on sample size
- **Time**: Hours to days depending on sample count and imputation queue

---

## Citation

If you use this pipeline, please cite:
- TOPMed Imputation Server
- All of Us Research Program (if used)
- MagicalRsq-X
- GRAF-anc
- ADMIXTURE
- RFMix (if LAI used)

---

## Support

For issues and feature requests, please open an issue on GitHub.

---

## License

[Your License Here]
