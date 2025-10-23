

# Genotyping Imputation and Ancestry Pipeline
### Complete Apptainer-Based Pipeline for HPC and Cloud Environments

## 📋 Overview

This Nextflow pipeline performs comprehensive genotyping data QC, imputation, merging, and ancestry estimation across multiple array platforms. It's designed to work seamlessly across HPC clusters (SLURM/PBS), cloud environments (Google Cloud, AWS), and research biobanks (All of Us, UK Biobank).

### Pipeline Modules

1. **Module 1**: Pre-imputation QC and file preparation
2. **Module 2**: Imputation via TOPMed and/or All of Us AnVIL (parallelized by platform)
3. **Module 3**: Post-imputation QC with MagicalRsq
4. **Module 4**: Platform merging with QC
5. **Module 5**: Re-imputation of merged data
6. **Module 6**: Post-merge QC (2nd pass) with GENESIS/PLINK
7. **Module 7**: Ancestry estimation (global and local ancestry inference)

---

## 🐳 Apptainer Containers

This pipeline uses **5 Apptainer containers** for reproducibility and portability:

### Container Summary

| Container | Tools Included | Used in Modules |
|-----------|---------------|-----------------|
| **plink.sif** | PLINK 1.9, PLINK 2.0 | 1, 3, 4, 6, 7 |
| **r_genetics.sif** | R 4.3.1, GENESIS, MagicalRsq, tidyverse, ggplot2 | 3, 6, 7 |
| **perl_vcftools.sif** | Perl, bcftools, vcftools, CrossMap, htslib | 1, 3, 4, 5, 6 |
| **ancestry_tools.sif** | ADMIXTURE, RFMix v1/v2, FLARE, g-nomix, Graf-anc, SHAPEIT4, Eagle, TRACTOR | 7 |
| **python_tools.sif** | Python 3.10, pandas, pysam, API libraries, terra-notebook-utils | 2, 5 |

### Software by Module

**Module 1** (Pre-imputation QC):
- PLINK 1.9/2.0 (basic QC, frequency calculations)
- Perl (Will Rayner strand checking scripts)
- CrossMap (liftover from hg19 to hg38)
- bcftools/vcftools (VCF manipulation)

**Module 2** (Imputation):
- Python (TOPMed API, All of Us AnVIL API)
- bcftools (VCF format conversion)

**Module 3** (Post-imputation QC):
- R (MagicalRsq filtering)
- bcftools (VCF filtering, call rate calculation)

**Module 4** (Merging):
- PLINK 1.9/2.0 (merging platforms)
- bcftools (VCF merging)

**Module 5** (Re-imputation):
- Python (API interactions)
- bcftools (file preparation)

**Module 6** (Post-merge QC):
- R with GENESIS (relatedness analysis, PCA)
- PLINK 1.9/2.0 (MAF filtering, HWE, heterozygosity)
- bcftools (VCF to PLINK conversion)

**Module 7** (Ancestry):
- **Global Ancestry**: ADMIXTURE, Graf-anc/GrafPop
- **Local Ancestry**: RFMix v1, RFMix v2, FLARE, g-nomix, TRACTOR
- **Phasing**: SHAPEIT4, Eagle
- **Analysis**: R (visualization, summary statistics)

---

## 📦 Repository Structure

genotyping-imputation-pipeline/
├── README.md
├── LICENSE
├── .gitignore
│
├── main.nf                          # Main workflow
├── nextflow.config                  # Configuration with Apptainer/Singularity
│
├── modules/                         # All .nf module files
│   ├── Module1_PreImputation.nf
│   ├── Module2_Imputation.nf       # DUAL WORKFLOW VERSION
│   ├── Module3_PostQC.nf           # DUAL WORKFLOW VERSION
│   ├── Module4_Merging.nf          # DUAL WORKFLOW VERSION
│   ├── Module5_Reimputation.nf     # DUAL WORKFLOW VERSION
│   ├── Module6_PostMergeQC.nf      # DUAL WORKFLOW VERSION
│   └── Module7_Ancestry.nf         # DUAL WORKFLOW VERSION
│
├── bin/                            # ⭐ Helper scripts
│   ├── check-topmed-strands.pl     # Will Rayner's HRC-1000G-check-bim.pl
│   ├── magicalrsqx_filter.R        # MagicalRsq-X filtering
│   ├── compare_imputation.R        # Compare TOPMed vs AnVIL
│   ├── pcrelate_qc.R               # GENESIS relatedness
│   ├── het_filter.R                # Heterozygosity filtering
│   └── imputation_api.py           # API helper functions
│
├── containers/                      # ⭐ Apptainer definitions
│   ├── plink.def
│   ├── r_analysis.def              # With MagicalRsq-X
│   ├── bcftools.def
│   ├── python_api.def
│   ├── admixture.def
│   └── rfmix.def
│
├── resources/                       # ⭐ Reference data
│   ├── references/                  # Main references directory
│   │   ├── hg38.fa                 # Human reference (download separately)
│   │   ├── hg38.fa.fai
│   │   ├── hg19ToHg38.over.chain.gz
│   │   ├── TOPMed_freq/            # TOPMed frequency files
│   │   │   └── PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab
│   │   ├── genetic_maps/           # For RFMix
│   │   │   ├── chr1.map
│   │   │   ├── chr2.map
│   │   │   └── ...
│   │   └── ancestry_panels/        # Ancestry reference panels
│   │       ├── admixture_ref/
│   │       └── rfmix_ref/
│   └── .gitkeep
│
├── submodules/                      # ⭐ Git submodules
│   └── graf-anc/                   # Graf-anc with its own references
│       └── (submodule with its own structure)
│
├── docs/                            # Documentation
│   ├── INSTALLATION.md
│   ├── USAGE.md
│   ├── API_SETUP.md
│   ├── DUAL_WORKFLOW_COMPARISON.md # NEW: Comparing TOPMed vs AnVIL
│   └── TROUBLESHOOTING.md
│
├── examples/                        # Example run scripts
│   ├── run_topmed_only.sh
│   ├── run_anvil_only.sh
│   ├── run_both_parallel.sh
│   └── sample_sheet.tsv
│
└── test/                            # Test data
    ├── test_data/
    └── run_test.sh

---

## 🚀 Quick Start

### 1. Clone Repository

```bash
git clone https://github.com/your-username/genotyping-imputation-pipeline.git
cd genotyping-imputation-pipeline
```

### 2. Build Apptainer Containers

```bash
cd containers
chmod +x build_containers.sh
./build_containers.sh
```

This creates 5 `.sif` files in `container_images/` directory.

### 3. Deploy Containers

**For HPC:**
```bash
# Copy to shared storage
cp container_images/*.sif /shared/containers/

# Update nextflow.config
# Change: container_path = '/shared/containers'
```

**For Google Cloud (All of Us):**
```bash
# Upload to Google Cloud Storage
gsutil -m cp container_images/*.sif gs://your-bucket/containers/

# Update nextflow.config
# Change: container_path = 'gs://your-bucket/containers'
```

**For AWS:**
```bash
# Upload to S3
aws s3 cp container_images/ s3://your-bucket/containers/ --recursive --include "*.sif"

# Update nextflow.config
# Change: container_path = 's3://your-bucket/containers'
```

### 4. Prepare Input Data

Create a tab-separated sample sheet (`platforms.tsv`):

```
platform1	/path/to/platform1.bed	/path/to/platform1.bim	/path/to/platform1.fam
platform2	/path/to/platform2.bed	/path/to/platform2.bim	/path/to/platform2.fam
```

### 5. Download Reference Files

```bash
mkdir -p references

# Download hg38 reference genome
wget -O references/hg38.fa.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip references/hg38.fa.gz

# Download liftOver chain file
wget -O references/hg19ToHg38.over.chain.gz \
  http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# Download genetic maps for RFMix
mkdir -p references/genetic_maps
# (See docs/SETUP_GUIDE.md for detailed instructions)
```

### 6. Run Pipeline

**On HPC (SLURM):**
```bash
nextflow run main.nf \
  --sample_sheet platforms.tsv \
  --run_topmed \
  --topmed_api_token $TOPMED_TOKEN \
  -profile slurm,apptainer \
  -resume
```

**On Google Cloud (All of Us AnVIL):**
```bash
nextflow run main.nf \
  --sample_sheet platforms.tsv \
  --run_anvil \
  --anvil_workspace my-workspace \
  --anvil_project my-gcp-project \
  -profile google,apptainer \
  -resume
```

**On PBS Cluster:**
```bash
nextflow run main.nf \
  --sample_sheet platforms.tsv \
  --run_topmed \
  --topmed_api_token $TOPMED_TOKEN \
  -profile pbs,apptainer \
  -resume
```

---

## 🔧 Configuration

### Key Parameters

Edit `nextflow.config` or use command-line flags:

```bash
# General
--input_build hg19         # Input genome build (hg19 or hg38)
--outdir results           # Output directory
--skip_modules "5,7"       # Skip specific modules

# Imputation Services
--run_topmed               # Use TOPMed imputation (default: true)
--run_anvil                # Use All of Us AnVIL (default: false)
--topmed_api_token TOKEN   # TOPMed API token
--anvil_workspace NAME     # AnVIL workspace name

# Module 6
--use_genesis              # Use GENESIS for relatedness (default: true)

# Module 7
--lai_methods "rfmix2,tractor"  # LAI methods to run
--admixture_k "5,6,7"      # ADMIXTURE K values
```

### Execution Profiles

Available profiles in `nextflow.config`:

- `apptainer` - Use Apptainer containers (recommended)
- `docker` - Use Docker containers (local testing)
- `slurm` - SLURM cluster execution
- `pbs` - PBS/Torque cluster execution
- `google` - Google Cloud / All of Us AnVIL
- `aws` - AWS Batch execution
- `test` - Reduced resources for testing

**Combine profiles with commas:**
```bash
-profile slurm,apptainer
-profile google,apptainer
```

---

## 📊 Output Structure

```
results/
├── module1/                 # Pre-imputation QC
│   ├── 01_basic_qc/        # Per-platform QC
│   ├── 02_strand_check/    # Strand alignment
│   ├── 03_liftover/        # Build conversion (if needed)
│   └── 04_formatted/       # Ready for imputation
├── module2/                 # Imputation results
│   ├── topmed/             # TOPMed results (if run_topmed=true)
│   └── anvil/              # AnVIL results (if run_anvil=true)
├── module3/                 # Post-imputation QC
│   ├── 01_rsq_filtered/    # MagicalRsq R²>=0.3
│   └── 02_sample_filtered/ # Call rate >= 95%
├── module4/                 # Merged platforms
│   ├── 01_per_chr/         # Per-chromosome merges
│   └── 02_merged_genome/   # Genome-wide merged file
├── module5/                 # Re-imputation
│   └── (same structure as module2)
├── module6/                 # Final QC
│   ├── 03_relatedness/     # GENESIS/PLINK relatedness
│   ├── 04_hwe_filter/      # HWE filtering
│   ├── 05_het_filter/      # Heterozygosity filtering
│   └── 06_final_datasets/  # ✨ ANALYSIS-READY DATA
│       ├── final_maf01.*   # MAF >= 0.01
│       └── final_nomaf.*   # No MAF filter
└── module7/                 # Ancestry
    ├── 01_grafpop/         # Global ancestry
    ├── 02_admixture/       # ADMIXTURE results
    └── 03_local_ancestry/  # RFMix, FLARE, g-nomix, TRACTOR
```

**Analysis-Ready Datasets:**
- `results/module6/06_final_datasets/final_maf01.*` - MAF filtered (recommended for most analyses)
- `results/module6/06_final_datasets/final_nomaf.*` - No MAF filter (for rare variant studies)

---

## 📖 Documentation

- **SETUP_GUIDE.md** - Detailed setup instructions
- **MODULE_DETAILS.md** - In-depth module documentation
- **TROUBLESHOOTING.md** - Common issues and solutions
- **CONTAINER_GUIDE.md** - Container building and management

---

## 🔬 Citation

If you use this pipeline in your research, please cite:

```
Your Paper Title
Authors
Journal, Year
DOI: ...
```

And cite the tools used:
- **PLINK**: Purcell S, et al. (2007)
- **GENESIS**: Conomos MP, et al. (2016)
- **MagicalRsq**: Zhang J, et al. (2024)
- **ADMIXTURE**: Alexander DH, et al. (2009)
- **RFMix**: Maples BK, et al. (2013); Marron DL, et al. (2020)
- **TOPMed Imputation Server**: Taliun D, et al. (2021)

---

## 🤝 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

---

## 📧 Contact

- **Author**: Your Name
- **Email**: your.email@institution.edu
- **Issues**: [GitHub Issues](https://github.com/your-username/genotyping-imputation-pipeline/issues)

---

## 📄 License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- TOPMed Imputation Server
- All of Us Research Program
- NCBI Graf-anc team
- Nextflow community
- Bioconda and BioContainers projects

