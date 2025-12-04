# Benchmarking Framework

Comprehensive benchmarking to compare imputation strategies across ancestry groups for mixed-ancestry cohorts.

## Overview

This framework compares **6 imputation approaches** across **3 imputation servers** using publicly available data with WGS truth for validation.

## Approaches Compared

### Traditional Approaches (Shell Scripts)

| Approach | QC Timing | Merge Strategy | Description |
|----------|-----------|----------------|-------------|
| **A** | Before imputation | Intersect before | Traditional: Thorough QC → Intersect → Impute → R² filter |
| **B** | After imputation | Intersect after | Minimal prep → Impute → Intersect → Thorough QC after |
| **C** | Before imputation | Separate platforms | Thorough QC → Impute each platform separately → Merge after |
| **D** | After imputation | Separate platforms | Minimal prep → Impute each separately → Thorough QC → Merge |

### Our Pipeline Variants

| Approach | QC Timing | Merge Strategy | Description |
|----------|-----------|----------------|-------------|
| **E (Ours-1step)** | After imputation | Union within platform | Union merge → Impute → Intersect merge → MagicalRsq-X → Thorough QC |
| **F (Ours-2step)** | After imputation | Union + re-impute | Union merge → Impute → Intersect merge → Re-impute → MagicalRsq-X → Thorough QC |

### Key Differences in Our Pipeline

```
TRADITIONAL (Approach A):
  Intersect across platforms → Thorough QC → Impute → R² filter (0.3 or 0.8)
  ❌ Loses platform-specific variants before imputation
  ❌ HWE/MAF filters remove valid admixed population variants
  ❌ Static R² threshold not calibrated for ancestry

OUR PIPELINE (Approach E/F):
  Union merge within platform → Minimal QC → Impute → Intersect merge →
  [Optional: Re-impute] → MagicalRsq-X filter → Thorough QC
  ✓ Preserves all platform variants until after imputation
  ✓ QC after imputation catches true errors, not population structure
  ✓ MagicalRsq-X provides ancestry-calibrated quality scores
  ✓ 2-step imputation recovers variants lost in intersect merge
```

## Imputation Servers Tested

| Server | Reference Panel | Notes |
|--------|-----------------|-------|
| Michigan Imputation Server | 1000 Genomes Phase 3 | Baseline comparison |
| TOPMed Imputation Server | TOPMed Freeze 10 | ~97K samples, diverse |
| All of Us AnVIL | All of Us WGS | ~245K samples, diverse |

## Test Data Sources

### Genotype Array Data (Input)
- **1000 Genomes** Omni 2.5M array data (public)
- Multiple ancestry groups: AFR, EUR, EAS, AMR, SAS
- ~2500 samples with matched WGS

### WGS Truth Data (Validation)
- **HGDP-1KG High-Coverage WGS** (NYGC, ~30x)
- **1000 Genomes 30x WGS** (public)
- Matched samples for concordance calculation

## Metrics Calculated

### 1. Variant-Level Imputation Quality

| Metric | Description | Stratification |
|--------|-------------|----------------|
| INFO/R² Score | Imputation accuracy estimate | MAF, ancestry, functional category |
| MagicalRsq-X | Ancestry-calibrated quality | Ancestry-specific models |
| Variant Retention | % variants passing QC | MAF bin, ancestry |
| Imputation Rate | % variants successfully imputed | Platform, ancestry |

### 2. Concordance with WGS Truth

| Metric | Description | Expected (Common) | Expected (Rare) |
|--------|-------------|-------------------|-----------------|
| Genotype Concordance | Exact GT match rate | >98% | >90% |
| Dosage R² | Correlation of dosage with truth | >0.95 | >0.70 |
| Non-Reference Concordance | Alt allele accuracy | >95% | >85% |
| Heterozygote Concordance | Het call accuracy | >95% | >80% |

### 3. Cross-Platform Reproducibility
- Genotype concordance for samples typed on multiple arrays
- Dosage correlation across platforms post-imputation
- Platform-specific variant bias detection

### 4. Association/PRS Performance
| Metric | Description | Stratification |
|--------|-------------|----------------|
| Effect Size Bias | Deviation from expected β | Ancestry, MAF |
| Lambda (λ GC) | Genomic inflation factor | Ancestry |
| Power | % known GWAS hits recovered | Ancestry, trait |
| PRS R² | Variance explained by PRS | Ancestry group |
| PRS AUC | Discriminative ability | Ancestry group |

### 5. Ancestry-Aware Metrics
| Metric | Description | Notes |
|--------|-------------|-------|
| Global Ancestry Quality | INFO by GRAF-anc group | 8 major + 50 subgroups |
| Local Ancestry Quality | INFO by LAI tract ancestry | Per-chromosome segments |
| Ancestry-Specific Allele Freq | AF bias by ancestry | Compare to gnomAD |
| LAI-Aware PRS | Tractor framework performance | Ancestry-specific effects |

### 6. Stratification Dimensions
- **MAF bins:** 0-0.1%, 0.1-0.5%, 0.5-1%, 1-5%, 5-50%
- **Global ancestry:** AFR, EUR, EAS, SAS, AMR, MEN, OCN, MIX (GRAF-anc)
- **Local ancestry:** AFR, EUR, NAT tracts (from RFMix2/FLARE)
- **Functional category:** Exonic, intronic, intergenic, GWAS catalog
- **Chromosome:** Per-chromosome consistency check

### 7. Computational Metrics
| Metric | Description |
|--------|-------------|
| Wall-clock time | Total runtime |
| CPU hours | Computational cost |
| Peak memory | Maximum RAM usage |
| Storage | Disk space required |
| I/O operations | Read/write load |

## Publication Figures

### Main Figures

| Figure | Content | Type |
|--------|---------|------|
| **Fig 1** | Pipeline overview flowchart | Schematic |
| **Fig 2** | Concordance by MAF × Ancestry (heatmap) | Heatmap |
| **Fig 3** | INFO/R² distributions by approach | Violin/Box |
| **Fig 4A** | Local ancestry–specific imputation quality | Heatmap (ancestry × MAF) |
| **Fig 4B** | INFO difference by LAI tracts | Violin plot |
| **Fig 4C** | PRS R²/AUC by ancestry group | Barplot |
| **Fig 5A** | GWAS Manhattan plots comparison | Manhattan |
| **Fig 5B** | QQ plots showing inflation | QQ plot |
| **Fig 5C** | % known GWAS hits by ancestry | Barplot |
| **Fig 6** | Computational efficiency comparison | Scatter/Bar |

### Supplementary Tables

| Table | Content |
|-------|---------|
| **S1** | Variant inclusion/exclusion summary per array |
| **S2** | Per-chromosome variant counts and R² distributions |
| **S3** | Per-ancestry imputation quality metrics |
| **S4** | Correlation between global ancestry and imputation accuracy |
| **S5** | Standard INFO vs MagicalRsq-X filtering comparison |
| **S6** | Computational resource summary by approach |
| **S7** | Cross-platform reproducibility metrics |
| **S8** | Functional annotation enrichment by ancestry |

## Directory Structure

```
benchmarking/
├── README.md                      # This file
├── download_benchmark_data.sh     # Download public test data + WGS truth
├── run_all_benchmarks.sh          # Master script to run all approaches
│
├── alternative_approaches/        # Shell script pipelines (Traditional)
│   ├── approach_a_michigan.sh     # A: QC before + Intersect + Michigan/1KG
│   ├── approach_a_topmed.sh       # A: QC before + Intersect + TOPMed
│   ├── approach_a_allofus.sh      # A: QC before + Intersect + All of Us
│   ├── approach_b_topmed.sh       # B: QC after + Intersect + TOPMed
│   ├── approach_c_topmed.sh       # C: QC before + Separate platforms + TOPMed
│   ├── approach_d_topmed.sh       # D: QC after + Separate platforms + TOPMed
│   └── common_functions.sh        # Shared QC/processing functions
│
├── our_pipeline_variants/         # Our pipeline test configurations
│   ├── ours_1step_topmed.config   # E: Union + 1 imputation + MagicalRsq-X
│   ├── ours_1step_allofus.config  # E: Union + 1 imputation + MagicalRsq-X
│   ├── ours_2step_topmed.config   # F: Union + 2-step imputation + MagicalRsq-X
│   └── ours_2step_allofus.config  # F: Union + 2-step imputation + MagicalRsq-X
│
├── bench-helper-scripts/          # Comparison and metrics scripts
│   ├── calculate_concordance.R    # Calculate concordance with WGS
│   ├── stratify_by_ancestry.R     # Stratify metrics by GRAF-anc groups
│   ├── compare_approaches.R       # Compare all approaches
│   ├── generate_benchmark_report.R # Generate final report
│   ├── plot_metrics.R             # Visualization scripts
│   └── extract_timing_metrics.sh  # Extract computational metrics
│
├── modules/                       # Nextflow modules for benchmarking
│   ├── Benchmark_Concordance.nf   # Calculate concordance metrics
│   ├── Benchmark_Compare.nf       # Compare approaches
│   └── Benchmark_Report.nf        # Generate reports
│
├── test_data/                     # Downloaded test data
│   ├── genotypes/                 # Array genotype data
│   ├── wgs_truth/                 # WGS truth data
│   └── sample_info/               # Sample metadata + ancestry
│
└── results/                       # Benchmark results
    ├── approach_a_michigan/
    ├── approach_a_topmed/
    ├── approach_a_allofus/
    ├── ...
    ├── our_pipeline_topmed/
    ├── our_pipeline_allofus/
    └── comparison_report/
```

## Quick Start

```bash
# 1. Download benchmark data (genotypes + WGS truth)
./benchmarking/download_benchmark_data.sh -o ./benchmarking/test_data

# 2. Run GRAF-anc to classify samples by ancestry
nextflow run modules/Module7_Ancestry.nf \
    --input_plink benchmarking/test_data/genotypes/1kg_omni \
    --run_grafanc_only true \
    --outdir benchmarking/test_data/sample_info

# 3. Run all benchmark approaches
./benchmarking/run_all_benchmarks.sh

# 4. Generate comparison report
Rscript benchmarking/bench-helper-scripts/generate_benchmark_report.R
```

## Data Requirements

| Data Type | Size | Source |
|-----------|------|--------|
| 1KG Omni genotypes | ~5GB | EBI/NCBI |
| 1KG 30x WGS (chr subset) | ~50GB | IGSR |
| HGDP-1KG WGS (if available) | ~100GB | NYGC |

## Expected Outputs

1. **Concordance tables** per approach × server × ancestry
2. **Timing comparison** across approaches
3. **Variant yield comparison**
4. **Publication-ready figures**:
   - Concordance by MAF × ancestry (heatmap)
   - Approach comparison (bar plots)
   - Computational efficiency (scatter)

## GWAS Validation (Optional)

For GWAS replication to demonstrate improvement:
1. Use UK Biobank summary statistics (public)
2. Compare effect size concordance
3. Focus on:
   - Rare variant associations
   - Trans-ancestry meta-analysis traits
   - Admixed population-specific signals

## References

- TOPMed Imputation Server: https://imputation.biodatacatalyst.nhlbi.nih.gov/
- Michigan Imputation Server: https://imputationserver.sph.umich.edu/
- All of Us AnVIL: https://anvil.terra.bio/
- 1000 Genomes: https://www.internationalgenome.org/
- HGDP-1KG: https://gnomad.broadinstitute.org/news/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/
