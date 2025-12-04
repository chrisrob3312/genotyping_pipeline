# Benchmarking Framework

Comprehensive benchmarking to compare imputation strategies across ancestry groups for mixed-ancestry cohorts.

## Overview

This framework compares **6 imputation approaches** across **3 imputation servers** using publicly available data with WGS truth for validation.

## Approaches Compared

| Approach | QC Timing | Merge Strategy | Description |
|----------|-----------|----------------|-------------|
| **A** | Before imputation | Intersect before | Traditional: QC → Intersect → Impute → R² filter |
| **B** | After imputation | Intersect after | Minimal prep → Impute → QC after |
| **C** | Before imputation | Separate platforms | Thorough QC → Impute each platform → Merge |
| **D** | After imputation | Separate platforms | Minimal prep → Impute each → QC after → Merge |
| **E (Ours)** | Minimal before, thorough after | Union merge | Our pipeline with MagicalRsq-X 2-step filtering |
| **F (Ours+)** | Same as E | Union + re-impute | Our pipeline with optional re-imputation step |

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

### Per-Ancestry Group (GRAF-anc 8 Major Groups)

| Metric | Description | Expected (Common) | Expected (Rare) |
|--------|-------------|-------------------|-----------------|
| Dosage R² | Correlation with WGS | >0.95 | >0.70 |
| Genotype Concordance | Exact GT match | >98% | >90% |
| Non-Reference Concordance | Alt allele accuracy | >95% | >85% |
| Variant Yield | Variants passing QC | - | - |
| Imputation Rate | % variants imputed | - | - |

### Stratified By
- MAF bins: 0-0.5%, 0.5-1%, 1-5%, 5-50%
- Ancestry group: AFR, EUR, EAS, SAS, AMR, MEN, OCN, MIX
- Chromosome (to check consistency)

### Computational Metrics
- Wall-clock time
- CPU hours
- Peak memory usage
- Storage requirements

## Directory Structure

```
benchmarking/
├── README.md                      # This file
├── download_benchmark_data.sh     # Download public test data + WGS truth
├── run_all_benchmarks.sh          # Master script to run all approaches
│
├── alternative_approaches/        # Shell script pipelines (non-Nextflow)
│   ├── approach_a_michigan.sh     # Approach A + Michigan
│   ├── approach_a_topmed.sh       # Approach A + TOPMed
│   ├── approach_a_allofus.sh      # Approach A + All of Us
│   ├── approach_b_michigan.sh     # Approach B + Michigan
│   ├── approach_b_topmed.sh       # Approach B + TOPMed
│   ├── approach_b_allofus.sh      # Approach B + All of Us
│   ├── approach_c_michigan.sh     # Approach C + Michigan
│   ├── approach_c_topmed.sh       # Approach C + TOPMed
│   ├── approach_c_allofus.sh      # Approach C + All of Us
│   ├── approach_d_michigan.sh     # Approach D + Michigan
│   ├── approach_d_topmed.sh       # Approach D + TOPMed
│   └── approach_d_allofus.sh      # Approach D + All of Us
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
