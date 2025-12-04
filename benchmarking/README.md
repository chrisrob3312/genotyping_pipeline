# Benchmarking Imputation Approaches

## Literature Background

Our benchmarking framework is based on key findings from these studies:

### Key References

1. **Southam et al. 2011** (Eur J Hum Genet)
   - "The effect of genome-wide association scan quality control on imputation outcome"
   - Found that stringent pre-imputation QC can REDUCE imputation quality
   - Recommends minimal QC before imputation, thorough QC after
   - DOI: [10.1038/ejhg.2010.242](https://doi.org/10.1038/ejhg.2010.242)

2. **Charon et al. 2021** (Scientific Reports)
   - "Impact of pre- and post-variant filtration strategies on imputation"
   - Systematic comparison of QC timing strategies
   - Confirms post-imputation filtering often superior
   - DOI: [10.1038/s41598-021-85333-z](https://doi.org/10.1038/s41598-021-85333-z)

3. **Preprint 2024** (bioRxiv)
   - Updates on imputation strategies for diverse populations
   - DOI: [10.1101/2024.04.19.24306081](https://doi.org/10.1101/2024.04.19.24306081)

4. **Verma et al. 2014** (Frontiers in Genetics)
   - "Imputation and quality control steps for combining multiple genome-wide datasets"
   - Early workflow for multi-platform merging
   - **Note:** Pre-dates diversity-aware tools (MagicalRsq-X, diverse panels)
   - DOI: [10.3389/fgene.2014.00370](https://doi.org/10.3389/fgene.2014.00370)

---

## Benchmark Matrix

### 6 Approaches × 3 Servers = 18 Combinations

```
              TOPMed    AllOfUs   Michigan_1KG
Approach A      ✓         ✓           ✓
Approach B      ✓         ✓           ✓
Approach C      ✓         ✓           ✓
Approach D      ✓         ✓           ✓
Approach E      ✓         ✓           ✓
Approach F      ✓         ✓           ✓
```

---

## Approach Definitions

### Traditional Approaches (A-D)

| Approach | Name | QC Timing | Merge Strategy | Based On |
|----------|------|-----------|----------------|----------|
| **A** | Stringent QC-Before | Before imputation | Per-platform, then intersect | Traditional practice |
| **B** | Minimal QC-Before | Before imputation | Per-platform, then intersect | Southam 2011 |
| **C** | Intersect-First | Before imputation | Intersect first, then QC | Common practice |
| **D** | Intersect → QC-After | After merge/imputation | Intersect, merge, then QC after | Charon 2021 / Verma 2014 |

### Our Pipeline Approaches (E-F)

| Approach | Name | Strategy | Key Features |
|----------|------|----------|--------------|
| **E** | Ours 1-Step | Union → Impute → MagicalRsq-X → QC | Single imputation pass |
| **F** | Ours 2-Step | Union → Impute → Intersect → Re-impute → MagicalRsq-X → QC | Two imputation passes |

---

## Imputation Servers

| Server | Panel | Samples | Best For | Automation |
|--------|-------|---------|----------|------------|
| **TOPMed** | TOPMed r2 | ~97K | Diverse/admixed | imputationbot |
| **All of Us** | TOPMed-based | ~245K | Latino, AFR-AMR | terralab |
| **Michigan 1KG** | 1000G Phase 3 | ~2.5K | Quick baseline | imputationbot |

**Note:** Michigan HRC excluded as EUR-focused panel inappropriate for diverse cohorts.

---

## Detailed Approach Workflows

### Approach A: Stringent QC-Before (Traditional)

**Philosophy:** Clean data thoroughly before imputation
**Critique:** May remove valid variants in admixed populations

```
Per Platform:
├── MAF > 1% filter (removes rare variants)
├── HWE p < 1e-6 filter (may remove admixed structure)
├── Variant call rate > 98%
├── Sample call rate > 98%
├── Heterozygosity ± 3 SD
└── Relatedness filter (pi-hat > 0.25)
        ↓
[Intersect across platforms]
        ↓
[Submit to Imputation Server]
        ↓
[Simple R² > 0.3 filter]
        ↓
Final Output
```

### Approach B: Minimal QC-Before (Southam 2011 Recommendation)

**Philosophy:** Let imputation handle variant filtering; QC after
**Based on:** Southam finding that stringent pre-QC hurts imputation quality

```
Per Platform:
├── Variant call rate > 90% only
├── Sample call rate > 90% only
└── Extreme outliers only
        ↓
[Intersect across platforms]
        ↓
[Submit to Imputation Server]
        ↓
[Thorough Post-Imputation QC]
├── R² > 0.3 filter
├── MAF filter (if desired)
├── HWE filter
├── Call rate filters
└── Relatedness
        ↓
Final Output
```

### Approach C: Intersect-First

**Philosophy:** Ensure all platforms have identical variants before any processing
**Critique:** Loses platform-specific variants that could be informative

```
Platform 1 ─┐
Platform 2 ─┼─→ [Find Common Variants First]
Platform 3 ─┘           ↓
                [Apply QC to common variants]
                ├── Standard QC filters
                └── Per-platform sample QC
                        ↓
                [Merge all platforms]
                        ↓
                [Submit to Imputation Server]
                        ↓
                [Post-QC]
                        ↓
                Final Output
```

### Approach D: Intersect → Merge → QC-After (Verma/Charon style)

**Philosophy:** Intersect for consistency, but delay QC until after imputation
**Based on:** Verma 2014 workflow + Charon 2021 QC timing findings

```
Platform 1 ─┐
Platform 2 ─┼─→ [Find Common Variants]
Platform 3 ─┘           ↓
                [Minimal per-platform QC]
                ├── Call rate only
                └── Extreme outliers
                        ↓
                [Merge all samples]
                        ↓
                [Submit to Imputation Server]
                        ↓
                [THOROUGH Post-Imputation QC] ← Key from Charon 2021
                ├── R² or MagicalRsq-X filter
                ├── MAF filter
                ├── HWE filter (cautious in admixed)
                ├── Call rate filters
                └── Relatedness (GENESIS preferred)
                        ↓
                Final Output
```

### Approach E: Our Pipeline (1-Step)

**Philosophy:** Preserve all variants, use ancestry-aware tools, QC after imputation
**Key innovations:** Union merge, MagicalRsq-X, skip HWE for admixed

```
Platform 1 ─┐
Platform 2 ─┼─→ [UNION Merge - Keep ALL Variants]
Platform 3 ─┘           ↓
                [Minimal QC]
                ├── Basic call rate
                └── Strand check
                        ↓
                [Submit to Imputation Server]
                        ↓
                [MagicalRsq-X Filter] ← Ancestry-calibrated
                        ↓
                [Thorough Post-QC]
                ├── Call rate filters
                ├── Skip HWE for admixed populations
                └── GENESIS PCRelate (handles admixture)
                        ↓
                Final Output
```

### Approach F: Our Pipeline (2-Step)

**Philosophy:** Same as E, plus re-imputation to recover variants lost in merge
**Rationale:** Second imputation pass on merged data fills gaps

```
Platform 1 ─┐
Platform 2 ─┼─→ [UNION Merge]
Platform 3 ─┘       ↓
              [Minimal QC]
                    ↓
              [Submit to Server - Pass 1]
                    ↓
              [Intersect Merge across platforms]
                    ↓
              [Submit to Server - Pass 2] ← Re-imputation
                    ↓
              [MagicalRsq-X Filter]
                    ↓
              [Thorough Post-QC]
                    ↓
              Final Output
```

---

## Why Traditional Approaches May Fail for Diverse Populations

| Issue | Traditional Impact | Our Solution |
|-------|-------------------|--------------|
| **HWE filtering before imputation** | Removes real variants in admixed populations | Skip HWE or apply after, cautiously |
| **MAF filtering before imputation** | Loses rare population-specific variants | No MAF filter before imputation |
| **Static R² threshold** | Penalizes non-EUR ancestry | MagicalRsq-X ancestry-calibrated |
| **Intersect merge** | Loses platform-unique variants | Union merge preserves all |
| **PLINK relatedness** | Assumes homogeneous population | GENESIS PCRelate handles admixture |
| **EUR-focused panels** | Poor imputation for AFR/AMR | TOPMed/All of Us diverse panels |

---

## Expected Outcomes

Based on literature and panel composition:

| Ancestry | Best Traditional | Our Pipeline Advantage |
|----------|-----------------|----------------------|
| **EUR** | B or D | Slight improvement |
| **AFR** | B or D | Significant improvement |
| **AMR** | B or D | Major improvement (Latino/NAT) |
| **Admixed** | D | Major improvement |
| **Rare variants** | Poor all | Better with TOPMed |

---

## Running the Full Matrix

```bash
# Setup credentials (see SETUP_IMPUTATIONBOT.md)
export TOPMED_TOKEN="your_token"
export MICHIGAN_TOKEN="your_token"
terralab login  # For All of Us

# Run full matrix (18 combinations)
./run_full_benchmark_matrix.sh \
    --input benchmarking/test_data \
    --output benchmark_results \
    --servers "topmed,allofus,michigan_1kg"

# Or run specific approaches
./run_full_benchmark_matrix.sh \
    --input test_data \
    --approaches "A,B,E,F" \
    --servers "topmed,allofus"
```

---

## Directory Structure

```
benchmarking/
├── README.md                          # This file
├── SETUP_IMPUTATIONBOT.md            # Credential setup guide
├── GWAS_SIMULATION_GUIDE.md          # Phenotype simulation
│
├── download_benchmark_data.sh         # Download 1KG test data
├── download_gwas_sumstats.sh          # Download GWAS data
├── simulate_gwas_phenotypes.sh        # Simulate phenotypes
├── run_all_benchmarks.sh              # Simple benchmark runner
├── run_full_benchmark_matrix.sh       # Full 18-combination matrix
│
├── alternative_approaches/            # Traditional approach scripts
│   ├── approach_a_topmed.sh          # A × TOPMed
│   ├── approach_a_allofus.sh         # A × All of Us
│   ├── approach_a_michigan_1kg.sh    # A × Michigan 1KG
│   ├── approach_b_topmed.sh          # B × TOPMed
│   ├── approach_b_allofus.sh         # B × All of Us
│   ├── approach_b_michigan_1kg.sh    # B × Michigan 1KG
│   ├── approach_c_*.sh               # C variants
│   ├── approach_d_*.sh               # D variants
│   └── common_functions.sh           # Shared functions
│
├── our_pipeline_variants/             # Our pipeline configs
│   ├── ours_1step_topmed.config      # E × TOPMed
│   ├── ours_1step_allofus.config     # E × All of Us
│   ├── ours_1step_michigan_1kg.config # E × Michigan 1KG
│   ├── ours_2step_topmed.config      # F × TOPMed
│   ├── ours_2step_allofus.config     # F × All of Us
│   └── ours_2step_michigan_1kg.config # F × Michigan 1KG
│
├── bench-helper-scripts/              # Analysis scripts
│   ├── calculate_concordance.R
│   ├── compare_approaches.R
│   ├── compare_rvas_approaches.R
│   ├── benchmark_rare_variants.R
│   ├── generate_publication_figures.R
│   ├── summarize_benchmark_timings.R
│   ├── simulate_phenotypes.R
│   ├── simulate_ancestry_aware_phenotypes.R
│   ├── simulate_lai_phenotypes.R
│   ├── check_hit_recovery.R
│   ├── calculate_lai_metrics.R
│   └── extract_timing_metrics.sh
│
└── test_data/                         # Downloaded data
    ├── genotypes/
    ├── wgs_truth/
    ├── phenotypes/
    └── sample_info/
```

---

## Metrics Calculated

### Imputation Quality
- INFO/R² by MAF bin and ancestry
- MagicalRsq-X scores
- Variant retention rates

### Concordance with WGS Truth
- Genotype concordance (exact match)
- Dosage R² (correlation)
- Non-reference concordance
- By MAF × ancestry stratification

### GWAS/PRS Performance
- GWAS hit recovery rate
- Effect size correlation with published
- PRS R² by ancestry
- Power comparison

### Computational
- Wall-clock time per step
- Total pipeline time
- Variants per minute efficiency

---

## Citations

If using this framework, please cite:

1. Southam L et al. (2011) Eur J Hum Genet. DOI: 10.1038/ejhg.2010.242
2. Charon C et al. (2021) Sci Rep. DOI: 10.1038/s41598-021-85333-z
3. Verma SS et al. (2014) Front Genet. DOI: 10.3389/fgene.2014.00370
4. Sun Q et al. (2024) AJHG - MagicalRsq-X
5. TOPMed Imputation Reference Panel
6. [Your pipeline citation]
