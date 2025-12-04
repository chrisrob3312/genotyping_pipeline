# COMPLETE PIPELINE WORKFLOW DOCUMENTATION

**Last Updated:** December 2025
**Version:** 1.0

---

## TABLE OF CONTENTS

1. [Pipeline Overview](#pipeline-overview)
2. [Module 0: Container Build](#module-0-container-build)
3. [Module 1: Pre-Imputation QC](#module-1-pre-imputation-qc)
4. [Module 2: Imputation](#module-2-imputation)
5. [Module 3: Post-Imputation QC](#module-3-post-imputation-qc)
6. [Module 4: Platform Merging](#module-4-platform-merging)
7. [Module 5: Re-Imputation](#module-5-re-imputation)
8. [Module 6: Post-Merge QC](#module-6-post-merge-qc)
9. [Module 7: Ancestry Estimation](#module-7-ancestry-estimation)
10. [Local Ancestry Inference Tools](#local-ancestry-inference-tools)

---

## PIPELINE OVERVIEW

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                          COMPLETE PIPELINE FLOW                              │
└─────────────────────────────────────────────────────────────────────────────┘

                              ┌──────────────────┐
                              │   Sample Sheet   │
                              │      (CSV)       │
                              └────────┬─────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│ MODULE 1: PRE-IMPUTATION QC                                                   │
│ ─────────────────────────────────────────────────────────────────────────── │
│ • Discover input files (4 structure types)                                   │
│ • Merge samples → batches → platforms (UNION)                               │
│ • Align to reference (bcftools +fixref)                                     │
│ • Liftover hg19 → hg38 (CrossMap)                                           │
│ • Light QC (call rates, duplicates, monomorphic)                            │
│ • Create service-specific VCFs                                              │
└──────────────────────────────────────────────────────────────────────────────┘
                          │                    │
                   TOPMed VCFs           AnVIL VCF
                   (22 per platform)     (1 per platform)
                          │                    │
                          ▼                    ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│ MODULE 2: IMPUTATION                                                          │
│ ─────────────────────────────────────────────────────────────────────────── │
│ TOPMed Server:                     │  All of Us AnVIL:                       │
│ • imputationbot submit             │  • terralab login (OAuth)               │
│ • Auto-monitor & download          │  • terralab submit array_imputation     │
│ • --refpanel topmed-r3             │  • Auto-monitor & download              │
│ • --autoDownload                   │  • --multiSampleVcf                     │
└──────────────────────────────────────────────────────────────────────────────┘
                                       │
                          ┌────────────┴────────────┐
                          ▼                         ▼
                   TOPMed Results            AnVIL Results
                          │                         │
                          └────────────┬────────────┘
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│ MODULE 3: POST-IMPUTATION QC                                                  │
│ ─────────────────────────────────────────────────────────────────────────── │
│ • MagicalRsq-X imputation quality filtering (R² ≥ 0.3)                       │
│ • Sample call rate filter (≥ 95%)                                            │
│ • Variant call rate filter (≥ 95%)                                           │
│ • Per-platform QC reports                                                    │
└──────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│ MODULE 4: PLATFORM MERGING                                                    │
│ ─────────────────────────────────────────────────────────────────────────── │
│ • Cross-platform variant harmonization                                       │
│ • INTERSECTION merge (shared variants only)                                  │
│ • Multi-allelic resolution                                                   │
│ • Strand alignment verification                                              │
└──────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│ MODULE 5: RE-IMPUTATION (OPTIONAL)                                            │
│ ─────────────────────────────────────────────────────────────────────────── │
│ • Same services as Module 2                                                  │
│ • Re-impute merged dataset                                                   │
│ • Fill gaps from intersection merging                                        │
│ • Optional: Can skip if not needed                                           │
└──────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│ MODULE 6: POST-MERGE QC (FINAL)                                               │
│ ─────────────────────────────────────────────────────────────────────────── │
│ • MagicalRsq-X second-pass filtering                                         │
│ • Hardy-Weinberg Equilibrium filter (p < 1e-6)                               │
│ • MAF filtering (separate MAF/noMAF outputs)                                 │
│ • Heterozygosity outlier detection (3 SD)                                    │
│ • Relatedness detection (GENESIS PCRelate or PLINK)                          │
│ • PCA for population structure                                               │
└──────────────────────────────────────────────────────────────────────────────┘
                                       │
                          ┌────────────┴────────────┐
                          ▼                         ▼
                   MAF-filtered             No MAF filter
                   (for GWAS)               (for ancestry)
                          │                         │
                          └────────────┬────────────┘
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│ MODULE 7: ANCESTRY ESTIMATION                                                 │
│ ─────────────────────────────────────────────────────────────────────────── │
│ GLOBAL ANCESTRY:                                                             │
│ • GRAF-anc (8 major groups + 50+ subgroups)                                  │
│ • ADMIXTURE (K=2-12, cross-validation)                                       │
│                                                                               │
│ LOCAL ANCESTRY (OPTIONAL):                                                   │
│ • RFMix v2 (default) - Random forest + CRF                                   │
│ • RFMix v1           - PopPhased variant                                     │
│ • FLARE              - Fast HMM-based                                        │
│ • G-NOMIX            - Neural network-based                                  │
│                                                                               │
│ VISUALIZATIONS:                                                              │
│ • PCA plots colored by ancestry                                              │
│ • ADMIXTURE bar plots (sorted by ancestry)                                   │
│ • Summary-level LAI statistics (cohort-level)                                │
└──────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
                           ┌───────────────────────┐
                           │   FINAL OUTPUTS       │
                           ├───────────────────────┤
                           │ • Analysis-ready PLINK│
                           │ • Ancestry estimates  │
                           │ • QC reports          │
                           │ • Visualizations      │
                           └───────────────────────┘
```

---

## MODULE 0: CONTAINER BUILD

### Purpose
Builds all Apptainer/Singularity containers required by the pipeline.

### Containers Built

| Container | Contents | Used By |
|-----------|----------|---------|
| `plink_1.9.sif` | PLINK 1.9, PLINK 2.0 | Modules 1, 3, 4, 6 |
| `bcftools.sif` | bcftools 1.18+, tabix | Modules 1, 2, 3 |
| `perl_crossmap.sif` | Perl, CrossMap | Module 1 |
| `python_api.sif` | Python 3.11, imputationbot, terralab-cli | Module 2, 5 |
| `r_genetics.sif` | R 4.x, genetics packages | Modules 3, 6, 7 |
| `ancestry_suite.sif` | GRAF-anc, ADMIXTURE, RFMix, FLARE | Module 7 |

### Usage

```bash
nextflow run modules/Module0_Apptainer_Build.nf -profile apptainer
```

---

## MODULE 1: PRE-IMPUTATION QC

Module 1 is responsible for harmonizing diverse input data formats into standardized, imputation-ready VCFs. The pipeline is designed to handle heterogeneous data sources while preserving maximum genetic diversity through UNION merge strategies.

### Input Flexibility & Data Harmonization

The pipeline accepts 4 different input file structures, automatically detecting and handling each:

| File Structure | Description | Example |
|----------------|-------------|---------|
| `individual_samples` | Separate files per sample | `/path/sample001.bed`, `sample002.bed`, etc. |
| `individual_chr_split` | Per-sample, chromosome-split | `/path/sample001_chr*.bed` |
| `merged_batch` | Pre-merged batch file | `/path/batch_A.bed` (all samples combined) |
| `merged_chr_split` | Pre-merged, chromosome-split | `/path/batch_A_chr*.bed` |

**Format Support:** Both PLINK (bed/bim/fam) and VCF inputs are supported with automatic conversion.

**Build Detection:** Genome build (hg19/hg38) is detected per-batch, not globally, allowing mixed-build inputs.

### Summary Flow

```
INPUT                    PROCESS                           OUTPUT
─────                    ───────                           ──────
Sample Sheet  ───────►  Discover Files (4 structure types)
                              │
                              ▼
                        Prepare Samples (VCF→PLINK, concat chr-split)
                              │
                              ▼
                        Merge to Batch (UNION - preserves all variants)
                              │
                              ▼
                        Align to Reference (bcftools +fixref)
                              │
                              ▼
              ┌───────────────┴───────────────┐
              ▼                               ▼
        hg19 batches                    hg38 batches
              │                               │
              ▼                               │
        Liftover (CrossMap VCF mode)          │
              │                               │
              └───────────────┬───────────────┘
                              ▼
                        All hg38
                              │
                              ▼
                        Merge to Platform (UNION - ~99.999% retention)
                              │
                              ▼
              ┌───────────────┴───────────────┐
              ▼                               ▼
        TOPMed Validation              AnVIL Validation
        (Will Rayner script)           (passthrough)
              │                               │
              ▼                               ▼
        Light QC                        Light QC
        (biallelic, dups, call rate)   (biallelic, dups, call rate)
              │                               │
              ▼                               ▼
        Create TOPMed VCFs        Create AnVIL VCF  ───► Single VCF
        (22 per platform)                               (all autosomes)
```

### Reference Alignment Process

1. **Convert to VCF**: `plink2 --export vcf-4.2 bgz`
2. **Split multiallelics**: `bcftools norm -m-any`
3. **Check REF**: `bcftools norm --check-ref w`
4. **Fix mismatches**: `bcftools +fixref` (strand flips, REF/ALT swaps)
5. **Sort and index**: `bcftools sort && bcftools index -c`

### Liftover Process (hg19 batches only)

Uses CrossMap VCF mode (not BED mode) for accurate coordinate conversion:
- Updates genomic coordinates hg19 → hg38
- Validates REF alleles against hg38 reference
- Handles strand orientation automatically
- Creates `.unmap` file for failed lifts

### Light QC Applied

| Filter | Command | Purpose |
|--------|---------|---------|
| Biallelic SNPs | `--snps-only just-acgt --max-alleles 2` | Remove indels, multiallelic |
| Remove duplicates | `--rm-dup exclude-all` | Eliminate duplicate positions |
| Remove monomorphic | `--maf 0.000001` | Only removes truly invariant sites |
| Variant call rate | `--geno 0.05` | ≥95% call rate per variant |
| Sample call rate | `--mind 0.05` | ≥95% call rate per sample |

**Explicitly NOT applied in Module 1:**
- HWE filtering (Module 6 - optional)
- MAF filtering (Module 6)
- Heterozygosity filtering (Module 6)
- Relatedness filtering (Module 6)

### Service-Specific Outputs

| Service | Output Format | Files | Notes |
|---------|--------------|-------|-------|
| TOPMed | VCF 4.2, CSI indexed | 22 per platform (chr1-22) | All chromosomes submitted as single job |
| All of Us AnVIL | VCF 4.3, CSI indexed | 1 per platform (all autosomes) | 95% quota savings vs per-chr |

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sample_call_rate` | 0.95 | Minimum sample call rate |
| `variant_call_rate` | 0.95 | Minimum variant call rate |

See `documentation/Module_1_Workflow.md` for complete process-level details.

---

## MODULE 2: IMPUTATION

### Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           MODULE 2: IMPUTATION                               │
└─────────────────────────────────────────────────────────────────────────────┘

FROM MODULE 1:
  ├── TOPMed VCFs (22 per platform)
  └── AnVIL VCF (1 per platform)
              │
              ▼
┌═════════════════════════════════════════════════════════════════════════════┐
║ PROCESS: groupChromosomesByPlatform                                          ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ Group all chromosome VCFs per platform for single submission                 ║
║ Output: Manifest file listing all chr1-22 VCFs                              ║
└═════════════════════════════════════════════════════════════════════════════┘
              │
              ├─────────────────────────────────────────┐
              ▼                                         ▼
┌═══════════════════════════════════════┐ ┌═══════════════════════════════════┐
║ PROCESS: submitToTOPMed               ║ ║ PROCESS: submitToAllOfUs          ║
╠───────────────────────────────────────╢ ╠───────────────────────────────────╢
║ Container: python_api.sif             ║ ║ Container: python_api.sif         ║
║                                       ║ ║                                   ║
║ Commands:                             ║ ║ Prerequisites:                    ║
║  imputationbot impute \               ║ ║  - terralab login (OAuth)         ║
║    --files chr1.vcf.gz ... chr22.vcf ║ ║  - Run on login node first        ║
║    --refpanel topmed-r3 \             ║ ║                                   ║
║    --autoDownload \                   ║ ║ Commands:                         ║
║    --password $DECRYPT_PASSWORD \     ║ ║  terralab submit array_imputation ║
║    --r2Filter 0                       ║ ║    --multiSampleVcf input.vcf.gz  ║
║                                       ║ ║    --outputBasename platform_anvil║
║ Limits:                               ║ ║    --description "..."            ║
║  - maxForks 3 (TOPMed limit)         ║ ║                                   ║
║  - Auto-monitors job status           ║ ║ Monitoring:                       ║
║  - Auto-downloads on completion       ║ ║  terralab jobs details --id $ID   ║
║                                       ║ ║  terralab download --id $ID       ║
╠───────────────────────────────────────╢ ╠───────────────────────────────────╢
║ Output:                               ║ ║ Output:                           ║
║  - Imputed VCFs (chr1-22)            ║ ║  - Imputed VCFs                   ║
║  - INFO scores                        ║ ║  - INFO scores                    ║
║  - QC metrics                         ║ ║  - QC metrics                     ║
└═══════════════════════════════════════┘ └═══════════════════════════════════┘
              │                                         │
              └─────────────────┬───────────────────────┘
                                ▼
                         TO MODULE 3
```

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `run_topmed` | true | Enable TOPMed imputation |
| `run_anvil` | false | Enable All of Us imputation |
| `topmed_api_token` | null | TOPMed API token |
| `topmed_password` | null | Password for auto-decrypt |
| `monitor_interval_minutes` | 45 | Check interval for job status |

---

## MODULE 3: POST-IMPUTATION QC

### Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                      MODULE 3: POST-IMPUTATION QC                            │
└─────────────────────────────────────────────────────────────────────────────┘

FROM MODULE 2:
  └── Imputed VCFs + INFO files
              │
              ▼
┌═════════════════════════════════════════════════════════════════════════════┐
║ PROCESS: INDEX_VCF                                                           ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ bcftools index -c ${vcf}                                                     ║
║ Create CSI index for each imputed VCF                                        ║
└═════════════════════════════════════════════════════════════════════════════┘
              │
              ▼
┌═════════════════════════════════════════════════════════════════════════════┐
║ PROCESS: MAGICALRSQ_FILTER                                                   ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ Container: r_genetics.sif                                                    ║
║                                                                              ║
║ Script: magicalrsq_filter.R                                                  ║
║                                                                              ║
║ Steps:                                                                       ║
║   1. Parse INFO scores from imputation                                       ║
║   2. Apply MagicalRsq-X model (population-specific)                         ║
║   3. Filter variants with R² < threshold (default 0.3)                      ║
║   4. Generate QC report                                                      ║
║                                                                              ║
║ Parameters:                                                                  ║
║   --threshold ${params.magicalrsq_threshold}                                ║
║   --models-dir ${params.magicalrsq_models_dir}                              ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ Output: Filtered VCFs + QC reports                                           ║
└═════════════════════════════════════════════════════════════════════════════┘
              │
              ▼
┌═════════════════════════════════════════════════════════════════════════════┐
║ PROCESS: SAMPLE_VARIANT_QC                                                   ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ plink2 --geno 0.05 --mind 0.05 ...                                          ║
║                                                                              ║
║ Filters:                                                                     ║
║   - Variant call rate ≥ 95%                                                 ║
║   - Sample call rate ≥ 95%                                                  ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ Output: QC'd VCFs ready for Module 4                                         ║
└═════════════════════════════════════════════════════════════════════════════┘
              │
              ▼
                         TO MODULE 4
```

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `magicalrsq_threshold` | 0.3 | MagicalRsq-X R² threshold |
| `sample_call_rate` | 0.95 | Minimum sample call rate |
| `variant_call_rate` | 0.95 | Minimum variant call rate |

---

## MODULE 4: PLATFORM MERGING

### Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                       MODULE 4: PLATFORM MERGING                             │
└─────────────────────────────────────────────────────────────────────────────┘

FROM MODULE 3:
  ├── Platform A VCFs (QC'd)
  ├── Platform B VCFs (QC'd)
  └── Platform N VCFs (QC'd)
              │
              ▼
┌═════════════════════════════════════════════════════════════════════════════┐
║ PROCESS: VCF_TO_PLINK                                                        ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ plink2 --vcf ${vcf} dosage=DS --make-pgen                                   ║
║ Convert VCFs to PLINK2 format (pgen/pvar/psam)                              ║
║ Preserve dosage information                                                  ║
└═════════════════════════════════════════════════════════════════════════════┘
              │
              ▼
┌═════════════════════════════════════════════════════════════════════════════┐
║ PROCESS: FIRST_PASS_MERGE                                                    ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ plink2 --pmerge-list platforms.txt --merge-mode 1 (INTERSECTION)            ║
║                                                                              ║
║ NOTE: INTERSECTION merge - keeps only variants present in ALL platforms     ║
║ (Different from Module 1 which uses UNION)                                   ║
║                                                                              ║
║ Handles:                                                                     ║
║   - Variant ID harmonization                                                 ║
║   - Strand conflicts detection                                               ║
║   - Multi-allelic resolution                                                 ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ May fail with: "variant position collisions" or "strand issues"             ║
└═════════════════════════════════════════════════════════════════════════════┘
              │
              ├── Success ────────────────────────────────┐
              │                                           │
              ▼ (if conflicts)                            │
┌═════════════════════════════════════════════════════════════════════════════┐
║ PROCESS: FIX_MERGE_CONFLICTS                                                 ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ Parse .missnp file from failed merge                                         ║
║ Flip strands or exclude unresolvable conflicts                              ║
║ Re-attempt merge                                                             ║
└═════════════════════════════════════════════════════════════════════════════┘
              │                                           │
              └───────────────────────────────────────────┤
                                                          ▼
┌═════════════════════════════════════════════════════════════════════════════┐
║ PROCESS: PLINK_TO_VCF                                                        ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ plink2 --export vcf-4.2 bgz                                                  ║
║ Convert merged PLINK to VCF format                                           ║
║ bcftools index -c                                                            ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ Output: Merged VCF (all platforms combined)                                  ║
└═════════════════════════════════════════════════════════════════════════════┘
              │
              ▼
                    TO MODULE 5 (or MODULE 6 if skipping)
```

---

## MODULE 5: RE-IMPUTATION

### Purpose
Re-impute the merged dataset to fill gaps created by intersection merging.

### Flow
Uses same processes as Module 2, but with:
- Single merged VCF input (instead of per-platform)
- May use different reference panel settings

### When to Skip
Set `--skip_modules "5"` if:
- Platforms have high variant overlap
- Imputation quota is limited
- Speed is priority over completeness

---

## MODULE 6: POST-MERGE QC

### Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                      MODULE 6: POST-MERGE QC (FINAL)                         │
└─────────────────────────────────────────────────────────────────────────────┘

FROM MODULE 5 (or MODULE 4):
  └── Merged/Re-imputed VCFs
              │
              ▼
┌═════════════════════════════════════════════════════════════════════════════┐
║ PROCESS: SECOND_MAGICALRSQ                                                   ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ Second-pass MagicalRsq-X filtering                                           ║
║ Re-evaluate imputation quality after merging                                 ║
└═════════════════════════════════════════════════════════════════════════════┘
              │
              ▼
┌═════════════════════════════════════════════════════════════════════════════┐
║ PROCESS: HWE_FILTER (OPTIONAL - skipped by default)                          ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ plink2 --hwe ${params.hwe_pvalue} midp                                       ║
║ Remove variants out of HWE (p < 1e-6)                                        ║
║ Using mid-p correction for small samples                                     ║
║                                                                              ║
║ SKIPPED BY DEFAULT (params.skip_hwe = true) for admixed populations         ║
║ Set --skip_hwe false to enable HWE filtering                                ║
└═════════════════════════════════════════════════════════════════════════════┘
              │
              ▼
┌═════════════════════════════════════════════════════════════════════════════┐
║ PROCESS: HETEROZYGOSITY_FILTER                                               ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ plink2 --het                                                                 ║
║ Remove samples with heterozygosity > mean ± 3 SD                            ║
║ Identifies sample contamination or quality issues                            ║
└═════════════════════════════════════════════════════════════════════════════┘
              │
              ▼
┌═════════════════════════════════════════════════════════════════════════════┐
║ PROCESS: RELATEDNESS_FILTER                                                  ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ IF params.use_genesis:                                                       ║
║   GENESIS PCRelate (GRM-based, handles admixture)                           ║
║ ELSE:                                                                        ║
║   plink2 --king-cutoff ${params.kinship_threshold}                          ║
║                                                                              ║
║ Identify related pairs (kinship > 0.177 = 1st degree)                       ║
║ Remove one from each pair (maximize sample retention)                        ║
└═════════════════════════════════════════════════════════════════════════════┘
              │
              ▼
┌═════════════════════════════════════════════════════════════════════════════┐
║ PROCESS: CALCULATE_PCS                                                       ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ plink2 --pca 20 approx                                                       ║
║ Calculate top 20 principal components                                        ║
║ Use approximate algorithm for large samples                                  ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ Output: .eigenvec (PCs) + .eigenval (variance explained)                    ║
└═════════════════════════════════════════════════════════════════════════════┘
              │
              ▼
┌═════════════════════════════════════════════════════════════════════════════┐
║ PROCESS: MAF_FILTER_SPLIT                                                    ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ Create TWO outputs:                                                          ║
║                                                                              ║
║ 1. MAF-filtered (for GWAS):                                                  ║
║    plink2 --maf 0.01 --make-bed --out final_maf_filtered                    ║
║                                                                              ║
║ 2. No MAF filter (for ancestry):                                            ║
║    plink2 --make-bed --out final_no_maf                                     ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ Output: Two PLINK filesets                                                   ║
└═════════════════════════════════════════════════════════════════════════════┘
              │
              ▼
                         TO MODULE 7
```

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `skip_hwe` | true | Skip HWE filtering (recommended for admixed populations) |
| `hwe_pvalue` | 1e-6 | HWE filter p-value threshold (if skip_hwe=false) |
| `het_sd_threshold` | 3 | Heterozygosity SD threshold |
| `use_genesis` | true | Use GENESIS PCRelate (vs PLINK) |
| `kinship_threshold` | 0.177 | 1st-degree relatedness cutoff |

### Imputation Server Comparison

When both TOPMed and All of Us imputation are run, Module 6 generates comparative metrics:

```
┌═════════════════════════════════════════════════════════════════════════════┐
║ PROCESS: COMPARE_IMPUTATION_PERFORMANCE                                      ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ Compares TOPMed vs All of Us AnVIL results                                   ║
║                                                                              ║
║ Metrics Calculated:                                                          ║
║   • Total variant counts per server                                          ║
║   • Shared variants (intersection)                                           ║
║   • Server-unique variants                                                   ║
║   • MAF distribution comparison                                              ║
║   • INFO/R² score distributions                                              ║
║                                                                              ║
║ Output:                                                                      ║
║   • HTML comparison report                                                   ║
║   • PDF visualization plots                                                  ║
║   • Tab-delimited metrics file                                               ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ Location: ${outdir}/module6/13_imputation_comparison/                        ║
└═════════════════════════════════════════════════════════════════════════════┘
```

---

## MODULE 7: ANCESTRY ESTIMATION

### Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                     MODULE 7: ANCESTRY ESTIMATION                            │
└─────────────────────────────────────────────────────────────────────────────┘

FROM MODULE 6:
  ├── Final PLINK (no MAF filter)
  └── PCA results
              │
              ▼
┌═════════════════════════════════════════════════════════════════════════════┐
║ PROCESS: PREPARE_FOR_ANCESTRY                                                ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ LD pruning for ancestry analysis:                                            ║
║   plink --indep-pairwise 200 50 0.25                                        ║
║   plink --extract prune.in --make-bed                                       ║
║                                                                              ║
║ Moderate pruning to retain ancestry-informative SNPs                        ║
└═════════════════════════════════════════════════════════════════════════════┘
              │
              ├─────────────────────────────────────────────────────┐
              ▼                                                     ▼
┌═════════════════════════════════════════┐ ┌═════════════════════════════════┐
║ PROCESS: GRAFANC_ANCESTRY               ║ ║ PROCESS: RUN_ADMIXTURE          ║
╠─────────────────────────────────────────╢ ╠─────────────────────────────────╢
║ Container: ancestry_suite.sif           ║ ║ Container: ancestry_suite.sif   ║
║                                         ║ ║                                 ║
║ Command:                                ║ ║ For K in {2..12}:               ║
║   grafanc ${prefix} \                   ║ ║   admixture ${bed} ${K} --cv=10 ║
║     ${output}.txt \                     ║ ║                                 ║
║     --threads ${cpus}                   ║ ║ Cross-validation to find best K ║
║                                         ║ ║ Plot CV error curve             ║
║ Outputs:                                ║ ║                                 ║
║ • 8 major ancestry groups:              ║ ║ Outputs:                        ║
║   AFR, MEN, EUR, SAS, EAS, AMR, OCN, MIX║ ║ • .Q files (ancestry fractions) ║
║ • 50+ subcontinental groups             ║ ║ • .P files (allele frequencies) ║
║ • GD1-GD6 genetic distances             ║ ║ • CV error log                  ║
║ • Probability estimates (Pe, Pf, Pa)    ║ ║ • Bar plots per K               ║
╠─────────────────────────────────────────╢ ╠─────────────────────────────────╢
║ Visualization:                          ║ ║ Visualization:                  ║
║ • GD scatter plots                      ║ ║ • Structure-style bar plots     ║
║ • Ancestry distribution histograms      ║ ║ • Sorted by dominant ancestry   ║
║ • Summary statistics by group           ║ ║ • Grouped by GRAF-anc groups    ║
└═════════════════════════════════════════┘ └═════════════════════════════════┘
              │                                         │
              └─────────────────┬───────────────────────┘
                                │
                                ▼
              ┌─────────────────────────────────────────┐
              │     IF params.run_lai == true           │
              └─────────────────┬───────────────────────┘
                                │
                                ▼
┌═════════════════════════════════════════════════════════════════════════════┐
║ LOCAL ANCESTRY INFERENCE (OPTIONAL)                                          ║
╠═════════════════════════════════════════════════════════════════════════════╣
║                                                                              ║
║ Choice based on params.lai_methods:                                          ║
║                                                                              ║
║ ┌─────────────────────────────────────────────────────────────────────────┐ ║
║ │ RFMIX v2 (default)                                                      │ ║
║ │ rfmix -f query.vcf -r reference.vcf -m sample_map -g genetic_map       │ ║
║ │       -o output --chromosome=${chr} -e 5 --reanalyze-reference         │ ║
║ │ Output: .msp.tsv, .rfmix.Q, .fb.tsv                                    │ ║
║ └─────────────────────────────────────────────────────────────────────────┘ ║
║                                                                              ║
║ ┌─────────────────────────────────────────────────────────────────────────┐ ║
║ │ RFMIX v1                                                                │ ║
║ │ RFMix_PopPhased -a alleles.txt -p classes.txt -m markers.txt -o output │ ║
║ │ Note: Requires different input format (binary alleles)                 │ ║
║ └─────────────────────────────────────────────────────────────────────────┘ ║
║                                                                              ║
║ ┌─────────────────────────────────────────────────────────────────────────┐ ║
║ │ FLARE                                                                   │ ║
║ │ java -Xmx${mem}g -jar flare.jar ref=ref.vcf ref-panel=panels.txt      │ ║
║ │      gt=query.vcf map=genetic_map.txt out=output                       │ ║
║ │ Output: .anc.vcf.gz, .global.anc.gz, .model, .log                     │ ║
║ └─────────────────────────────────────────────────────────────────────────┘ ║
║                                                                              ║
║ ┌─────────────────────────────────────────────────────────────────────────┐ ║
║ │ G-NOMIX                                                                 │ ║
║ │ python3 gnomix.py query.vcf output_dir ${chr} False model.pkl          │ ║
║ │ Or train: python3 gnomix.py query.vcf genetic_map output_dir ${chr}    │ ║
║ │           False reference.vcf sample_map.txt                           │ ║
║ │ Output: .msp, .fb, ancestry plots                                      │ ║
║ └─────────────────────────────────────────────────────────────────────────┘ ║
║                                                                              ║
└═════════════════════════════════════════════════════════════════════════════┘
              │
              ▼
┌═════════════════════════════════════════════════════════════════════════════┐
║ PROCESS: SUMMARIZE_LAI                                                       ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ Combine chromosome-level LAI results                                         ║
║ Generate genome-wide ancestry tracts                                         ║
║ Create summary-level LAI visualizations (no per-sample karyograms)          ║
║ Calculate per-ancestry statistics                                            ║
║                                                                              ║
║ NOTE: Per-sample karyograms are NOT generated to support large cohorts      ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ Output: Summary files + cohort-level visualizations                          ║
└═════════════════════════════════════════════════════════════════════════════┘
              │
              ▼
                         FINAL OUTPUTS
```

### Imputation Performance by Ancestry

Module 7 includes ancestry-stratified imputation quality analysis:

```
┌═════════════════════════════════════════════════════════════════════════════┐
║ PROCESS: IMPUTATION_PERFORMANCE_BY_ANCESTRY                                  ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ Stratifies imputation quality metrics by GRAF-anc ancestry groups            ║
║                                                                              ║
║ Metrics Tracked:                                                             ║
║   • INFO score (from imputation server)                                      ║
║   • R² score (imputation accuracy estimate)                                  ║
║   • MagicalRsq-X calibrated scores (if enabled)                              ║
║   • Variant counts by MAF bin                                                ║
║                                                                              ║
║ Stratification Levels:                                                       ║
║   • 8 Major Groups: AFR, MEN, EUR, SAS, EAS, AMR, OCN, MIX                  ║
║   • 50+ Subgroups: Detailed subcontinental populations                       ║
║                                                                              ║
║ Output:                                                                      ║
║   • Summary statistics per ancestry group                                    ║
║   • PDF plots comparing INFO/R² distributions                               ║
║   • Tab-delimited metrics files                                              ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ Location: ${outdir}/module7/06_imputation_performance/${server}/             ║
└═════════════════════════════════════════════════════════════════════════════┘
```

### WGS Truth Data Comparison (Optional)

For validation studies, the pipeline outputs can be compared against WGS truth data:

```
Recommended Validation Workflow:
────────────────────────────────
1. Obtain WGS truth data for subset of samples
   • HGDP-1KG high-coverage WGS (~4000 samples)
   • Study-specific WGS (if available)

2. Match variants between imputed and WGS data
   • bcftools isec for intersection
   • Ensure consistent variant IDs (chr:pos:ref:alt)

3. Calculate concordance metrics
   • Genotype concordance (GT match rate)
   • Dosage R² (correlation of imputed dosage vs truth)
   • Non-reference concordance (for rare variants)

4. Stratify by:
   • MAF bins (0-0.01, 0.01-0.05, 0.05-0.5)
   • Ancestry group (from GRAF-anc)
   • Imputation server (TOPMed vs All of Us)

5. Compare MagicalRsq-X calibration
   • Pre-filter INFO/R² threshold selection
   • Post-filter variant retention
```

Key metrics for publication:
| Metric | Description | Expected Range |
|--------|-------------|----------------|
| Dosage R² | Correlation with WGS dosage | >0.9 (common), >0.7 (low-freq) |
| Concordance | Exact genotype match | >95% overall |
| NRC | Non-reference concordance | >90% for MAF >1% |
| Variant Yield | Variants passing QC | ~8-40M depending on MAF |

---

## LOCAL ANCESTRY INFERENCE TOOLS

### Comparison

| Tool | Algorithm | Speed | Accuracy | Input Format |
|------|-----------|-------|----------|--------------|
| RFMix v2 | Random Forest + CRF | Medium | High | VCF/BCF |
| RFMix v1 | Random Forest | Slower | High | Binary alleles |
| FLARE | HMM + EM | Fast | High | VCF |
| G-NOMIX | Neural Network | Fastest | High | VCF |

### Required Reference Data

All LAI tools require:
1. **Reference VCF**: Phased haplotypes from reference populations
2. **Sample Map**: Assigns reference samples to ancestral populations
3. **Genetic Map**: Recombination rates (cM/Mb)

### RFMix v2 (Default)

```bash
rfmix \
    -f query.vcf.gz \
    -r reference.vcf.gz \
    -m sample_map.txt \
    -g genetic_map.txt \
    -o output_prefix \
    --chromosome=22 \
    -e 5 \
    --reanalyze-reference \
    --n-threads=8
```

**Sample Map Format (tab-separated):**
```
sample1    EUR
sample2    AFR
sample3    EAS
```

**Genetic Map Format:**
```
chromosome    position    rate    cM
22           16050075    0       0.0
22           16050115    0.0001  0.00001
```

### FLARE

```bash
java -Xmx16g -jar flare.jar \
    ref=reference.vcf.gz \
    ref-panel=panels.txt \
    gt=query.vcf.gz \
    map=genetic_map.txt \
    out=output_prefix \
    nthreads=8 \
    em=true
```

**Panel File Format:**
```
sample1    EUR
sample2    AFR
```

### G-NOMIX

```bash
# Using pre-trained model
python3 gnomix.py \
    query.vcf.gz \
    output_dir/ \
    22 \
    False \
    pretrained_model.pkl

# Training new model
python3 gnomix.py \
    query.vcf.gz \
    genetic_map.txt \
    output_dir/ \
    22 \
    False \
    reference.vcf.gz \
    sample_map.txt
```

---

## CONFIGURATION REFERENCE

### Complete Parameter List

```groovy
params {
    // Input/Output
    sample_sheet                = null
    outdir                      = 'results'
    skip_modules                = ''

    // Reference Files
    hg19_fasta                  = 'resources/references/hg19.fa'
    hg38_fasta                  = 'resources/references/hg38.fa'
    liftover_chain              = 'resources/references/hg19ToHg38.over.chain.gz'
    topmed_reference            = 'resources/references/PASS.Variants.TOPMed_freeze10_hg38.tab.gz'

    // Imputation Services
    run_topmed                  = true
    run_anvil                   = false
    topmed_api_token            = null
    topmed_password             = null
    monitor_interval_minutes    = 45

    // MagicalRsq-X Configuration
    magicalrsqx_dir             = '/opt/MagicalRsqX'          // Container path
    magicalrsqx_models          = '/opt/MagicalRsqX/models'   // Pre-trained XGBoost models
    primary_ancestry            = 'EUR'                       // Model selection: EUR, AFR, AMR, EAS
    magicalrsq_threshold        = 0.3                         // R² filter threshold

    // QC Thresholds
    sample_call_rate            = 0.95
    variant_call_rate           = 0.95
    skip_hwe                    = true       // Skip HWE filtering (recommended for admixed)
    hwe_pvalue                  = 1e-6       // HWE threshold if skip_hwe=false
    het_sd_threshold            = 3
    kinship_threshold           = 0.177
    use_genesis                 = true
    n_pcs                       = 15         // Number of PCs to calculate
    maf_threshold               = 0.01       // MAF filter for GWAS output

    // Ancestry - Global
    admixture_k_min             = 5          // ADMIXTURE K range start
    admixture_k_max             = 7          // ADMIXTURE K range end
    admixture_cv_folds          = 5          // Cross-validation folds

    // Ancestry - Local (LAI)
    run_lai                     = false
    lai_methods                 = 'rfmix2'   // Options: rfmix2, rfmix1, flare, gnomix
    lai_reference               = 'resources/ancestry_references/lai_reference'
    genetic_map_dir             = 'resources/genetic_maps'
}
```

---

## TROUBLESHOOTING

### Common Issues

| Issue | Cause | Solution |
|-------|-------|----------|
| "terralab not authenticated" | OAuth tokens expired | Run `terralab login` again |
| "TOPMed job queue full" | Server limits | Wait or reduce concurrent jobs |
| Merge conflicts | Strand issues | Module 4 auto-fixes most cases |
| Low liftover rate | Wrong input build | Verify `build` column in sample sheet |
| LAI crashes | Insufficient memory | Increase memory allocation |

### Log Files

- Pipeline logs: `results/pipeline_info/`
- Per-process logs: `work/*/` directories
- Execution trace: `results/pipeline_info/execution_trace.txt`

---

**Document Version:** 1.0
**Last Updated:** December 2025
