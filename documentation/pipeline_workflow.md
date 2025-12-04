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
│ • Karyogram plots (for LAI)                                                  │
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

See `documentation/Module_1_Workflow.md` for detailed process flow.

### Summary

```
INPUT                    PROCESS                           OUTPUT
─────                    ───────                           ──────
Sample Sheet  ───────►  Discover Files
                              │
                              ▼
                        Prepare Samples
                              │
                              ▼
                        Merge to Batch (UNION)
                              │
                              ▼
                        Align to Reference
                              │
                              ▼
              ┌───────────────┴───────────────┐
              ▼                               ▼
        hg19 batches                    hg38 batches
              │                               │
              ▼                               │
        Liftover (CrossMap)                   │
              │                               │
              └───────────────┬───────────────┘
                              ▼
                        All hg38
                              │
                              ▼
                        Merge to Platform (UNION)
                              │
                              ▼
              ┌───────────────┴───────────────┐
              ▼                               ▼
        TOPMed Validation              AnVIL Validation
              │                               │
              ▼                               ▼
        Light QC                        Light QC
              │                               │
              ▼                               ▼
        Create TOPMed VCFs        Create AnVIL VCF  ───► Single VCF
        (22 per platform)                               (all autosomes)
```

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sample_call_rate` | 0.95 | Minimum sample call rate |
| `variant_call_rate` | 0.95 | Minimum variant call rate |

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
║ PROCESS: HWE_FILTER                                                          ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ plink2 --hwe ${params.hwe_pvalue} midp                                       ║
║ Remove variants out of HWE (p < 1e-6)                                        ║
║ Using mid-p correction for small samples                                     ║
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
| `hwe_pvalue` | 1e-6 | HWE filter p-value threshold |
| `het_sd_threshold` | 3 | Heterozygosity SD threshold |
| `use_genesis` | true | Use GENESIS PCRelate (vs PLINK) |
| `kinship_threshold` | 0.177 | 1st-degree relatedness cutoff |

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
║ Create karyogram visualizations                                              ║
║ Calculate per-ancestry statistics                                            ║
╠─────────────────────────────────────────────────────────────────────────────╢
║ Output: Summary files + visualizations                                       ║
└═════════════════════════════════════════════════════════════════════════════┘
              │
              ▼
                         FINAL OUTPUTS
```

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

    // QC Thresholds
    magicalrsq_threshold        = 0.3
    sample_call_rate            = 0.95
    variant_call_rate           = 0.95
    hwe_pvalue                  = 1e-6
    het_sd_threshold            = 3
    kinship_threshold           = 0.177
    use_genesis                 = true

    // Ancestry
    admixture_k                 = '5,6,7'
    run_lai                     = false
    lai_methods                 = 'rfmix2'
    reference_panel             = 'resources/ancestry_references/reference_panel.rds'
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
