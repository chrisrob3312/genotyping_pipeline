# Module 7: Ancestry Analysis - Visual Workflow Guide

## 📊 Quick Reference

| Component | Status | Purpose |
|-----------|--------|---------|
| GRAF-anc |   Global & subcontinental ancestry assignment |
| ADMIXTURE | Admixture modeling (K=2-12) |
| ADMIXTURE by Ancestry  | K=6 & K=9 organized by GRAF-anc groups |
| PCA Colored | Visualize PCA with ancestry colors |
| QC by Ancestry | Compare quality metrics across groups |
| Imputation Performance | Compare imputation quality by ancestry |
| Local Ancestry (LAI) |  Chromosome-level ancestry (RFMix, etc.) |

---

## 🔄 Complete Workflow Diagram

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         MODULE 7 INPUTS                                   │
│                                                                           │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐                  │
│  │  QC'd Data   │  │  PCA Results │  │ Imputed VCFs │                  │
│  │  (Module 6)  │  │  (Module 6)  │  │  (Servers)   │                  │
│  │ .bed/.bim/.fam│  │ .eigenvec/val│  │    .vcf.gz   │                  │
│  └──────┬───────┘  └──────┬───────┘  └──────┬───────┘                  │
└─────────┼──────────────────┼──────────────────┼───────────────────────────┘
          │                  │                  │
          ▼                  │                  │
┌─────────────────────┐     │                  │
│  PREPARE_FOR_       │     │                  │
│  ANCESTRY           │     │                  │
│  • LD pruning       │     │                  │
│  • Format data      │     │                  │
└─────┬───────────┬───┘     │                  │
      │           │         │                  │
      │           └─────────┼──────────┐       │
      │                     │          │       │
      ▼                     ▼          ▼       │
┌─────────────┐   ┌─────────────┐  ┌──────────────┐
│  GRAFANC_   │   │ RUN_        │  │ PCA_BY_      │
│  ANCESTRY   │   │ ADMIXTURE   │  │ ANCESTRY     │
│             │   │             │  │              │
│ Continental │   │ K=2 to K=12 │  │ Color by     │
│ & Subcon-   │   │ (parallel)  │  │ GRAF-anc     │
│ tinental    │   │             │  │              │
│ assignment  │   │ Cross-      │  │ PC1/PC2      │
│             │   │ validation  │  │ PC3/PC4      │
│ Major: 1-8  │   └──────┬──────┘  └──────────────┘
│ Sub: 101-800│          │
└──────┬──────┘          │
       │                 │
       │                 ▼
       │        ┌─────────────────┐
       │        │ SUMMARIZE_      │
       │        │ ADMIXTURE       │
       │        │                 │
       │        │ • CV plot       │
       │        │ • Standard      │
       │        │   barplots      │
       │        │ • Optimal K     │
       │        └────────┬────────┘
       │                 │
       ├─────────────────┤
       │                 │
       ▼                 ▼
┌──────────────────────────────────┐
│  ADMIXTURE_BY_GRAFANC            │
│  🆕 NEW FEATURE                  │
│                                  │
│  K=6 & K=9 plots organized by:   │
│  • GRAF-anc major groups         │
│  • Sorted by dominant component  │
│  • Creates "waterfall" effect    │
│                                  │
│  [AFR] [EUR] [EAS] [SAS] [AMR]  │
│   ████  ████  ████  ████  ████   │
└──────────────────────────────────┘
       │
       ├────────────────────────────────────────────┐
       │                                            │
       ▼                                            ▼
┌──────────────────┐                    ┌──────────────────────┐
│ QC_BY_           │                    │ IMPUTATION_          │
│ ANCESTRY         │                    │ PERFORMANCE_BY_      │
│                  │                    │ ANCESTRY             │
│ • Call rates     │                    │ 🆕 NEW FEATURE       │
│ • Missingness    │                    │                      │
│ • Sample counts  │                    │ • INFO/R² scores     │
│   by group       │                    │ • Quality by chr     │
│                  │                    │ • Sample distribution│
└──────────────────┘                    └──────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│              OPTIONAL: LOCAL ANCESTRY INFERENCE             │
│              (Enable with params.run_lai = true)            │
│                                                             │
│  ┌────────────────┐                                        │
│  │ PREPARE_FOR_   │                                        │
│  │ LAI            │                                        │
│  │ • Check phasing│                                        │
│  └────────┬───────┘                                        │
│           │                                                │
│           ▼                                                │
│  ┌────────────────┐    Run for each chromosome            │
│  │ RFMIX_V2 or    │    (chr 1-22)                         │
│  │ other LAI tool │    ════════════════                   │
│  │                │                                        │
│  │ • RFMix v2     │◄── Parallel execution                 │
│  │ • RFMix v1     │                                        │
│  │ • FLARE        │                                        │
│  │ • G-NOMIX      │                                        │
│  └────────┬───────┘                                        │
│           │                                                │
│           ▼                                                │
│  ┌────────────────┐                                        │
│  │ SUMMARIZE_LAI  │                                        │
│  │ • Combine chrs │                                        │
│  │ • Plot results │                                        │
│  └────────────────┘                                        │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│                        MODULE 7 OUTPUTS                     │
│                                                             │
│  📁 01_grafanc/          📁 02_admixture/                  │
│  📁 04_pca_colored/      📁 05_qc_by_ancestry/             │
│  📁 06_imputation_performance/ (NEW)                        │
│  📁 03_local_ancestry/   (if LAI enabled)                  │
└─────────────────────────────────────────────────────────────┘
```

---

##  Process Execution Order

### Phase 1: Data Preparation (Sequential)
```
PREPARE_FOR_ANCESTRY
├── Input:  QC'd genotypes (bed/bim/fam)
├── Action: LD pruning (200 50 0.25)
└── Output: Pruned data for ancestry analysis
    Time: ~5-15 minutes
```

### Phase 2: Core Ancestry Analysis (Parallel)

#### Branch A: GRAF-anc (Critical Path)
```
GRAFANC_ANCESTRY
├── Input:  Pruned genotypes
├── Action: Compare to 282K ancestry SNP panel
│           Calculate GD1-GD6 scores
│           Assign major (1-8) and subgroup (101-800) codes
├── Output: • grafanc_results.txt (full data)
│           • by_major_group.txt (continental counts)
│           • by_subgroup.txt (subcontinental counts)
│           • GD_plots.pdf (visualizations)
└── Time: ~10-30 minutes (depends on sample size)
```

#### Branch B: ADMIXTURE (Parallel Execution)
```
RUN_ADMIXTURE (K=2 to K=12)
├── Runs in parallel for each K
├── K=2  ─┐
├── K=3   │
├── K=4   ├─> Parallel execution
├── K=5   │   (up to 12 concurrent jobs)
├── ...   │
└── K=12 ─┘
    Time: ~30 minutes to 2 hours per K
          (all K run simultaneously)
    
    ↓
    
SUMMARIZE_ADMIXTURE
├── Collect all K results
├── Generate CV plot (identify optimal K)
└── Create standard barplots for each K
    Time: ~5-10 minutes
```

### Phase 3: Integration & Organization (Sequential, depends on Phases 1-2)

#### Process A: ADMIXTURE by Ancestry (NEW)
```
ADMIXTURE_BY_GRAFANC
├── Requires: GRAF-anc results + ADMIXTURE K=6,9 results
├── Action:  1. Group samples by GRAF-anc major ancestry
│            2. Sort by dominant ADMIXTURE component
│            3. Create organized "waterfall" plots
├── Output:  • K6 plot organized by ancestry
│            • K9 plot organized by ancestry
│            • Average proportions per group
└── Time: ~5 minutes
```

#### Process B: PCA Coloring
```
PCA_BY_ANCESTRY
├── Requires: PCA results + GRAF-anc results
├── Action:  Color PC1/PC2, PC3/PC4 by ancestry
├── Output:  pca_by_ancestry.pdf
└── Time: ~2 minutes
```

#### Process C: QC Comparison
```
QC_BY_ANCESTRY
├── Requires: Original genotypes + GRAF-anc results
├── Action:  Calculate call rates per ancestry group
├── Output:  • qc_by_ancestry.txt
│            • qc_by_ancestry.pdf
└── Time: ~5-10 minutes
```

#### Process D: Imputation Performance (NEW)
```
IMPUTATION_PERFORMANCE_BY_ANCESTRY
├── Requires: Imputed VCFs + GRAF-anc results
├── Action:  1. Extract INFO/R² scores
│            2. Calculate stats by chromosome
│            3. Count samples per ancestry
├── Output:  • imputation_by_ancestry_summary.txt
│            • imputation_by_ancestry_plots.pdf
│            • metrics by major/subgroup
└── Time: ~10-20 minutes
```

### Phase 4: Local Ancestry (Optional, Parallel)

```
IF params.run_lai = true:
    
    PREPARE_FOR_LAI
    ├── Check if VCF is phased
    └── Format for LAI tools
        Time: ~5 minutes
        
        ↓
        
    RFMIX_V2 (or other LAI tool)
    ├── Chr 1  ─┐
    ├── Chr 2   │
    ├── Chr 3   ├─> Parallel execution
    ├── ...     │   (22 chromosomes)
    └── Chr 22 ─┘
        Time: ~1-4 hours per chromosome
              (all run simultaneously)
        
        ↓
        
    SUMMARIZE_LAI
    ├── Combine chromosome results
    └── Create genome-wide plots
        Time: ~10-15 minutes
```

---

##  Output Structure (Annotated)

```
results/module7/
│
├── 00_prep/                              [Phase 1]
│   └── [server]/
│       └── *_pruned.{bed,bim,fam}       # LD-pruned data
│
├── 01_grafanc/                           [Phase 2A - CRITICAL]
│   └── [server]/
│       ├── *_grafanc_results.txt          #Main ancestry file
│       │                                  # (Sample, #SNPs, GD1-6, Pe/Pf/Pa, AncGroupID)
│       ├── *_grafanc_by_major_group.txt   #Continental counts
│       │                                  # (AFR=1XX, EUR=3XX, EAS=5XX, etc.)
│       ├── *_grafanc_by_subgroup.txt     #Subcontinental counts
│       │                                   #(101=Nigeria, 305=NE_Europe, etc.)
│       ├── *_grafanc_summary.txt         # Human-readable summary
│       └── *_grafanc_GD_plots.pdf         #Population structure plots
│
├── 02_admixture/                         [Phase 2B + 3A]
│   └── [server]/
│       ├── K2/ ... K12/                  # Raw ADMIXTURE outputs
│       │   ├── *.Q                       # Ancestry proportions
│       │   └── *.P                       # Allele frequencies
│       │
│       ├── *_admixture_summary.txt        #CV errors & optimal K
│       ├── *_admixture_cv_plot.pdf        #Cross-validation plot
│       ├── *_admixture_K2.pdf             #Standard barplot K=2
│       ├── *_admixture_K3.pdf             #Standard barplot K=3
│       ├── ...                           # Through K=12
│       │
│       └── by_ancestry/                 
│           ├── *_by_ancestry_K6.pdf      #K=6 organized by GRAF-anc
│           │                               # [AFR][EUR][EAS] groups
│           ├── *_by_ancestry_K9.pdf      #K=9 organized by GRAF-anc
│           └── *_ancestry_summary.txt    #Avg proportions per group
│
├── 03_local_ancestry/                    [Phase 4 - OPTIONAL]
│   ├── 00_prep/
│   ├── rfmix_v2/                        # If LAI enabled
│   └── *_lai_summary.txt
│
├── 04_pca_colored/                       [Phase 3B]
│   └── [server]/
│       └── *_pca_by_ancestry.pdf         # PC1/PC2, PC3/PC4 colored
│
├── 05_qc_by_ancestry/                    [Phase 3C]
│   └── [server]/
│       ├── *_qc_by_ancestry.txt          # Call rates by group
│       └── *_qc_by_ancestry.pdf          # QC visualizations
│
└── 06_imputation_performance/            [Phase 3D - NEW]
    └── [server]/
        ├── *_imputation_by_ancestry_summary.txt   Quality summary
        ├── *_imputation_by_ancestry_plots.pdf     Quality plots
        ├── *_metrics_by_major_group.txt           Sample counts
        └── *_metrics_by_subgroup.txt              Detailed counts


```

---

##  Branching & Dependencies

### Critical Path (Must Complete)
```
QC Data → PREPARE_FOR_ANCESTRY → GRAFANC_ANCESTRY → [All downstream processes]
```

### Parallel Branches (Can run simultaneously after PREPARE_FOR_ANCESTRY)
```
Branch 1: GRAFANC_ANCESTRY
Branch 2: RUN_ADMIXTURE (K=2,3,4...12 all parallel)
```

### Integration Branches (Require multiple inputs)
```
ADMIXTURE_BY_GRAFANC
├── Requires: GRAFANC_ANCESTRY output
└── Requires: RUN_ADMIXTURE K=6,9 outputs

PCA_BY_ANCESTRY
├── Requires: GRAFANC_ANCESTRY output
└── Requires: PCA input from Module 6

QC_BY_ANCESTRY
├── Requires: GRAFANC_ANCESTRY output
└── Requires: Original genotype data

IMPUTATION_PERFORMANCE_BY_ANCESTRY
├── Requires: GRAFANC_ANCESTRY output
└── Requires: Imputed VCF files
```

---

## ⏱️ Runtime Estimates

| Process | Small Dataset<br>(100 samples) | Medium Dataset<br>(1,000 samples) | Large Dataset<br>(10,000 samples) |
|---------|---------|----------|---------|
| **PREPARE_FOR_ANCESTRY** | 5 min | 10 min | 30 min |
| **GRAFANC_ANCESTRY** | 10 min | 20 min | 60 min |
| **RUN_ADMIXTURE (per K)** | 30 min | 2 hours | 8 hours |
| **SUMMARIZE_ADMIXTURE** | 2 min | 5 min | 10 min |
| **ADMIXTURE_BY_GRAFANC** | 2 min | 5 min | 10 min |
| **PCA_BY_ANCESTRY** | 1 min | 2 min | 5 min |
| **QC_BY_ANCESTRY** | 5 min | 10 min | 30 min |
| **IMPUTATION_PERFORMANCE** | 5 min | 15 min | 45 min |
| **RFMIX_V2 (per chr)** | 30 min | 2 hours | 6 hours |
| **Total (without LAI)** | ~2 hours | ~6 hours | ~24 hours |
| **Total (with LAI)** | ~12 hours | ~48 hours | ~144 hours |

*Note: ADMIXTURE K values run in parallel, so total time ≈ time for single K*

---

##  Configuration Options

### Required Parameters
```groovy
params {
    outdir = 'results'                # Output directory
    admixture_k_min = 2              # Minimum K to test
    admixture_k_max = 12             # Maximum K to test
}
```

### Optional Parameters
```groovy
params {
    admixture_cv_folds = 5           # Cross-validation folds
    run_lai = false                  # Enable Local Ancestry Inference
    lai_tool = 'rfmix_v2'           # LAI tool to use
}
```

### Process Resources (adjust for your system)
```groovy
process {
    withLabel: 'process_medium' {
        cpus = 4
        memory = 16.GB
    }
    withLabel: 'process_high' {
        cpus = 8
        memory = 32.GB
    }
}
```

---

##  Decision Points

### 1. Run Local Ancestry Inference?
```
params.run_lai = true  → Execute LAI branch (Phase 4)
params.run_lai = false → Skip LAI (default)
```

### 2. Which LAI Tool?
```
params.lai_tool = 'rfmix_v2'  → Use RFMix v2 (recommended)
params.lai_tool = 'rfmix_v1'  → Use RFMix v1
params.lai_tool = 'flare'     → Use FLARE
params.lai_tool = 'gnomix'    → Use G-NOMIX
params.lai_tool = 'all'       → Run all available tools
```

---

##  Process Details

### GRAFANC_ANCESTRY
```
Input:   Pruned PLINK files
Tool:    grafanc executable
Command: grafanc input_prefix output.txt --threads N --maxmem M
Output:  Tab-delimited file with columns:
         Sample, #SNPs, GD1-GD3, EA1-EA4, AF1-AF3, EU1-EU3,
         SA1-SA2, IC1-IC3, Pe, Pf, Pa, RawPe, RawPf, RawPa,
         AncGroupID (3-digit code)
```

### RUN_ADMIXTURE
```
Input:   Pruned PLINK .bed file
Tool:    ADMIXTURE
Command: admixture --cv=5 input.bed K -jN
Output:  • .Q file (ancestry proportions)
         • .P file (allele frequencies)
         • CV error for this K
Runs:    Once per K value (in parallel)
```

### ADMIXTURE_BY_GRAFANC (New)
```
Input:   • GRAF-anc results
         • ADMIXTURE .Q files for K=6 and K=9
         • FAM file (sample IDs)
Process: 1. Merge GRAF-anc ancestry with ADMIXTURE
         2. Group by major ancestry (AFR, EUR, EAS, etc.)
         3. Within each group, sort by dominant component
         4. Create faceted barplots
Output:  • Organized K=6 plot
         • Organized K=9 plot
         • Summary statistics
```

### IMPUTATION_PERFORMANCE_BY_ANCESTRY (New)
```
Input:   • Imputed VCF files
         • GRAF-anc major/subgroup files
Tool:    bcftools query + R
Extract: INFO or R² scores from VCF
Analyze: • Overall quality distribution
         • Quality by chromosome
         • Sample counts by ancestry
Output:  • Summary text file
         • Quality plots (histograms, boxplots)
         • Metrics by major/subgroup
```

---

##  Key Files for Downstream Analysis

### For Stratified GWAS
```bash
# Use GRAF-anc major groups
*_grafanc_by_major_group.txt

# Extract specific ancestry
awk '$2 == "EUR"' *_grafanc_results.txt > european_samples.txt
```

### For Multi-Ancestry Meta-Analysis
```bash
# Sample counts per group
*_grafanc_by_major_group.txt
*_grafanc_by_subgroup.txt

# Ancestry proportions as covariates
awk '{print $1,$1,$14,$15,$16}' *_grafanc_results.txt > ancestry_pcs.txt
# Columns: FID, IID, Pe, Pf, Pa
```

### For Quality Control
```bash
# Check imputation quality by ancestry
*_imputation_by_ancestry_summary.txt

# Identify poor quality samples by ancestry
*_qc_by_ancestry.txt
```

### For Admixture Studies
```bash
# ADMIXTURE proportions organized by ancestry
*_admixture_by_ancestry_K6.pdf
*_admixture_ancestry_summary.txt

# Individual-level admixture
*.6.Q  # K=6 proportions for each sample
```

---

##  Execution Flow Summary

```
START
  │
  ├─> Phase 1: Prepare (Sequential)
  │   └─> LD prune data
  │
  ├─> Phase 2: Core Analysis (Parallel)
  │   ├─> GRAF-anc (Continental + Subcontinental)
  │   └─> ADMIXTURE (K=2 to K=12, all parallel)
  │
  ├─> Phase 3: Integration (Sequential, needs Phase 2)
  │   ├─> ADMIXTURE organized by ancestry (NEW)
  │   ├─> PCA colored by ancestry
  │   ├─> QC metrics by ancestry
  │   └─> Imputation performance by ancestry (NEW)
  │
  └─> Phase 4: LAI (Optional, Parallel)
      └─> Local ancestry per chromosome
  
END → All outputs in results/module7/
```

---

##  Quick Troubleshooting

| Issue | Check | Solution |
|-------|-------|----------|
| grafanc fails | `which grafanc` | Add to PATH or set GRAFPATH |
| Few ancestry SNPs | Check log for SNP count | Update variant IDs to rsIDs |
| ADMIXTURE diverges | Check log files | Reduce K range or CV folds |
| Missing imputation metrics | VCF header | Ensure INFO/R² field present |
| Memory errors | Process logs | Increase memory allocation |

---

##  Quick Reference Commands

```bash
# Check grafanc is available
grafanc --help

# View ancestry distribution
cat results/module7/01_grafanc/*/​*_by_major_group.txt

# Find optimal K
grep -h "CV error" results/module7/02_admixture/*/​*.txt | sort -k3 -n | head -1

# Extract European samples
awk 'NR==1 || $NF ~ /^3/' *_grafanc_results.txt > europeans.txt

# Check imputation quality
cat results/module7/06_imputation_performance/*/​*_summary.txt
```

---

**Last Updated:** January 2025  
**Version:** 1.0 (Visual Workflow Edition)
