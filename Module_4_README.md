# Module 4: Platform Merging - Visual Workflow Guide

## Overview

**Goal:** Merge multi-platform genotype data with 3-pass strategy  
**Input:** Module 3 QC'd VCFs (multiple platforms)  
**Output:** Single merged dataset per reference panel  
**QC:** SNPs only, bi-allelic, no monomorphic, call rates >95%

---


## 🔬 Process Summary

| # | Process | Input | Output | Time | Memory |
|---|---------|-------|--------|------|--------|
| 0 | VCF→PLINK | VCF.gz | .pgen/.pvar/.psam | 2 min | 8 GB |
| 1a | Pass 1: Identify | PLINK files | Error log | 3 min | 16 GB |
| 1b | Pass 2: Analyze | Error log | Fix lists | 1 min | 4 GB |
| 1c | Pass 2: Apply | PLINK + fixes | Fixed PLINK | 2 min | 8 GB |
| 1d | Pass 3: Merge | Fixed PLINK | Merged | 5 min | 24 GB |
| 2 | SNPs only | Merged | Biallelic | 2 min | 8 GB |
| 3 | Remove mono | Biallelic | MAF>0 | 2 min | 8 GB |
| 4 | Variant CR | MAF>0 | CR>95% | 2 min | 8 GB |
| 5 | Sample CR | Variant CR | CR>95% | 2 min | 8 GB |
| 6 | PLINK→VCF | QC'd PLINK | VCF.gz | 3 min | 8 GB |
| 7 | QC Report | All logs | HTML | 1 min | 4 GB |
| **Total** | | | | **25 min** | **24 GB peak** |

---

## 🛠️ QC Steps Detail

```
┌─────────────────────────────────────────┐
│ 1. SNPs Only (Biallelic)               │
│    ├─ Remove indels                    │
│    ├─ Remove multi-allelic             │
│    └─ Keep: A/T, C/G only              │
│                                         │
│ 2. Polymorphic (MAF > 0)               │
│    ├─ Calculate allele frequencies     │
│    └─ Remove monomorphic variants      │
│                                         │
│ 3. Variant Call Rate (>95%)            │
│    ├─ Calculate per-variant missingness│
│    └─ Remove variants with <95% calls  │
│                                         │
│ 4. Sample Call Rate (>95%)             │
│    ├─ Calculate per-sample missingness │
│    └─ Remove samples with <95% calls   │
└─────────────────────────────────────────┘
```
---

## 🔄 3-Pass Strategy Explained

```
┌──────────────────────────────────────────────────────┐
│  PASS 1: IDENTIFY                                    │
├──────────────────────────────────────────────────────┤
│  Platform A: rs123  REF=A  ALT=T                     │
│  Platform B: rs123  REF=T  ALT=A  ← MISMATCH!        │
│                                                       │
│  Action: Let merge FAIL, capture error message       │
│  Result: List of 1,500 problematic variants          │
└──────────────────────────────────────────────────────┘
                         ↓
┌──────────────────────────────────────────────────────┐
│  PASS 2: FIX                                         │
├──────────────────────────────────────────────────────┤
│  Analyze each mismatch:                              │
│                                                       │
│  Type 1: Strand Flip (85%)                           │
│    A ↔ T → Apply --flip                              │
│                                                       │
│  Type 2: REF/ALT Swap (10%)                          │
│    REF=A/ALT=T vs REF=T/ALT=A                        │
│    → Apply --ref-allele                              │
│                                                       │
│  Type 3: Both (4%)                                   │
│    → Apply --flip + --ref-allele                     │
│                                                       │
│  Type 4: Unfixable (1%)                              │
│    → Apply --exclude                                 │
│                                                       │
│  Result: 1,450 fixed, 50 excluded                    │
└──────────────────────────────────────────────────────┘
                         ↓
┌──────────────────────────────────────────────────────┐
│  PASS 3: CLEAN MERGE                                 │
├──────────────────────────────────────────────────────┤
│  Merge with corrected platforms                      │
│                                                       │
│  Input:  2,000,000 variants                          │
│  Merged: 1,999,950 variants (99.9975%)               │
│  Saved:  1,450 variants vs traditional approach!     │
│                                                       │
│  Result: ✅ SUCCESS                                  │
└──────────────────────────────────────────────────────┘
```

---

## 🚀 Resource Scaling

### Memory Allocation (Dynamic)

```groovy
// Base memory + scaling per retry attempt
memory { BASE_GB + (SCALE_GB * Math.min(task.attempt, 3)) }
```

| Process | Base | Scale | Max (3 retries) |
|---------|------|-------|-----------------|
| Convert | 8 GB | +2 GB | 14 GB |
| Identify | 16 GB | +4 GB | 28 GB |
| Analyze | 4 GB | +1 GB | 7 GB |
| Fix | 8 GB | +2 GB | 14 GB |
| Merge | 24 GB | +8 GB | 48 GB |
| QC | 8 GB | +2 GB | 14 GB |

### CPU Allocation

```
Standard: 4 cores
Merge:    8 cores
Report:   2 cores
```

### Sample Size Scaling

| Samples | Platforms | Peak Memory | Recommended |
|---------|-----------|-------------|-------------|
| 100 | 2 | 16 GB | 32 GB |
| 500 | 4 | 24 GB | 48 GB |
| 1,000 | 4 | 32 GB | 64 GB |
| 5,000 | 6 | 48 GB | 96 GB |
| 10,000 | 10 | 64 GB | 128 GB |

---
## 📁 Output Structure

```
module4/
├── topmed_ref/
│   ├── 00_vcf_to_plink/        Platform files
│   ├── 01_merge_pass1/         Mismatch logs
│   ├── 02_mismatch_analysis/   Fix lists
│   ├── 03_fixed_platforms/     Corrected files
│   ├── 04_merged_final/        Pass 3 merge
│   ├── 05_snps_only/           Biallelic SNPs
│   ├── 06_polymorphic/         MAF > 0
│   ├── 07_variant_callrate/    Variant CR >95%
│   ├── 08_sample_callrate/     Sample CR >95%
│   ├── 09_final_vcf/           ⭐ READY FOR MODULE 5
│   └── 10_qc_reports/          HTML reports
│
└── allofus_ref/                (Same structure)
```

---

## ⚙️ Configuration

### nextflow.config

```groovy
params {
    outdir = 'results'
    variant_call_rate = 0.05    // 95% call rate
    sample_call_rate = 0.05     // 95% call rate
}

process {
    withLabel: 'plink2' {
        container = 'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'
        cpus = 4
        memory = { 8.GB + (2.GB * task.attempt) }
        errorStrategy = 'retry'
        maxRetries = 3
    }
    
    withLabel: 'plink_merge' {
        container = 'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'
        cpus = 8
        memory = { 16.GB + (4.GB * task.attempt) }
        errorStrategy = 'retry'
        maxRetries = 3
    }
    
    withLabel: 'plink_large' {
        container = 'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'
        cpus = 8
        memory = { 24.GB + (8.GB * task.attempt) }
        errorStrategy = 'retry'
        maxRetries = 3
    }
    
    withLabel: 'R' {
        container = 'rocker/tidyverse:4.3'
        cpus = 2
        memory = 4.GB
    }
}

apptainer {
    enabled = true
    autoMounts = true
    cacheDir = "$HOME/.apptainer/cache"
}
```

---

##  Expected Results

### Variant Retention Through Pipeline

```
┌────────────────────────┬─────────────┬──────────┐
│ Stage                  │ Variants    │ Retained │
├────────────────────────┼─────────────┼──────────┤
│ Input (4 platforms)    │ 2,000,000   │ 100.0%   │
│ Pass 3: Merged         │ 1,999,950   │ 99.998%  │
│ QC: SNPs only          │ 1,999,900   │ 99.995%  │
│ QC: MAF > 0            │ 1,999,850   │ 99.993%  │
│ QC: Variant CR >95%    │ 1,999,800   │ 99.990%  │
│ QC: Sample CR >95%     │ 1,999,750   │ 99.988%  │
│ FINAL                  │ 1,999,750   │ 99.99%   │
└────────────────────────┴─────────────┴──────────┘

 Traditional (no fixes): 1,998,500 variants (99.5%)
   You SAVE: 1,250 variants per chromosome!
```

---

## 🏃 Quick Start

### 1. Pre-Pull Containers

```bash
apptainer pull docker://quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0
apptainer pull docker://rocker/tidyverse:4.3
```

### 2. Test Single Chromosome

```bash
nextflow run Module4_Merging_3PASS_DUAL_CORRECTED.nf \
    --qc_data "results/module3/*/chr22*" \
    --outdir "test_module4" \
    -with-apptainer
```

### 3. Production Run

```bash
nextflow run Module4_Merging_3PASS_DUAL_CORRECTED.nf \
    --qc_data "results/module3/*" \
    --outdir "results/module4" \
    -with-apptainer \
    -resume
```

---

## Quality Checks

### 1. Verify Merge Success

```bash
grep -r "SUCCESS" results/module4/*/04_merged_final/*/chr*.log
```

**Expected:** All chromosomes show SUCCESS

### 2. Check Final Counts

```bash
for vcf in results/module4/*/09_final_vcf/*/*.vcf.gz; do
    echo "$vcf:"
    echo "  Variants: $(bcftools view -H $vcf | wc -l)"
    echo "  Samples: $(bcftools query -l $vcf | wc -l)"
done
```

**Expected:** >99.9% variant retention

### 3. Review QC Summary

```bash
cat results/module4/topmed_ref/10_qc_reports/chr1/*_qc_summary.txt
```

**Expected:**
- SNPs removed: <1,000 per chromosome
- Monomorphic: <100 per chromosome  
- Low CR variants: <200 per chromosome
- Low CR samples: <5% of total

---

## !!! Troubleshooting

### Issue: Out of Memory

```bash
# Increase memory in nextflow.config
withLabel: 'plink_large' {
    memory = { 32.GB + (8.GB * task.attempt) }
}

# Or process fewer chromosomes at once
nextflow run Module4_Merging.nf \
    --qc_data "results/module3/*/chr{1..5}*"
```

### Issue: Pass 3 Still Fails

```bash
# Check error details
less results/module4/*/04_merged_final/chr*/chr*_pass3.log

# Verify genome builds match
bcftools view -h results/module3/*/*.vcf.gz | grep "^##reference"

# Check variant ID consistency
bcftools query -f '%ID\n' results/module3/*/*.vcf.gz | head -20
```

### Issue: Too Many Variants Removed

```bash
# Check QC thresholds
cat results/module4/*/10_qc_reports/chr1/*_qc_summary.txt

# Adjust if needed (in params)
--variant_call_rate 0.10  # Allow 90% instead of 95%
--sample_call_rate 0.10
```

---

## Performance Monitoring

```bash
# Run with monitoring
nextflow run Module4_Merging.nf \
    -with-trace \
    -with-timeline timeline.html \
    -with-report report.html

# View results
firefox timeline.html  # See timing per process
firefox report.html    # See resource usage
```

---

## 🎓 Key Points

| Concept | Detail |
|---------|--------|
| **Strategy** | 3-pass: Identify → Fix → Merge |
| **QC in M4** | SNPs, MAF>0, call rates |
| **Retention** | >99.99% with fixes |
| **Time** | ~25 min per chromosome |
| **Memory** | 24 GB peak (48 GB for large datasets) |
| **Scaling** | Dynamic retry with +memory |

---

## 🔗 Next Steps

After Module 4 completes:

1.  Verify all QC reports look good
2.  Check final variant and sample counts
3.  Proceed to **Module 5** (additional post-merge QC)
4.  Then **Module 6** (ancestry)

---

**Module 4 provides a script for large-scale multi-platform merging with automatic quality control and intelligent resource scaling!** 
