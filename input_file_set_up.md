# Input File Setup Guide v2.0

## Overview

This pipeline is **flexible** and handles multiple real-world data structures:
- Individual samples (all chromosomes together OR split by chromosome)
- Pre-merged batch files
- Mixed file formats (PLINK, VCF)
- Different naming conventions
- Different genome builds per batch

---

## Sample Sheet Format (CSV)

### Basic Columns (Always Required)

```csv
platform_id,batch_id,input_path,file_type,build,file_structure
```

### Column Definitions

| Column | Description | Required | Values |
|--------|-------------|----------|--------|
| `platform_id` | Genotyping platform/array | Yes | GSAv1, GSAv2, Omni25, MEGA |
| `batch_id` | Batch identifier | Yes | batch_2020_01, cohort_A, phase1 |
| `input_path` | Path to directory OR file prefix | Yes | /data/samples/ OR /data/merged_batch |
| `file_type` | File format | Yes | plink, vcf |
| `build` | Genome build for this batch | Yes | hg19, hg38 |
| `file_structure` | How files are organized | Yes | See below |

### File Structure Options

| `file_structure` | Description | Example Files |
|------------------|-------------|---------------|
| `individual_samples` | One file set per sample, all chr together | sample_001.bed, sample_002.bed |
| `individual_chr_split` | One file per sample per chromosome | sample_001_chr1.bed, sample_001_chr2.bed |
| `merged_batch` | Already merged batch file | batch1_all_samples.bed |
| `merged_chr_split` | Merged batch, split by chromosome | batch1_chr1.bed, batch1_chr2.bed |

---

## Scenario 1: Individual Samples (All Chromosomes Together)

**Most common format** - Each sample is separate, contains all chromosomes.

### Directory Structure:
```
/data/gsa_v1/batch_2020_01/
├── sample_001.bed/bim/fam
├── sample_002.bed/bim/fam
├── sample_003.bed/bim/fam
└── ... (850 samples)
```

### Sample Sheet Entry:
```csv
platform_id,batch_id,input_path,file_type,build,file_structure
GSAv1,batch_2020_01,/data/gsa_v1/batch_2020_01,plink,hg19,individual_samples
```

### What Pipeline Does:
```
1. Finds all *.bed files in directory
2. Processes each sample in parallel (850 parallel jobs)
3. Merges samples into batch
4. Continues with QC, liftover, etc.
```

---

## Scenario 2: Individual Samples Split by Chromosome

**Common after pre-processing** - Each sample has separate files per chromosome.

### Directory Structure:
```
/data/gsa_v2/batch_2021_03/
├── sample_001_chr1.bed/bim/fam
├── sample_001_chr2.bed/bim/fam
├── ...
├── sample_001_chr22.bed/bim/fam
├── sample_002_chr1.bed/bim/fam
├── sample_002_chr2.bed/bim/fam
└── ... (920 samples × 22 chromosomes = 20,240 files)
```

### Sample Sheet Entry:
```csv
platform_id,batch_id,input_path,file_type,build,file_structure
GSAv2,batch_2021_03,/data/gsa_v2/batch_2021_03,plink,hg38,individual_chr_split
```

### What Pipeline Does:
```
1. Finds all files matching pattern: *_chr*.bed
2. Groups by sample ID (extracts from filename)
3. For each sample: Concatenates chr1-22 into single file
4. Processes concatenated samples in parallel
5. Merges samples into batch
```

---

## Scenario 3: Pre-Merged Batch File

**Time-saver if you've already merged** - One file contains all samples.

### File Structure:
```
/data/omni25/
├── cohort_A_merged.bed
├── cohort_A_merged.bim
└── cohort_A_merged.fam
```

### Sample Sheet Entry:
```csv
platform_id,batch_id,input_path,file_type,build,file_structure
Omni25,cohort_A,/data/omni25/cohort_A_merged,plink,hg19,merged_batch
```

### What Pipeline Does:
```
1. Detects merged batch (no sample-level processing needed)
2. Skips sample discovery and merging steps
3. Goes directly to QC → liftover → strand check
```

---

## Scenario 4: Pre-Merged Batch, Split by Chromosome

**After imputation prep elsewhere** - Merged batch already split by chromosome.

### File Structure:
```
/data/mega/phase2/
├── phase2_all_samples_chr1.bed/bim/fam
├── phase2_all_samples_chr2.bed/bim/fam
├── ...
└── phase2_all_samples_chr22.bed/bim/fam
```

### Sample Sheet Entry:
```csv
platform_id,batch_id,input_path,file_type,build,file_structure
MEGA,phase2,/data/mega/phase2/phase2_all_samples,plink,hg38,merged_chr_split
```

### What Pipeline Does:
```
1. Recognizes pre-split chromosomes
2. Processes each chromosome independently
3. Skips sample-level steps, goes to QC
```

---

## Scenario 5: Mixed Structures in One Study

**Real-world complexity** - Different batches organized differently.

### Sample Sheet:
```csv
platform_id,batch_id,input_path,file_type,build,file_structure
GSAv1,early_batch,/data/gsa/early,plink,hg19,individual_samples
GSAv1,mid_batch,/data/gsa/mid,plink,hg19,individual_chr_split
GSAv1,late_batch,/data/gsa/late_merged,plink,hg38,merged_batch
```

### What Pipeline Does:
```
1. Processes each batch according to its file_structure
2. Brings all to same state (merged batch, hg38)
3. Merges all batches into single GSAv1 platform
```

---

## Scenario 6: VCF Format (Individual or Merged)

**VCF instead of PLINK** - Same logic applies.

### Individual VCF Samples:
```csv
platform_id,batch_id,input_path,file_type,build,file_structure
Array_v3,batch1,/vcf/samples,vcf,hg38,individual_samples
```

Directory contains: `sample_001.vcf.gz`, `sample_002.vcf.gz`, etc.

### Pre-Merged VCF:
```csv
platform_id,batch_id,input_path,file_type,build,file_structure
Array_v3,batch2,/vcf/merged/batch2_all,vcf,hg38,merged_batch
```

File: `batch2_all.vcf.gz`

---

## File Naming Pattern Detection

### How Pipeline Identifies Files

For `individual_samples`:
```bash
# PLINK: Finds *.bed, removes .bed to get prefix
find ${input_path} -name "*.bed" | sed 's/.bed$//'

# VCF: Finds *.vcf.gz
find ${input_path} -name "*.vcf.gz"
```

For `individual_chr_split`:
```bash
# Detects pattern: *_chr*.bed or *_chr*.vcf.gz
# Groups by sample: sample_001_chr1, sample_001_chr2 → sample_001
```

For `merged_batch`:
```bash
# Uses input_path as file prefix directly
# Expects: ${input_path}.bed/bim/fam OR ${input_path}.vcf.gz
```

---

## Custom File Patterns (Advanced)

If your files don't match standard patterns, add `file_pattern` column:

```csv
platform_id,batch_id,input_path,file_type,build,file_structure,file_pattern
GSAv1,special,/data/special,plink,hg19,individual_samples,*_QC_filtered_final.bed
```

Pipeline will use: `find ${input_path} -name "${file_pattern}"`

---

## Complete Example Sample Sheets

### Example 1: Simple Study (All Same Structure)
```csv
platform_id,batch_id,input_path,file_type,build,file_structure
Omni25,cohort_A,/study1/omni/cohort_a,plink,hg19,individual_samples
Omni25,cohort_B,/study1/omni/cohort_b,plink,hg19,individual_samples
Omni25,cohort_C,/study1/omni/cohort_c,plink,hg38,individual_samples
```

### Example 2: Multi-Platform Meta-Analysis
```csv
platform_id,batch_id,input_path,file_type,build,file_structure
GSAv1,study_A,/arrays/gsa1/study_a,plink,hg19,individual_samples
GSAv2,study_B,/arrays/gsa2/study_b,plink,hg38,individual_chr_split
MEGA,study_C,/arrays/mega/merged_c,plink,hg38,merged_batch
Omni25,study_D,/arrays/omni/study_d,vcf,hg19,individual_samples
CoreExome,study_E,/exome/study_e/merged,vcf,hg38,merged_chr_split
```

### Example 3: Incremental Additions
```csv
platform_id,batch_id,input_path,file_type,build,file_structure
GSAv3,wave1,/data/wave1/samples,plink,hg38,individual_samples
GSAv3,wave2,/data/wave2/samples,plink,hg38,individual_samples
GSAv3,wave3,/data/wave3/samples,plink,hg38,individual_samples
```

---

## Validation Script

Before running pipeline, validate your sample sheet:

```bash
#!/bin/bash
# validate_samplesheet.sh

while IFS=',' read -r platform batch path type build structure; do
  # Skip header
  [[ "$platform" == "platform_id" ]] && continue
  
  echo "Checking: $platform / $batch"
  
  # Check if path exists
  if [[ "$structure" == "individual_samples" ]] || [[ "$structure" == "individual_chr_split" ]]; then
    # Should be a directory
    if [[ ! -d "$path" ]]; then
      echo "  ✗ ERROR: Directory not found: $path"
      continue
    fi
    
    # Count files
    if [[ "$type" == "plink" ]]; then
      n_files=$(find "$path" -name "*.bed" | wc -l)
    else
      n_files=$(find "$path" -name "*.vcf.gz" | wc -l)
    fi
    
    echo "  ✓ Found $n_files files in $path"
    
  elif [[ "$structure" == "merged_batch" ]]; then
    # Should be a file prefix
    if [[ "$type" == "plink" ]]; then
      if [[ -f "${path}.bed" ]]; then
        echo "  ✓ Found merged batch: ${path}.bed"
      else
        echo "  ✗ ERROR: File not found: ${path}.bed"
      fi
    else
      if [[ -f "${path}.vcf.gz" ]]; then
        echo "  ✓ Found merged VCF: ${path}.vcf.gz"
      else
        echo "  ✗ ERROR: File not found: ${path}.vcf.gz"
      fi
    fi
    
  elif [[ "$structure" == "merged_chr_split" ]]; then
    # Should find chr-split files
    n_chr=$(find "$(dirname "$path")" -name "$(basename "$path")_chr*.bed" | wc -l)
    if [[ $n_chr -gt 0 ]]; then
      echo "  ✓ Found $n_chr chromosome files"
    else
      echo "  ✗ ERROR: No chromosome files found"
    fi
  fi
  
done < sample_sheet.csv

echo ""
echo "Validation complete! Check for any ✗ errors above."
```

Usage:
```bash
chmod +x validate_samplesheet.sh
./validate_samplesheet.sh
```

---

## Module 1 Processing Logic

Here's how Module 1 handles each `file_structure`:

```groovy
workflow PRE_IMPUTATION_QC {
    take:
    sample_sheet_ch
    
    main:
    // Branch by file_structure
    sample_sheet_ch.branch {
        individual: it[5] == 'individual_samples'
        chr_split: it[5] == 'individual_chr_split'
        merged: it[5] == 'merged_batch'
        merged_chr: it[5] == 'merged_chr_split'
    }.set { branched }
    
    // Path 1: Individual samples (all chr together)
    individual_samples = DISCOVER_SAMPLES(branched.individual)
    processed_indiv = PROCESS_INDIVIDUAL_SAMPLE(individual_samples)
    batch1 = MERGE_BATCH_SAMPLES(processed_indiv)
    
    // Path 2: Individual samples (chr split)
    chr_split_samples = DISCOVER_CHR_SPLIT_SAMPLES(branched.chr_split)
    concatenated = CONCATENATE_CHROMOSOMES(chr_split_samples)
    processed_concat = PROCESS_INDIVIDUAL_SAMPLE(concatenated)
    batch2 = MERGE_BATCH_SAMPLES(processed_concat)
    
    // Path 3: Pre-merged batch
    batch3 = VALIDATE_MERGED_BATCH(branched.merged)
    
    // Path 4: Pre-merged, chr split
    batch4 = MERGE_CHR_SPLIT_BATCH(branched.merged_chr)
    
    // Combine all paths
    all_batches = batch1.mix(batch2, batch3, batch4)
    
    // Continue with liftover, strand check, etc.
    // ... rest of pipeline
}
```

---

## Decision Tree: Which file_structure Do I Use?

```
START: Look at your data
    ↓
Q1: Is data already merged at batch level?
    YES → Q2: Are chromosomes split into separate files?
              YES → Use: merged_chr_split
              NO → Use: merged_batch
    NO → Q3: Are samples split by chromosome?
              YES → Use: individual_chr_split
              NO → Use: individual_samples
```

### Quick Reference:

| What You Have | Use This |
|---------------|----------|
| `sample_001.bed` (all chr together) | `individual_samples` |
| `sample_001_chr1.bed`, `sample_001_chr2.bed`, ... | `individual_chr_split` |
| `batch1_all_samples.bed` (one file, all chr) | `merged_batch` |
| `batch1_chr1.bed`, `batch1_chr2.bed`, ... | `merged_chr_split` |

---

## Troubleshooting

### Issue: "No samples found"
**Check**: 
- Is `input_path` correct?
- Does `file_structure` match your data?
- Run validation script

### Issue: "Duplicate sample IDs"
**Solution**: 
- If legitimate duplicates: Set `--allow_duplicate_samples true`
- If error: Check file naming patterns

### Issue: "Mixed chromosome counts"
**Cause**: Some samples missing chromosomes
**Solution**: Pipeline will warn but continue (QC will filter)

### Issue: "Files in both PLINK and VCF format"
**Solution**: 
- Create separate rows in sample sheet
- One row per format type

---

## Advanced: Custom Naming Patterns

If your files have non-standard names, add `file_pattern`:

```csv
platform_id,batch_id,input_path,file_type,build,file_structure,file_pattern
Special,batch1,/data/special,plink,hg19,individual_samples,*_filtered_v2_QC.bed
Special,batch2,/data/special,plink,hg19,individual_chr_split,*_chr*_final.bed
```

Or use regex patterns:
```csv
platform_id,batch_id,input_path,file_type,build,file_structure,file_regex
Complex,batch1,/data,plink,hg19,individual_samples,^[A-Z0-9]+_sample_[0-9]+\.bed$
```

---

## Summary

### Supported Input Types:
✅ Individual samples, all chromosomes together  
✅ Individual samples, split by chromosome  
✅ Pre-merged batches  
✅ Pre-merged batches, split by chromosome  
✅ PLINK format (.bed/.bim/.fam)  
✅ VCF format (.vcf.gz)  
✅ Mixed formats within study  
✅ Mixed genome builds within platform  
✅ Custom file naming patterns  

### Key Points:
1. **Flexibility** - Pipeline adapts to your data structure
2. **Per-batch specification** - Each batch can be different
3. **Automatic detection** - Smart file discovery
4. **Validation** - Check your setup before running
5. **Parallelization** - Optimized for each structure type

---

**Version**: 2.0 - Multi-structure support  
**Updated**: November 6, 2024  
**Status**: Production-ready
