# Sample Sheet Guide

## Overview
The sample sheet is a CSV file that tells the pipeline where to find your genotyping data and how to process it. Each row represents one "batch" of samples from a particular genotyping platform.

---

## Column Descriptions

### 1. `platform_id` (REQUIRED)
**What it is:** A unique identifier for this batch of samples

**Rules:**
- Must be unique for each row
- No spaces (use underscores instead)
- Should be descriptive
- Alphanumeric characters only

**Examples:**
- `Platform1_Affy6`
- `GSAv1_batch1`
- `Illumina650K_cohort1`

**Purpose:** This ID will be used to name output files, so make it meaningful!

---

### 2. `file_path` (REQUIRED)
**What it is:** Path to your genotype files (WITHOUT file extension)

**For PLINK files:**
```
If you have:
  /data/genotypes/platform1.bed
  /data/genotypes/platform1.bim
  /data/genotypes/platform1.fam

Use: /data/genotypes/platform1
```

**For VCF files (chromosome-split):**
```
If you have:
  /data/vcfs/platform2_chr1.vcf.gz
  /data/vcfs/platform2_chr2.vcf.gz
  ...
  /data/vcfs/platform2_chr22.vcf.gz

Use: /data/vcfs/platform2_chr

The pipeline will automatically find chr1-22
```

**For single VCF (all chromosomes):**
```
If you have:
  /data/vcfs/platform3_all.vcf.gz

Use: /data/vcfs/platform3_all
```

**Path Types:**
- Absolute paths (recommended): `/full/path/to/file`
- Relative paths (from where you run pipeline): `data/file`

---

### 3. `file_type` (REQUIRED)
**What it is:** Format of your genotype files

**Options:**
- `plink` - PLINK binary format (.bed/.bim/.fam)
- `vcf` - VCF format (.vcf or .vcf.gz)
- `vcf_split` - VCF files split by chromosome

**How the pipeline handles each:**
- **plink**: Reads .bed/.bim/.fam directly
- **vcf**: Converts to PLINK if needed, splits by chromosome
- **vcf_split**: Processes each chromosome separately

---

### 4. `build` (REQUIRED)
**What it is:** Genome build/assembly version

**Options:**
- `hg19` (also known as GRCh37)
- `hg38` (also known as GRCh38)

**Why it matters:**
- TOPMed imputation server uses hg38
- The pipeline will perform liftOver if needed
- Build mismatches can cause major issues!

**How to check your build:**
```bash
# Look at chromosome positions in .bim file
head -n 5 yourfile.bim

# hg19 example (chr1 SNP positions):
# rs12345  1  0  752566  A  G

# hg38 example (chr1 SNP positions, slightly different):
# rs12345  1  0  817186  A  G

# Or check your array manufacturer documentation
```

---

### 5. `collection_batch` (OPTIONAL)
**What it is:** Groups samples collected at the same time

**Purpose:**
- Helps track batch effects
- Useful for QC reporting
- Can be used for meta-data

**Examples:**
- Date-based: `2020_01`, `2021_Q1`
- Sequential: `batch1`, `batch2`, `batch3`
- Study-based: `pilot`, `main_cohort`
- Lab-based: `lab_A`, `lab_B`

---

## Complete Examples

### Example 1: Simple Setup (2 platforms, all samples together)
```csv
platform_id,file_path,file_type,build,collection_batch
Affy6_cohort1,/data/geno/affy6_all,plink,hg19,cohort1
Illumina_cohort2,/data/geno/illumina_all,plink,hg38,cohort2
```

### Example 2: Multiple Batches of Same Platform
```csv
platform_id,file_path,file_type,build,collection_batch
GSAv1_jan2020,/data/GSA/2020_01/samples,plink,hg19,2020_01
GSAv1_jun2020,/data/GSA/2020_06/samples,plink,hg19,2020_06
GSAv1_jan2021,/data/GSA/2021_01/samples,plink,hg19,2021_01
GSAv2_jun2021,/data/GSA/2021_06/samples,plink,hg38,2021_06
GSAv2_dec2021,/data/GSA/2021_12/samples,plink,hg38,2021_12
```

### Example 3: Mixed File Types and Builds
```csv
platform_id,file_path,file_type,build,collection_batch
Platform1,/project/arrays/platform1,plink,hg19,batch1
Platform2,/project/arrays/platform2,plink,hg38,batch1
Platform3_VCF,/project/vcfs/platform3_all,vcf,hg38,batch2
Platform4_splitVCF,/project/vcfs/platform4_chr,vcf_split,hg19,batch2
```

---

## How the Pipeline Uses This Information

### Step 1: Input Reading
The pipeline reads this CSV and creates a "channel" (Nextflow concept) with one element per row.

### Step 2: Parallel Processing
Each row is processed **in parallel** through Modules 1-3:
```
Row 1 (Platform1) → Module 1 → Module 2 → Module 3
Row 2 (Platform2) → Module 1 → Module 2 → Module 3  } Parallel!
Row 3 (GSAv1_b1)  → Module 1 → Module 2 → Module 3
```

### Step 3: Grouping for Merge
In Module 4, all platforms are merged together.

### Step 4: Naming Convention
Output files are named using `platform_id`:
```
Platform1_qc.bed
Platform1_topmed_chr1.vcf.gz
Platform1_allofus_filtered.bed
```

---

## Common Questions

### Q: Can I have the same samples across multiple rows?
**A:** No! Each sample should appear in only ONE row. If you have:
- Same samples on different platforms → keep them separate, merge happens in Module 4
- Duplicate samples (for QC) → handle separately before the pipeline

### Q: What if I have samples split across multiple PLINK files?
**A:** First merge them with PLINK, then add to sample sheet:
```bash
# Merge multiple PLINK files
plink --bfile file1 --merge-list merge_list.txt --make-bed --out combined

# Then use 'combined' in your sample sheet
```

### Q: Can I process just some platforms?
**A:** Yes! Just remove rows from the sample sheet or create a new CSV with only the platforms you want.

### Q: Do I need collection_batch?
**A:** No, it's optional. But it's helpful for:
- Tracking where samples came from
- Understanding batch effects
- Organizing your results

### Q: What if my VCF files are split differently (e.g., regions not chromosomes)?
**A:** You'll need to restructure them first:
```bash
# Example: split by chromosome
for chr in {1..22}; do
  bcftools view -r ${chr} input.vcf.gz -Oz -o output_chr${chr}.vcf.gz
done
```

### Q: Can I mix hg19 and hg38 in the same analysis?
**A:** Yes! The pipeline handles this automatically with liftOver, but it's possible some SNPs may be lost or undergo position changes during liftOver via the imputation servers. Performing liftOver before imputation may help accuracy but isn't essential as both All of Us and TOPMed Servers incoporate liftOver as part of their process.

---

## Validation

### Before Running Pipeline
Check your sample sheet for:

1. **No duplicate platform_ids**
```bash
# Should return "OK" if no duplicates
cut -d',' -f1 sample_sheet.csv | sort | uniq -d | wc -l
# Output: 0 (means no duplicates)
```

2. **All files exist**
```bash
# For PLINK files, check .bed files exist:
while IFS=',' read -r id path type build batch; do
  if [[ "$type" == "plink" ]]; then
    if [[ ! -f "${path}.bed" ]]; then
      echo "ERROR: ${path}.bed not found!"
    fi
  fi
done < sample_sheet.csv
```

3. **Correct CSV format**
- No extra commas
- No spaces (except in file paths if quoted)
- Header row present
- Consistent number of columns

---

## Troubleshooting

### "File not found" errors
- Check paths are correct (try absolute paths)
- Check file extensions match file_type
- Check permissions (can you read the files?)

### "Duplicate platform_id" errors
- Make each platform_id unique
- Add batch number or date to differentiate

### "Invalid build" errors
- Use only 'hg19' or 'hg38'
- Check for typos (hg 19, GRCh37, etc. won't work)



---

## Template for Your Use

Save this as `my_sample_sheet.csv`:
```csv
platform_id,file_path,file_type,build,collection_batch
CHANGEME1,/path/to/your/files,plink,hg19,batch1
CHANGEME2,/path/to/your/files,plink,hg38,batch2
```

Then edit with your actual information!

---


```bash
nextflow run main.nf --input sample_sheet.csv --outdir results/
```
