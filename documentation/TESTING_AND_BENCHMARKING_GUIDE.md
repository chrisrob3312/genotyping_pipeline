# Genotyping Pipeline Testing Guide


---


## Module Testing Checklist

### Pre-Testing Setup

- [ ] Install Nextflow (v23.04+)
- [ ] Install Apptainer/Singularity
- [ ] Clone repository
- [ ] Download/build all containers
- [ ] Download reference files:
  - [ ] Reference genome FASTA (hg19 and hg38)
  - [ ] Chain files
  - [ ] BRAVO Freeze 10 data
  - [ ] Genetic maps (from TractorWorkflow)
  - [ ] GRAF-anc data files
- [ ] Create test dataset (use `helper_scripts/create_test_data.sh`)

### Test Data Matrix

Test all 8 input format combinations:

| # | Format | Build | Structure | Status |
|---|--------|-------|-----------|--------|
| 1 | PLINK  | hg19  | individual_samples | [ ] |
| 2 | PLINK  | hg19  | combined | [ ] |
| 3 | PLINK  | hg38  | individual_samples | [ ] |
| 4 | PLINK  | hg38  | combined | [ ] |
| 5 | VCF    | hg19  | individual_samples | [ ] |
| 6 | VCF    | hg19  | combined | [ ] |
| 7 | VCF    | hg38  | individual_samples | [ ] |
| 8 | VCF    | hg38  | combined | [ ] |

---

### Module 0: Container Build

**Purpose**: Build Apptainer/Singularity containers

- [ ] All container definitions exist
- [ ] Containers build successfully
- [ ] Containers contain expected tools:
  - [ ] PLINK 1.9 and 2.0
  - [ ] BCFtools
  - [ ] SAMtools
  - [ ] Eagle2
  - [ ] ADMIXTURE
  - [ ] GRAF-anc
  - [ ] Python with pandas, numpy
  - [ ] R with required packages

**Test Command**:
```bash
nextflow run main.nf -profile test,singularity -entry MODULE0_BUILD_CONTAINERS
```

---

### Module 1: Input Validation & Standardization

**Purpose**: Validate inputs, convert formats, coordinate liftover

**Tests**:
- [ ] **Input validation**
  - [ ] Detects invalid sample sheet
  - [ ] Detects missing input files
  - [ ] Validates genome build specification
  - [ ] Validates file format specification

- [ ] **Format conversion**
  - [ ] VCF → PLINK conversion works
  - [ ] PLINK → VCF conversion works
  - [ ] Preserves all variants
  - [ ] Preserves sample IDs

- [ ] **Coordinate liftover**
  - [ ] hg19 → hg38 liftover works
  - [ ] hg38 → hg19 liftover works (if needed)
  - [ ] Unmapped variants logged
  - [ ] Strand flips handled correctly

- [ ] **Chromosome splitting**
  - [ ] Creates per-chromosome files
  - [ ] Handles chrX/chrY correctly
  - [ ] Outputs proper naming convention

**Test Command**:
```bash
nextflow run main.nf -profile test,singularity -entry MODULE1_VALIDATE \
  --sample_sheet test_data/sample_sheet.csv
```

**Verification**:
```bash
# Check output structure
ls -la results/module1/
# Verify variant counts match
bcftools stats results/module1/standardized/*.vcf.gz
```

---

### Module 2: Imputation Server Preparation

**Purpose**: QC and format data for imputation servers

**Tests**:
- [ ] **Pre-imputation QC**
  - [ ] Rayner/McCarthy HRC check runs
  - [ ] Strand flip correction applied
  - [ ] Allele frequency comparison works
  - [ ] Multi-allelic variants handled

- [ ] **Per-server formatting**
  - [ ] Michigan format correct
  - [ ] TOPMed format correct
  - [ ] NYGC format correct
  - [ ] All of Us format correct

- [ ] **Output validation**
  - [ ] VCF files are valid
  - [ ] Chromosome naming correct for each server
  - [ ] Files ready for upload

**Test Command**:
```bash
nextflow run main.nf -profile test,singularity -entry MODULE2_PREP_IMPUTATION \
  --imputation_servers "michigan,topmed"
```

**Verification**:
```bash
# Validate VCF for Michigan server
bcftools view -h results/module2/michigan/chr22.vcf.gz | head
# Check variant counts
bcftools stats results/module2/*/chr22.vcf.gz
```

---

### Module 3: Imputation Server Submission

**Purpose**: Submit to and download from imputation servers

**Tests**:
- [ ] **Server connectivity**
  - [ ] Michigan API authentication
  - [ ] TOPMed API authentication
  - [ ] NYGC/All of Us (manual process documented)

- [ ] **Job submission**
  - [ ] Files upload successfully
  - [ ] Job parameters set correctly
  - [ ] Job ID captured

- [ ] **Status monitoring**
  - [ ] Polls at correct interval
  - [ ] Handles server errors gracefully
  - [ ] Logs progress

- [ ] **Result download**
  - [ ] Downloads all result files
  - [ ] Decryption works (for encrypted results)
  - [ ] Verifies file integrity

**Test Command**:
```bash
# Dry run first
nextflow run main.nf -profile test,singularity -entry MODULE3_SUBMIT \
  --dry_run true --imputation_servers "michigan"
```

---

### Module 4: Post-Imputation QC

**Purpose**: Quality control of imputed data

**Tests**:
- [ ] **Quality filtering**
  - [ ] Rsq/INFO filtering works
  - [ ] MAF filtering applied
  - [ ] Configurable thresholds

- [ ] **Metrics calculation**
  - [ ] Per-variant Rsq extracted
  - [ ] MAF calculated
  - [ ] Sample-level metrics

- [ ] **MagicalRsq-X integration** (optional)
  - [ ] Model files loaded
  - [ ] Calibrated Rsq calculated
  - [ ] Comparison plots generated

- [ ] **QC reports**
  - [ ] HTML report generated
  - [ ] Variant counts by quality
  - [ ] MAF distribution plots

**Test Command**:
```bash
nextflow run main.nf -profile test,singularity -entry MODULE4_POST_QC \
  --rsq_threshold 0.3 --maf_filter 0.01
```

---

### Module 5: Server Result Comparison & Re-imputation

**Purpose**: Compare results across servers, identify failures, re-impute

**Tests**:
- [ ] **Cross-server comparison**
  - [ ] Variant overlap calculated
  - [ ] Concordance metrics computed
  - [ ] Server-specific variants identified

- [ ] **Region identification**
  - [ ] Low-quality regions detected
  - [ ] Failed samples identified
  - [ ] Ancestry-specific issues flagged

- [ ] **Re-imputation workflow**
  - [ ] Problem regions extracted
  - [ ] Secondary server submission works
  - [ ] Results merged correctly

**Test Command**:
```bash
nextflow run main.nf -profile test,singularity -entry MODULE5_COMPARE \
  --servers_to_compare "michigan,topmed"
```

---

### Module 6: Merge & Final QC

**Purpose**: Merge server results, final quality control

**Tests**:
- [ ] **Data merging**
  - [ ] Best-server-per-region logic works
  - [ ] Consensus calling correct
  - [ ] No duplicate variants

- [ ] **Final QC**
  - [ ] HWE filtering
  - [ ] Missingness filtering
  - [ ] Population-specific QC

- [ ] **Output formats**
  - [ ] PLINK format correct
  - [ ] VCF format correct
  - [ ] Dosage preserved

**Test Command**:
```bash
nextflow run main.nf -profile test,singularity -entry MODULE6_MERGE
```

---

### Module 7: Ancestry Estimation

**Purpose**: Global and local ancestry inference

**Tests**:
- [ ] **GRAF-anc**
  - [ ] Runs successfully
  - [ ] Ancestry proportions calculated
  - [ ] Classification output correct

- [ ] **ADMIXTURE**
  - [ ] Multiple K values tested
  - [ ] Cross-validation works
  - [ ] Optimal K selected

- [ ] **Local Ancestry (if enabled)**
  - [ ] **RFMix v2** runs correctly
  - [ ] **FLARE** runs correctly (if selected)
  - [ ] **G-NOMIX** runs correctly (if selected)
  - [ ] **RFMix v1** runs correctly (if selected)
  - [ ] Only selected method(s) run

- [ ] **Visualization**
  - [ ] PCA plots generated
  - [ ] ADMIXTURE bar plots
  - [ ] LAI karyograms (if LAI enabled)

**Test Command**:
```bash
nextflow run main.nf -profile test,singularity -entry MODULE7_ANCESTRY \
  --run_lai true --lai_methods "rfmix2"
```

---

### Integration Testing

**Full Pipeline Run**:
```bash
# Small test dataset
nextflow run main.nf -profile test,singularity \
  --sample_sheet test_data/sample_sheet.csv \
  --outdir results/integration_test \
  -resume

# Check all outputs
ls -laR results/integration_test/
```

**End-to-End Validation**:
- [ ] Pipeline completes without errors
- [ ] All expected output files present
- [ ] Logs capture all steps
- [ ] Resource usage within limits
- [ ] Results reproducible with `-resume`

---

### Troubleshooting Checklist

**Common Issues**:

| Symptom | Likely Cause | Solution |
|---------|--------------|----------|
| Container not found | Missing .sif file | Build with Module 0 or download |
| Liftover fails | Missing chain file | Download from UCSC |
| Rayner check fails | Wrong reference format | Use BRAVO→Rayner converter |
| GRAF-anc error | Missing data files | Download AncSnpPopAFs.txt |
| LAI not running | Wrong `lai_methods` format | Use comma-separated, no spaces |
| Memory error | Process limit too low | Increase in nextflow.config |
| API timeout | Server busy | Increase `monitor_interval` |

**Debug Commands**:
```bash
# View detailed error
nextflow log <run_name> -f status,hash,name,exit,error

# Check specific task
cat work/<hash>/.command.log
cat work/<hash>/.command.err

# Resume from failure
nextflow run main.nf -resume
```

---


## Appendix A: Quick Reference Commands

```bash
# Full pipeline test
nextflow run main.nf -profile test,singularity --outdir test_results

# Single module test
nextflow run main.nf -profile test,singularity -entry MODULE1_VALIDATE

# Resume failed run
nextflow run main.nf -resume

# Generate execution report
nextflow run main.nf -with-report report.html -with-trace trace.txt

# Clean work directory
nextflow clean -f

# View run history
nextflow log
```

## Appendix B: File Checklist

```
Required before testing:
├── config/nextflow.config (parameters updated)
├── resources/
│   ├── containers/*.sif (built)
│   ├── references/
│   │   ├── hg38.fa
│   │   └── hg19ToHg38.over.chain.gz
│   ├── ancestry_references/
│   │   ├── lai_reference/ (populated)
│   │   └── grafanc_data/AncSnpPopAFs.txt
│   └── genetic_maps/ (populated)
└── test_data/
    ├── sample_sheet.csv
    └── input_files/
```
