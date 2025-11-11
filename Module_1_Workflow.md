# MODULE 1: PRE-IMPUTATION QC - DETAILED WORKFLOW v5.0
**Last Updated:** November 11, 2025  
**Status:** Reflects Module 1 implementation code based on Module 1 v5.0

---

##  CRITICAL WORKFLOW PRINCIPLES

### QC Timing
- **NO QC** until AFTER: samples merged → batches aligned → batches lifted (if hg19) → batches union merged to platform → service validation
- **FIRST QC APPLICATION:** Process 8 - Light QC only (no HWE, no MAF beyond monomorphic)

### Merge Strategies
- **Sample → Batch:** UNION merge (keeps all variants from all samples)
- **Batch → Platform:** UNION merge (keeps all variants from all batches)
- **Module 4 only:** INTERSECTION merge (keeps only shared variants across platforms)

### Reference Alignment Strategy
- **Per-batch genome build detection** (not global parameter)
- **hg19 batches:** Align to hg19 reference → CrossMap liftover to hg38
- **hg38 batches:** Align to hg38 reference → Skip liftover
- **Both paths converge:** All batches in hg38 before platform merge

---

##  COMPLETE MODULE 1 PROCESS FLOW

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         INPUT: SAMPLE SHEET                              │
│  Columns: platform_id, batch_id, input_path, file_type (plink/vcf),     │
│           build (hg19/hg38), file_structure (4 types)                    │
└─────────────────────────────────────────────────────────────────────────┘
                                    ↓
┌═════════════════════════════════════════════════════════════════════════┐
║ PROCESS 1: DISCOVER_AND_LOAD_SAMPLES                                    ║
║ Container: python:3.11                                                   ║
├─────────────────────────────────────────────────────────────────────────┤
║ Input:  Sample sheet row                                                ║
║ Action: Discover files based on file_structure:                         ║
║         - individual_samples: Find all .bed or .vcf.gz files            ║
║         - individual_chr_split: Find *_chr*.bed or *_chr*.vcf.gz        ║
║         - merged_batch: Single file path                                ║
║         - merged_chr_split: Single file path (chr-split)                ║
║                                                                          ║
║ Output: file_list.txt with discovered file paths                        ║
║ QC:     NONE                                                             ║
║ Parallel: One job per batch                                             ║
└═════════════════════════════════════════════════════════════════════════┘
                                    ↓
                        [flatMap to individual samples]
                                    ↓
┌═════════════════════════════════════════════════════════════════════════┐
║ PROCESS 2: PREPARE_SAMPLES                                              ║
║ Container: plink2 + bcftools (mulled-v2-b0f77c6e0af9c4ab...)           ║
├─────────────────────────────────────────────────────────────────────────┤
║ Input:  Individual sample files (4 possible structures)                 ║
║                                                                          ║
║ BRANCH A: individual_samples or individual_chr_split                    ║
║   IF chr_split:                                                          ║
║     • VCF: bcftools concat chr1-22 → plink2 convert to PLINK            ║
║     • PLINK: plink --merge-list chr1-22 → single PLINK                  ║
║   IF not chr_split:                                                      ║
║     • VCF: plink2 convert to PLINK                                      ║
║     • PLINK: Already correct format                                     ║
║                                                                          ║
║ BRANCH B: merged_batch or merged_chr_split                              ║
║   [Goes to PROCESS 3B instead]                                          ║
║                                                                          ║
║ Output: Sample-level PLINK files (unfiltered)                           ║
║ QC:     NONE (format conversion only)                                   ║
║ Parallel: One job per sample                                            ║
║ Logging: Input/output variant counts                                    ║
└═════════════════════════════════════════════════════════════════════════┘
                                    ↓
                      [Group samples by batch_id]
                                    ↓
        ┌─────────────────────┬─────────────────────┐
        │ Individual samples  │  Pre-merged batches │
        │ (Process 3)         │  (Process 3B)       │
        └─────────────────────┴─────────────────────┘
                ↓                           ↓
┌═════════════════════════════════════════════════════════════════════════┐
║ PROCESS 3: MERGE_SAMPLES_TO_BATCH (UNION)                              ║
║ Container: plink:1.90b6.21                                              ║
├─────────────────────────────────────────────────────────────────────────┤
║ Input:  All samples from same batch (BRANCH A only)                     ║
║                                                                          ║
║ Action:                                                                  ║
║   IF single sample:                                                      ║
║     • Rename to batch prefix                                            ║
║   IF multiple samples:                                                   ║
║     • plink --merge-list with --merge-mode 6 (UNION)                    ║
║     • Keeps ALL variants from ALL samples                               ║
║     • Missing genotypes filled as 0/0                                   ║
║                                                                          ║
║ Output: Batch-level PLINK (unfiltered, union of all variants)           ║
║ QC:     NONE                                                             ║
║ Parallel: One job per batch                                             ║
║ Logging: Variants per sample before merge                               ║
║          Total variants in union                                        ║
║          Retention rate per sample                                      ║
║                                                                          ║
║ WHEN:   Only runs for individual_samples and individual_chr_split       ║
└═════════════════════════════════════════════════════════════════════════┘
                                    ↓
┌═════════════════════════════════════════════════════════════════════════┐
║ PROCESS 3B: PASSTHROUGH_MERGED_BATCH                                    ║
║ Container: python:3.11                                                   ║
├─────────────────────────────────────────────────────────────────────────┤
║ Input:  Pre-merged batch files (BRANCH B)                               ║
║                                                                          ║
║ Action: Simple copy/rename to batch prefix                              ║
║         (batch already merged - no merge needed)                        ║
║                                                                          ║
║ Output: Batch-level PLINK (already merged)                              ║
║ QC:     NONE                                                             ║
║ Parallel: One job per batch                                             ║
║ Logging: Batch already merged confirmation                              ║
║                                                                          ║
║ WHEN:   Only runs for merged_batch and merged_chr_split                 ║
└═════════════════════════════════════════════════════════════════════════┘
                                    ↓
                         [All batches converge]
                                    ↓
                         [Split by genome build]
                                    ↓
              ┌─────────────────┬─────────────────┐
              │   hg19 batches  │  hg38 batches   │
              └─────────────────┴─────────────────┘
                      ↓                   ↓
┌═════════════════════════════════════════════════════════════════════════┐
║ PROCESS 4: ALIGN_TO_REFERENCE (Per Batch, Per Build)                   ║
║ Container: bcftools:1.18                                                 ║
║ CRITICAL: Runs TWICE - once for hg19 batches, once for hg38 batches    ║
├─────────────────────────────────────────────────────────────────────────┤
║ Input:  Batch PLINK + reference FASTA (hg19 or hg38)                   ║
║                                                                          ║
║ Workflow:                                                                ║
║   1. Convert PLINK → VCF (plink2)                                       ║
║   2. Split multiallelics: bcftools norm -m-any                          ║
║   3. Check REF alleles: bcftools norm --check-ref w                     ║
║   4. Fix REF mismatches: bcftools +fixref                               ║
║      • Identifies strand flips                                          ║
║      • Identifies REF/ALT swaps                                         ║
║      • Fixes resolvable mismatches                                      ║
║      • Excludes unresolvable conflicts                                  ║
║   5. Sort: bcftools sort                                                ║
║   6. Index: bcftools index -c                                           ║
║   7. Convert back to PLINK                                              ║
║                                                                          ║
║ Output: Aligned batch PLINK (hg19 or hg38)                              ║
║ QC:     NONE (alignment only)                                            ║
║ Parallel: One job per batch                                             ║
║ Logging: Variants before/after alignment                                ║
║          Strand flips applied                                           ║
║          REF/ALT swaps applied                                          ║
║          Conflicts excluded                                             ║
║          Retention percentage                                           ║
└═════════════════════════════════════════════════════════════════════════┘
                                    ↓
              ┌─────────────────┬─────────────────┐
              │ hg19 aligned    │  hg38 aligned   │
              │ (needs liftover)│  (skip liftover)│
              └─────────────────┴─────────────────┘
                      ↓                   ↓
┌═════════════════════════════════════════════════════════════════════════┐
║ PROCESS 5: LIFTOVER_HG19_TO_HG38                                        ║
║ Container: CrossMap + plink2 + bcftools (mulled-v2-27978155...)        ║
├─────────────────────────────────────────────────────────────────────────┤
║ Input:  hg19-aligned batch PLINK                                        ║
║         hg19ToHg38.over.chain.gz                                        ║
║         hg38 reference FASTA                                            ║
║                                                                          ║
║ Workflow:                                                                ║
║   1. Convert PLINK → VCF with chr prefix                                ║
║   2. Index VCF: tabix -p vcf                                            ║
║   3. CrossMap.py vcf (VCF mode - CRITICAL)                              ║
║      • Updates genomic coordinates hg19 → hg38                          ║
║      • Validates REF alleles against hg38 reference                     ║
║      • Handles strand orientation automatically                         ║
║      • Excludes variants with mismatched REF                            ║
║      • Creates .unmap file for failed lifts                             ║
║   4. Compress and sort lifted VCF                                       ║
║   5. Index: tabix -p vcf                                                ║
║   6. Convert lifted VCF → PLINK                                         ║
║                                                                          ║
║ Output: hg38-lifted batch PLINK                                         ║
║         Lifted VCF + index                                              ║
║         Unlifted variants list (optional)                               ║
║                                                                          ║
║ QC:     NONE (liftover only)                                             ║
║ Parallel: One job per hg19 batch                                        ║
║ Logging: Input variants (hg19)                                          ║
║          Successfully lifted                                            ║
║          Failed to lift                                                 ║
║          Liftover success rate                                          ║
║                                                                          ║
║ WHEN:   Only runs for hg19 batches                                      ║
║ NOTE:   hg38 batches skip this process entirely                         ║
└═════════════════════════════════════════════════════════════════════════┘
                                    ↓
              [hg38-aligned batches skip liftover]
                                    ↓
                    [All batches now in hg38]
                                    ↓
                  [Group batches by platform_id]
                                    ↓
┌═════════════════════════════════════════════════════════════════════════┐
║ PROCESS 6: UNION_MERGE_TO_PLATFORM                                      ║
║ Container: plink:1.90b6.21                                              ║
├─────────────────────────────────────────────────────────────────────────┤
║ Input:  All batches from same platform (all hg38 now)                   ║
║                                                                          ║
║ Action:                                                                  ║
║   IF single batch:                                                       ║
║     • Rename to platform prefix                                         ║
║   IF multiple batches:                                                   ║
║     • plink --merge-list with --merge-mode 6 (UNION)                    ║
║     • Keeps ALL variants from ALL batches                               ║
║     • Preserves population-specific variants                            ║
║     • Missing genotypes filled as 0/0                                   ║
║                                                                          ║
║ Output: Platform-level PLINK (hg38, unfiltered, union)                  ║
║ QC:     NONE (yet)                                                       ║
║ Parallel: One job per platform                                          ║
║ Logging: Variants per batch before merge                                ║
║          Total unique variants in union                                 ║
║          Retention rate per batch                                       ║
║          Overall retention: ~99.999%                                    ║
║                                                                          ║
║ CRITICAL: This is UNION merge (different from Module 4)                 ║
└═════════════════════════════════════════════════════════════════════════┘
                                    ↓
                  [Split into two parallel paths]
                                    ↓
        ┌─────────────────────────────────────────────┐
        │ TOPMed Path          │   AnVIL Path         │
        │ (Process 7A)         │   (Process 7B)       │
        └─────────────────────────────────────────────┘
                ↓                           ↓
┌═════════════════════════════════════════════════════════════════════════┐
║ PROCESS 7A: TOPMED_VALIDATION                                           ║
║ Container: plink + perl + bcftools (mulled-v2-b0f77c6e0af9c4ab...)     ║
├─────────────────────────────────────────────────────────────────────────┤
║ Input:  Platform PLINK                                                   ║
║         Will Rayner strand check script (adapted)                       ║
║         TOPMed Freeze 10 reference                                      ║
║                                                                          ║
║ Action:                                                                  ║
║   1. Calculate allele frequencies: plink --freq                         ║
║   2. Run Will Rayner script with -n and -v flags:                       ║
║      • -n: Keep ALL SNPs (no exclusions by AF differences)              ║
║      • -v: Verbose output for detailed logging                          ║
║      • Identifies strand flips                                          ║
║      • Identifies REF/ALT swaps                                         ║
║      • Flags conflicts (but doesn't exclude with -n)                    ║
║   3. Generate validation report files:                                   ║
║      • *-Strand-Flip.txt (variants needing flip)                        ║
║      • *-Force-Allele1.txt (REF/ALT swaps)                              ║
║      • *-remove.txt (flagged variants - info only with -n)              ║
║                                                                          ║
║ Output: Platform PLINK (unchanged - validation only)                    ║
║         Validation report files                                         ║
║         Service marker: "topmed"                                        ║
║                                                                          ║
║ QC:     NONE (validation/QC check only, no exclusions)                  ║
║ Parallel: One job per platform                                          ║
║ Logging: Number of strand flips noted                                   ║
║          Number of force allele noted                                   ║
║          Number flagged for removal (info only)                         ║
║                                                                          ║
║ NOTE:   With -n flag, these are QC checks ONLY                          ║
║         No variants excluded by AF differences                          ║
└═════════════════════════════════════════════════════════════════════════┘
                                    ↓
┌═════════════════════════════════════════════════════════════════════════┐
║ PROCESS 7B: ANVIL_VALIDATION                                            ║
║ Container: python:3.11                                                   ║
├─────────────────────────────────────────────────────────────────────────┤
║ Input:  Platform PLINK                                                   ║
║                                                                          ║
║ Action: Simple passthrough                                               ║
║         Data already aligned to hg38 in Process 4                       ║
║         No additional validation needed for AnVIL                        ║
║                                                                          ║
║ Output: Platform PLINK (unchanged)                                      ║
║         Service marker: "anvil"                                         ║
║                                                                          ║
║ QC:     NONE                                                             ║
║ Parallel: One job per platform                                          ║
║ Logging: Confirmation that data ready for AnVIL                          ║
║          Variant and sample counts                                      ║
└═════════════════════════════════════════════════════════════════════════┘
                                    ↓
                    [Both paths converge]
                                    ↓
┌═════════════════════════════════════════════════════════════════════════┐
║ PROCESS 8: LIGHT_QC *** FIRST TIME QC IS APPLIED ***                   ║
║ Container: plink:1.90b6.21                                              ║
├─────────────────────────────────────────────────────────────────────────┤
║ Input:  Platform PLINK (tagged with service: topmed or anvil)           ║
║                                                                          ║
║ QC Pipeline (applied in order):                                         ║
║                                                                          ║
║   Step 1: Biallelic SNPs only                                           ║
║      • plink --snps-only just-acgt --max-alleles 2                      ║
║      • Removes indels and non-ACGT variants                             ║
║      • Logging: Variants before/removed/after                           ║
║                                                                          ║
║   Step 2: Remove duplicates                                             ║
║      • plink --rm-dup exclude-all                                       ║
║      • Logging: Variants before/removed/after                           ║
║                                                                          ║
║   Step 3: Remove monomorphic (MAF > 0)                                  ║
║      • plink --maf 0.000001                                             ║
║      • NOT filtering by meaningful MAF threshold                        ║
║      • Logging: Variants before/removed/after                           ║
║                                                                          ║
║   Step 4: Variant call rate ≥95%                                        ║
║      • plink --geno 0.05                                                ║
║      • Logging: Variants before/removed/after                           ║
║                                                                          ║
║   Step 5: Sample call rate ≥95%                                         ║
║      • plink --mind 0.05                                                ║
║      • Logging: Samples before/removed/after                            ║
║                                                                          ║
║ EXPLICITLY NOT APPLIED:                                                 ║
║   ✗ HWE filtering (Module 6 only)                                       ║
║   ✗ MAF filtering beyond monomorphic (Module 6 only)                    ║
║   ✗ Sex checks (Module 6 only)                                          ║
║   ✗ Heterozygosity filtering (Module 6 only)                            ║
║   ✗ Relatedness filtering (Module 6 only)                               ║
║                                                                          ║
║ Output: Platform PLINK (light QC applied, tagged by service)            ║
║ Parallel: One job per platform per service (2 jobs per platform)        ║
║ Logging: Detailed QC log with counts at each step                       ║
║          Variant and sample retention percentages                       ║
└═════════════════════════════════════════════════════════════════════════┘
                                    ↓
                  [Split by service for VCF creation]
                                    ↓
        ┌─────────────────────────────────────────────┐
        │ TOPMed Path          │   AnVIL Path         │
        │ (Process 9A)         │   (Process 9B)       │
        └─────────────────────────────────────────────┘
                ↓                           ↓
┌═════════════════════════════════════════════════════════════════════════┐
║ PROCESS 9A: CREATE_TOPMED_VCFS (Separate chr1-22)                      ║
║ Container: bcftools:1.18                                                 ║
├─────────────────────────────────────────────────────────────────────────┤
║ Input:  Platform PLINK (QC'd) + chromosome number                       ║
║                                                                          ║
║ Action per chromosome (chr1-22):                                         ║
║   1. Extract chromosome: plink2 --chr ${chr}                            ║
║   2. Convert to VCF 4.2: plink2 --export vcf-4.2 bgz                    ║
║   3. Add chr prefix: --output-chr chr26                                 ║
║   4. Set variant IDs: bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT'║
║   5. Clean FORMAT: bcftools annotate -x INFO,^FORMAT/GT                 ║
║   6. Sort: bcftools sort                                                ║
║   7. Index: bcftools index -c (CSI index)                               ║
║                                                                          ║
║ Output: Separate VCF per chromosome                                     ║
║         platform_topmed_chr1.vcf.gz + .csi                              ║
║         platform_topmed_chr2.vcf.gz + .csi                              ║
║         ...                                                             ║
║         platform_topmed_chr22.vcf.gz + .csi                             ║
║                                                                          ║
║ Format: VCFv4.2, bgzipped, CSI indexed                                  ║
║ Parallel: 22 jobs per platform (one per chromosome)                     ║
║ Logging: Variants per chromosome                                        ║
║          Format specifications                                          ║
║                                                                          ║
║ PURPOSE: TOPMed server requires separate chromosome files               ║
║          All chr1-22 submitted together in single job                   ║
└═════════════════════════════════════════════════════════════════════════┘
                                    ↓
┌═════════════════════════════════════════════════════════════════════════┐
║ PROCESS 9B: CREATE_ANVIL_VCF (ALL chr1-22 concatenated)                ║
║ Container: bcftools:1.18                                                 ║
├─────────────────────────────────────────────────────────────────────────┤
║ Input:  Platform PLINK (QC'd)                                           ║
║                                                                          ║
║ Action:                                                                  ║
║   1. Create per-chromosome VCFs (chr1-22):                              ║
║      For each chr:                                                      ║
║        • plink2 --chr ${chr} --export vcf-4.3 bgz                       ║
║        • Add chr prefix: --output-chr chr26                             ║
║        • Set variant IDs: bcftools annotate --set-id                    ║
║        • Clean FORMAT: bcftools annotate -x INFO,^FORMAT/GT             ║
║        • Index: bcftools index -c                                       ║
║                                                                          ║
║   2. Concatenate ALL chromosomes:                                        ║
║      • bcftools concat --file-list (all chr1-22)                        ║
║      • Creates single VCF with all autosomes                            ║
║                                                                          ║
║   3. Sort and index final VCF:                                          ║
║      • bcftools sort                                                    ║
║      • bcftools index -c                                                ║
║                                                                          ║
║ Output: Single concatenated VCF for all autosomes                       ║
║         platform_anvil_all_autosomes.vcf.gz + .csi                      ║
║                                                                          ║
║ Format: VCFv4.3, bgzipped, CSI indexed                                  ║
║ Parallel: One job per platform                                          ║
║ Logging: Total variants across all chromosomes                          ║
║          Chromosomes included (chr1-chr22)                              ║
║          Format specifications                                          ║
║                                                                          ║
║ CRITICAL: All autosomes in ONE file = 95% quota savings!                ║
║           AnVIL charges per genome, not per chromosome                  ║
║           Single submission instead of 22 separate jobs                 ║
└═════════════════════════════════════════════════════════════════════════┘
                                    ↓
┌═════════════════════════════════════════════════════════════════════════┐
║                         MODULE 1 OUTPUTS                                 ║
├─────────────────────────────────────────────────────────────────────────┤
║                                                                          ║
║ FOR EACH PLATFORM:                                                       ║
║                                                                          ║
║ TOPMed-ready files:                                                      ║
║   • platform_topmed_chr1.vcf.gz + .csi                                  ║
║   • platform_topmed_chr2.vcf.gz + .csi                                  ║
║   • ...                                                                 ║
║   • platform_topmed_chr22.vcf.gz + .csi                                 ║
║   Total: 22 VCF files + 22 index files                                  ║
║   Format: VCFv4.2, bgzipped, CSI indexed                                ║
║   Ready for: TOPMed Imputation Server submission                        ║
║              (all chr1-22 submitted as single job)                      ║
║                                                                          ║
║ AnVIL-ready files:                                                       ║
║   • platform_anvil_all_autosomes.vcf.gz + .csi                          ║
║   Total: 1 VCF file + 1 index file                                      ║
║   Format: VCFv4.3, bgzipped, CSI indexed                                ║
║   Ready for: All of Us AnVIL submission                                 ║
║              (single file = optimal quota usage)                        ║
║                                                                          ║
║ QC Logs (per platform):                                                 ║
║   • Preparation logs                                                     ║
║   • Merge logs (sample→batch, batch→platform)                           ║
║   • Alignment logs (hg19 or hg38)                                       ║
║   • Liftover logs (if hg19)                                             ║
║   • Validation logs (TOPMed or AnVIL)                                   ║
║   • Light QC logs with detailed filtering steps                         ║
║   • VCF creation logs                                                   ║
║                                                                          ║
║ NEXT MODULE: Files proceed to Module 2 (Imputation)                     ║
└═════════════════════════════════════════════════════════════════════════┘
```

---

##  KEY MODULE 1 CHARACTERISTICS

### Input Flexibility
 4 file structures supported (individual_samples, individual_chr_split, merged_batch, merged_chr_split)  
 Mixed PLINK/VCF formats  
 Mixed genome builds (hg19/hg38) per batch  
 Automatic format detection and conversion

### QC Philosophy
 NO QC until after ALL merging and alignment  
 Light QC only (no HWE, no MAF beyond monomorphic)  
 Preserves maximum variant diversity for ancestry analysis

### Merge Strategy
 UNION merges (sample→batch, batch→platform)  
 Preserves population-specific variants  
 ~99.999% variant retention

### Reference Handling
 Per-batch genome build detection  
 Proper reference alignment (bcftools norm + fixref)  
 Safe liftover with REF validation (CrossMap VCF mode)

### Service-Specific Outputs
 TOPMed: Separate chr1-22 VCFs (VCF 4.2)  
 AnVIL: Concatenated all-autosome VCF (VCF 4.3)  
 Quota-efficient submission strategy

### Parallelization
 Process-level parallelization throughout  
 Sample-level: PREPARE_SAMPLES  
 Batch-level: Multiple processes  
 Platform-level: Final processes  
 Chromosome-level: TOPMed VCF creation only

---

##  CRITICAL NOTES FOR MODULE DEVELOPMENT

1. **Container Usage:**
   - Python containers: Simple scripting, file discovery, passthrough
   - PLINK containers: All merging operations
   - bcftools containers: Alignment, validation, VCF creation
   - Mulled containers: Multi-tool processes (prepare, liftover, validate)

2. **File Format Transitions:**
   - Stay in PLINK as long as possible (minimize conversions)
   - Convert to VCF only when needed:
     * For reference alignment (bcftools operations)
     * For liftover (CrossMap VCF mode)
     * For final service-specific outputs

3. **Build Handling:**
   - Per-batch detection (not global parameter)
   - Separate alignment processes for hg19 and hg38
   - All batches converge to hg38 before platform merge

4. **Logging Requirements:**
   - Count variants/samples before each step
   - Count removed items
   - Count retained items
   - Calculate retention percentages

5. **Critical Differences from Original Plan:**
   - Uses bcftools +fixref instead of CrossMap BED mode for alignment
   - CrossMap VCF mode for liftover (not BED mode)
   - Will Rayner script with -n flag (validation only, no exclusions)
   - AnVIL path simplified (no additional validation beyond alignment)

---

**Version:** 5.0 - Based on Module 1 Code from November 11, 2025 Updates
**Workflow Visual Last Updated:** November 11, 2025  

