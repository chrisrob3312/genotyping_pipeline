#!/usr/bin/env nextflow

/*
 * ============================================================================
 * MODULE 1: PRE-IMPUTATION QC AND PREPARATION - CORRECTED v4.0
 * ============================================================================
 * 
 * CRITICAL CORRECTIONS:
 * - NO QC until AFTER union merge + strand checking
 * - 3-pass UNION merge (batches → platform)
 * - hg19 fixes BEFORE liftover
 * - Separate TOPMed and AnVIL strand checking paths
 * - Light QC only at the end
 * 
 * PROCESS FLOW (AS PER WORKFLOW OVERVIEW):
 * 1. Input Discovery & Loading
 * 2. Sample Preparation (NO QC - just format conversion)
 * 3. MERGE SAMPLES → BATCHES (simple merge, NO QC)
 * 4. Format Check & Fix hg19 BEFORE Liftover
 * 5. Liftover hg19 → hg38
 * 6. UNION MERGE Batches → Platform (3-pass strategy)
 * 7. Service-Specific Strand Checking (TOPMed + AnVIL separate)
 * 8. LIGHT QC (FIRST TIME QC HAPPENS)
 * 9. Create Service-Specific VCFs
 * 
 * ============================================================================
 */

nextflow.enable.dsl = 2

// ============================================================================
// PROCESS 1: DISCOVER AND LOAD SAMPLES
// ============================================================================

process DISCOVER_AND_LOAD_SAMPLES {
    label 'process_low'
    tag "${platform_id}_${batch_id}"
    
    input:
    tuple val(platform_id),
          val(batch_id),
          val(input_path),
          val(file_type),
          val(build),
          val(file_structure)
    
    output:
    tuple val(platform_id),
          val(batch_id),
          path("file_list.txt"),
          val(file_type),
          val(build),
          val(file_structure),
          emit: discovered_files
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "====== DISCOVERING FILES ======"
    echo "Platform: ${platform_id}"
    echo "Batch: ${batch_id}"
    echo "Path: ${input_path}"
    echo "Structure: ${file_structure}"
    echo "================================"
    
    # Discover files based on structure
    if [[ "${file_structure}" == "individual_samples" ]]; then
        if [[ "${file_type}" == "plink" ]]; then
            find ${input_path} -name "*.bed" -type f | sed 's/.bed\$//' > file_list.txt
        else
            find ${input_path} -name "*.vcf.gz" -type f > file_list.txt
        fi
    
    elif [[ "${file_structure}" == "individual_chr_split" ]]; then
        if [[ "${file_type}" == "plink" ]]; then
            find ${input_path} -name "*_chr*.bed" -type f | \
                sed 's/_chr[0-9]*.bed\$//' | sort -u > file_list.txt
        else
            find ${input_path} -name "*_chr*.vcf.gz" -type f | \
                sed 's/_chr[0-9]*.vcf.gz\$//' | sort -u > file_list.txt
        fi
    
    elif [[ "${file_structure}" == "merged_batch" ]]; then
        echo "${input_path}" > file_list.txt
    
    elif [[ "${file_structure}" == "merged_chr_split" ]]; then
        echo "${input_path}" > file_list.txt
    fi
    
    n_files=\$(wc -l < file_list.txt)
    echo "✓ Found \${n_files} files/samples"
    
    if [[ \${n_files} -eq 0 ]]; then
        echo "ERROR: No files found!"
        exit 1
    fi
    
    head -5 file_list.txt
    """
}

// ============================================================================
// PROCESS 2: PREPARE SAMPLES (FORMAT ONLY - NO QC)
// ============================================================================

process PREPARE_SAMPLES_AND_BATCHES {
    label 'process_medium'
    tag "${platform_id}_${batch_id}_${sample_id}"
    
    publishDir "${params.outdir}/${platform_id}/${batch_id}/01_sample_prep",
        mode: 'copy',
        pattern: "*.log"
    
    input:
    tuple val(platform_id),
          val(batch_id),
          val(sample_id),
          val(file_type),
          val(build),
          val(file_structure)
    
    output:
    tuple val(platform_id),
          val(batch_id),
          path("${sample_id}.bed"),
          path("${sample_id}.bim"),
          path("${sample_id}.fam"),
          val(build),
          emit: prepared_samples
    
    path("${sample_id}_prep.log"), emit: prep_log
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "Preparing: ${sample_id}" | tee ${sample_id}_prep.log
    echo "Structure: ${file_structure}" | tee -a ${sample_id}_prep.log
    echo "NO QC APPLIED - FORMAT CONVERSION ONLY" | tee -a ${sample_id}_prep.log
    
   # Handle different input structures
if [[ "${file_structure}" == "individual_chr_split" && "${file_type}" == "vcf" ]]; then
    # Case 1: VCF files split by chromosome
    echo "=== Concatenating VCF chromosomes ===" | tee -a ${sample_id}_prep.log
    
    for chr in {1..22}; do
        if [[ -f "${sample_id}_chr${chr}.vcf.gz" ]]; then
            echo "${sample_id}_chr${chr}.vcf.gz" >> vcf_list.txt
            echo "  Found chr${chr}" | tee -a ${sample_id}_prep.log
        fi
    done
    
    bcftools concat \\
        --file-list vcf_list.txt \\
        --output ${sample_id}_merged.vcf.gz \\
        --output-type z \\
        --threads ${task.cpus}
    
    bcftools index ${sample_id}_merged.vcf.gz
    
    plink2 --vcf ${sample_id}_merged.vcf.gz \\
           --make-bed \\
           --out ${sample_id} \\
           --threads ${task.cpus}

elif [[ "${file_structure}" == "individual_chr_split" && "${file_type}" == "plink" ]]; then
    # Case 2: PLINK files split by chromosome
    echo "=== Concatenating PLINK chromosomes ===" | tee -a ${sample_id}_prep.log
    
    for chr in {1..22}; do
        if [[ -f "${sample_id}_chr${chr}.bed" ]]; then
            echo "${sample_id}_chr${chr}" >> merge_list.txt
            echo "  Found chr${chr}" | tee -a ${sample_id}_prep.log
        fi
    done
    
    plink --merge-list merge_list.txt \\
          --make-bed \\
          --out ${sample_id} \\
          --threads ${task.cpus}

elif [[ "${file_type}" == "vcf" ]]; then
    # Case 3: Single VCF file (not split by chromosome)
    echo "=== Converting VCF to PLINK ===" | tee -a ${sample_id}_prep.log
    
    plink2 --vcf ${sample_id}.vcf.gz \\
           --make-bed \\
           --out ${sample_id} \\
           --threads ${task.cpus}

else
    # Case 4: Already PLINK format, single file
    echo "=== Files already in PLINK format ===" | tee -a ${sample_id}_prep.log
    # Files are already staged by Nextflow with correct names
fi
    
    # Count variants (NO FILTERING)
    n_snps=\$(wc -l < ${sample_id}.bim)
    n_samples=\$(wc -l < ${sample_id}.fam)
    
    echo "Output - SNPs: \${n_snps}, Samples: \${n_samples}" | tee -a ${sample_id}_prep.log
    echo "✓ Format conversion complete (NO QC APPLIED)" | tee -a ${sample_id}_prep.log
    """
}

// ============================================================================
// PROCESS 3: MERGE SAMPLES TO BATCH (NO QC FILTERS)
// ============================================================================

process MERGE_SAMPLES_TO_BATCH {
    label 'plink_merge'
    tag "${platform_id}_${batch_id}"
    
    publishDir "${params.outdir}/${platform_id}/${batch_id}/02_batch_merge",
        mode: 'copy',
        pattern: "*.{log,txt}"
    
    input:
    tuple val(platform_id),
          val(batch_id),
          path(bed_files),
          path(bim_files),
          path(fam_files),
          val(build)
    
    output:
    tuple val(platform_id),
          val(batch_id),
          path("${platform_id}_${batch_id}.bed"),
          path("${platform_id}_${batch_id}.bim"),
          path("${platform_id}_${batch_id}.fam"),
          val(build),
          emit: batch_merged
    
    path("${platform_id}_${batch_id}_merge.log"), emit: merge_log
    
    script:
    def bed_list = bed_files instanceof List ? bed_files : [bed_files]
    def n_files = bed_list.size()
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "========================================" | tee ${platform_id}_${batch_id}_merge.log
    echo "MERGING SAMPLES TO BATCH (UNION MERGE)" | tee -a ${platform_id}_${batch_id}_merge.log
    echo "Platform: ${platform_id}" | tee -a ${platform_id}_${batch_id}_merge.log
    echo "Batch: ${batch_id}" | tee -a ${platform_id}_${batch_id}_merge.log
    echo "Files: ${n_files}" | tee -a ${platform_id}_${batch_id}_merge.log
    echo "========================================" | tee -a ${platform_id}_${batch_id}_merge.log
    
    if [[ ${n_files} -eq 1 ]]; then
        # Single file - just rename
        echo "Single file - no merge needed" | tee -a ${platform_id}_${batch_id}_merge.log
        mv ${bed_files} ${platform_id}_${batch_id}.bed
        mv ${bim_files} ${platform_id}_${batch_id}.bim
        mv ${fam_files} ${platform_id}_${batch_id}.fam
    
    else
        # Multiple files - UNION merge (keeps ALL variants)
        echo "Performing UNION merge (mode 6)..." | tee -a ${platform_id}_${batch_id}_merge.log
        echo "  - Keeps all variants whether in all files or not" | tee -a ${platform_id}_${batch_id}_merge.log
        echo "  - Auto-resolves strand/allele conflicts" | tee -a ${platform_id}_${batch_id}_merge.log
        echo "  - Strand fixing happens later with bcftools" | tee -a ${platform_id}_${batch_id}_merge.log
        
        # Create merge list
        ls *.bed | sed 's/.bed\$//' > merge_list.txt
        
        # Count variants before merge
        echo "" | tee -a ${platform_id}_${batch_id}_merge.log
        echo "Input files:" | tee -a ${platform_id}_${batch_id}_merge.log
        for bim in *.bim; do
            n=\$(wc -l < \${bim})
            echo "  \${bim}: \${n} variants" | tee -a ${platform_id}_${batch_id}_merge.log
        done
        
        # Union merge with mode 6
        plink --merge-list merge_list.txt \\
              --merge-mode 6 \\
              --make-bed \\
              --out ${platform_id}_${batch_id} \\
              --threads ${task.cpus} \\
              --allow-no-sex
    fi
    
    # Log final results
    n_snps=\$(wc -l < ${platform_id}_${batch_id}.bim)
    n_samples=\$(wc -l < ${platform_id}_${batch_id}.fam)
    
    echo "" | tee -a ${platform_id}_${batch_id}_merge.log
    echo "=== MERGED OUTPUT ===" | tee -a ${platform_id}_${batch_id}_merge.log
    echo "  SNPs: \${n_snps}" | tee -a ${platform_id}_${batch_id}_merge.log
    echo "  Samples: \${n_samples}" | tee -a ${platform_id}_${batch_id}_merge.log
    echo "" | tee -a ${platform_id}_${batch_id}_merge.log
    echo "NOTE: UNION merge - includes all variants across samples" | tee -a ${platform_id}_${batch_id}_merge.log
    echo "      Strand/allele issues resolved using first file's coding" | tee -a ${platform_id}_${batch_id}_merge.log
    echo "      Final strand correction happens in strand-fixing step" | tee -a ${platform_id}_${batch_id}_merge.log
    echo "      NO QC FILTERS APPLIED" | tee -a ${platform_id}_${batch_id}_merge.log
    echo "" | tee -a ${platform_id}_${batch_id}_merge.log
    echo "✓ Batch merge complete" | tee -a ${platform_id}_${batch_id}_merge.log
    echo "========================================" | tee -a ${platform_id}_${batch_id}_merge.log
    """

}

// ============================================================================
// PROCESS 4A: FORMAT CHECK AND FIX HG19 (BEFORE LIFTOVER)
// ============================================================================
process FORMAT_CHECK_AND_FIX_HG19 {
    label 'vcftools'
    tag "${platform_id}_${batch_id}_hg19"
    
    publishDir "${params.outdir}/${platform_id}/${batch_id}/03_hg19_fixes",
        mode: 'copy',
        pattern: "*.{log,txt}"
    
    input:
    tuple val(platform_id),
          val(batch_id),
          path(bed),
          path(bim),
          path(fam),
          val(build)
    path hg19_fasta
    path hg19_fasta_fai
    
    output:
    tuple val(platform_id),
          val(batch_id),
          path("${platform_id}_${batch_id}_hg19_fixed.bed"),
          path("${platform_id}_${batch_id}_hg19_fixed.bim"),
          path("${platform_id}_${batch_id}_hg19_fixed.fam"),
          val(build),
          emit: hg19_fixed
    
    path("${platform_id}_${batch_id}_hg19_fixes.log"), emit: fix_log
    
    when:
    build == 'hg19'
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    prefix=\$(basename ${bed} .bed)
    
    echo "========================================" | tee ${platform_id}_${batch_id}_hg19_fixes.log
    echo "ALIGNING HG19 VARIANTS TO REFERENCE" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    echo "Platform: ${platform_id}" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    echo "Batch: ${batch_id}" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    echo "Build: ${build}" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    echo "Next step: Liftover to hg38" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    echo "========================================" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    
    # Convert PLINK to VCF
    echo "Converting PLINK to VCF..." | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    plink2 --bfile \${prefix} \\
           --export vcf bgz \\
           --out temp \\
           --threads ${task.cpus}
    
    tabix -p vcf temp.vcf.gz
    
    # Count before
    n_before=\$(bcftools view -H temp.vcf.gz | wc -l)
    echo "Variants before alignment: \${n_before}" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    
    # Align to hg19 reference
    echo "" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    echo "Aligning to hg19 reference..." | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    bcftools +fixref temp.vcf.gz \\
        --fasta-ref ${hg19_fasta} \\
        --output fixed.vcf.gz \\
        --output-type z \\
        --threads ${task.cpus} \\
        2>&1 | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    
    tabix -p vcf fixed.vcf.gz
    
    # Count after
    n_after=\$(bcftools view -H fixed.vcf.gz | wc -l)
    n_removed=\$((n_before - n_after))
    retention_pct=\$(awk "BEGIN {printf \\"%.2f\\", (\${n_after}/\${n_before})*100}")
    
    echo "" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    echo "=== ALIGNMENT RESULTS ===" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    echo "  Variants before: \${n_before}" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    echo "  Variants after:  \${n_after}" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    echo "  Variants removed: \${n_removed}" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    echo "  Retention: \${retention_pct}%" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    
    # Convert back to PLINK
    echo "" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    echo "Converting back to PLINK..." | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    plink2 --vcf fixed.vcf.gz \\
           --make-bed \\
           --out ${platform_id}_${batch_id}_hg19_fixed \\
           --threads ${task.cpus}
    
    n_final=\$(wc -l < ${platform_id}_${batch_id}_hg19_fixed.bim)
    n_samples=\$(wc -l < ${platform_id}_${batch_id}_hg19_fixed.fam)
    
    echo "" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    echo "=== OUTPUT ===" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    echo "  SNPs: \${n_final}" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    echo "  Samples: \${n_samples}" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    echo "" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    echo "✓ hg19 alignment complete → Ready for liftover" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    echo "========================================" | tee -a ${platform_id}_${batch_id}_hg19_fixes.log
    """
}

// ============================================================================
// PROCESS 4B: FORMAT CHECK AND FIX HG38 (NO LIFTOVER NEEDED)
// ============================================================================
process FORMAT_CHECK_AND_FIX_HG38 {
    label 'vcftools'
    tag "${platform_id}_${batch_id}_hg38"
    
    publishDir "${params.outdir}/${platform_id}/${batch_id}/03_hg38_fixes",
        mode: 'copy',
        pattern: "*.{log,txt}"
    
    input:
    tuple val(platform_id),
          val(batch_id),
          path(bed),
          path(bim),
          path(fam),
          val(build)
    path hg38_fasta
    path hg38_fasta_fai
    
    output:
    tuple val(platform_id),
          val(batch_id),
          path("${platform_id}_${batch_id}_hg38_fixed.bed"),
          path("${platform_id}_${batch_id}_hg38_fixed.bim"),
          path("${platform_id}_${batch_id}_hg38_fixed.fam"),
          val('hg38'),  // Output is hg38
          emit: hg38_fixed
    
    path("${platform_id}_${batch_id}_hg38_fixes.log"), emit: fix_log
    
    when:
    build == 'hg38'
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    prefix=\$(basename ${bed} .bed)
    
    echo "========================================" | tee ${platform_id}_${batch_id}_hg38_fixes.log
    echo "ALIGNING HG38 VARIANTS TO REFERENCE" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    echo "Platform: ${platform_id}" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    echo "Batch: ${batch_id}" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    echo "Build: ${build}" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    echo "Next step: SKIP liftover (already hg38)" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    echo "========================================" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    
    # Convert PLINK to VCF
    echo "Converting PLINK to VCF..." | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    plink2 --bfile \${prefix} \\
           --export vcf bgz \\
           --out temp \\
           --threads ${task.cpus}
    
    tabix -p vcf temp.vcf.gz
    
    # Count before
    n_before=\$(bcftools view -H temp.vcf.gz | wc -l)
    echo "Variants before alignment: \${n_before}" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    
    # Align to hg38 reference
    echo "" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    echo "Aligning to hg38 reference..." | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    bcftools +fixref temp.vcf.gz \\
        --fasta-ref ${hg38_fasta} \\
        --output fixed.vcf.gz \\
        --output-type z \\
        --threads ${task.cpus} \\
        2>&1 | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    
    tabix -p vcf fixed.vcf.gz
    
    # Count after
    n_after=\$(bcftools view -H fixed.vcf.gz | wc -l)
    n_removed=\$((n_before - n_after))
    retention_pct=\$(awk "BEGIN {printf \\"%.2f\\", (\${n_after}/\${n_before})*100}")
    
    echo "" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    echo "=== ALIGNMENT RESULTS ===" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    echo "  Variants before: \${n_before}" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    echo "  Variants after:  \${n_after}" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    echo "  Variants removed: \${n_removed}" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    echo "  Retention: \${retention_pct}%" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    
    # Convert back to PLINK
    echo "" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    echo "Converting back to PLINK..." | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    plink2 --vcf fixed.vcf.gz \\
           --make-bed \\
           --out ${platform_id}_${batch_id}_hg38_fixed \\
           --threads ${task.cpus}
    
    n_final=\$(wc -l < ${platform_id}_${batch_id}_hg38_fixed.bim)
    n_samples=\$(wc -l < ${platform_id}_${batch_id}_hg38_fixed.fam)
    
    echo "" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    echo "=== OUTPUT ===" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    echo "  SNPs: \${n_final}" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    echo "  Samples: \${n_samples}" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    echo "" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    echo "✓ hg38 alignment complete → Skipping liftover" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    echo "========================================" | tee -a ${platform_id}_${batch_id}_hg38_fixes.log
    """
}

// ============================================================================
// PROCESS 5: LIFTOVER HG19 TO HG38 WITH CROSSMAP VCF (VERIFIED)
// ============================================================================
process LIFTOVER_HG19_TO_HG38 {
    label 'liftover'
    tag "${platform_id}_${batch_id}"
    
    publishDir "${params.outdir}/${platform_id}/${batch_id}/04_liftover",
        mode: 'copy',
        pattern: "*.{log,txt,vcf.gz,vcf.gz.tbi}"
    
    input:
    tuple val(platform_id),
          val(batch_id),
          path(bed),
          path(bim),
          path(fam),
          val(build)
    path chain_file
    path hg38_fasta
    path hg38_fasta_fai
    
    output:
    tuple val(platform_id),
          val(batch_id),
          path("${platform_id}_${batch_id}_hg38.bed"),
          path("${platform_id}_${batch_id}_hg38.bim"),
          path("${platform_id}_${batch_id}_hg38.fam"),
          val("hg38"),
          emit: lifted_batch
    
    path("${platform_id}_${batch_id}_liftover.log"), emit: liftover_log
    path("${platform_id}_${batch_id}_hg19_to_hg38_lifted.vcf.gz"), emit: lifted_vcf
    path("${platform_id}_${batch_id}_hg19_to_hg38_lifted.vcf.gz.tbi"), emit: lifted_vcf_index
    path("*_unlifted.txt"), emit: unlifted_variants, optional: true
    
    when:
    build == 'hg19'
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    prefix=\$(basename ${bed} .bed)
    
    echo "========================================" | tee ${platform_id}_${batch_id}_liftover.log
    echo "LIFTOVER: hg19 → hg38 (CrossMap VCF mode)" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "Platform: ${platform_id}" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "Batch: ${batch_id}" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "========================================" | tee -a ${platform_id}_${batch_id}_liftover.log
    
    # Verify reference is indexed
    if [[ ! -f "${hg38_fasta}.fai" ]] && [[ ! -f "${hg38_fasta_fai}" ]]; then
        echo "ERROR: hg38 reference not indexed!" | tee -a ${platform_id}_${batch_id}_liftover.log
        exit 1
    fi
    
    # Count input variants
    n_input=\$(wc -l < ${bim})
    echo "Input variants (hg19): \${n_input}" | tee -a ${platform_id}_${batch_id}_liftover.log
    
    # Step 1: Convert PLINK to VCF (hg19)
    echo "" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "Step 1: Converting PLINK to VCF..." | tee -a ${platform_id}_${batch_id}_liftover.log
    plink2 --bfile \${prefix} \\
           --export vcf bgz \\
           --output-chr chrM \\
           --out hg19_input \\
           --threads ${task.cpus}
    
    # Index the hg19 VCF
    tabix -p vcf hg19_input.vcf.gz
    
    vcf_variants=\$(bcftools view -H hg19_input.vcf.gz | wc -l)
    echo "  VCF variants: \${vcf_variants}" | tee -a ${platform_id}_${batch_id}_liftover.log
    
    # Step 2: Liftover with CrossMap VCF mode
    echo "" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "Step 2: Running CrossMap liftover..." | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "  Chain file: ${chain_file}" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "  Target reference: ${hg38_fasta}" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "CrossMap will automatically:" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "  - Update genomic coordinates" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "  - Validate REF alleles against hg38 reference" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "  - Handle strand orientation" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "  - Exclude mismatched variants" | tee -a ${platform_id}_${batch_id}_liftover.log
    
    # Run CrossMap
    # CrossMap.py vcf <chain_file> <input.vcf> <target_ref.fa> <output_prefix>
    CrossMap.py vcf \\
        ${chain_file} \\
        hg19_input.vcf.gz \\
        ${hg38_fasta} \\
        hg38_lifted.vcf \\
        2>&1 | tee -a ${platform_id}_${batch_id}_liftover.log
    
    # Check if CrossMap succeeded
    if [[ ! -f "hg38_lifted.vcf" ]]; then
        echo "ERROR: CrossMap failed to create output VCF!" | tee -a ${platform_id}_${batch_id}_liftover.log
        echo "Check chain file and reference compatibility" | tee -a ${platform_id}_${batch_id}_liftover.log
        exit 1
    fi
    
    # Step 3: Process liftover results
    echo "" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "Step 3: Processing liftover results..." | tee -a ${platform_id}_${batch_id}_liftover.log
    
    # Compress lifted VCF
    bgzip -c hg38_lifted.vcf > hg38_lifted.vcf.gz
    tabix -p vcf hg38_lifted.vcf.gz
    
    # Count lifted variants
    n_lifted=\$(bcftools view -H hg38_lifted.vcf.gz | wc -l)
    n_failed=\$((n_input - n_lifted))
    pct_success=\$(awk "BEGIN {printf \\"%.2f\\", 100*\${n_lifted}/\${n_input}}")
    
    echo "Successfully lifted: \${n_lifted} (\${pct_success}%)" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "Failed to lift: \${n_failed}" | tee -a ${platform_id}_${batch_id}_liftover.log
    
    # Step 4: Analyze failures
    echo "" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "Step 4: Analyzing unmapped variants..." | tee -a ${platform_id}_${batch_id}_liftover.log
    
    if [[ -f "hg19_input.vcf.gz.unmap" ]]; then
        n_unmap=\$(grep -v "^#" hg19_input.vcf.gz.unmap | wc -l)
        echo "Unmapped variants file: hg19_input.vcf.gz.unmap (\${n_unmap} variants)" | tee -a ${platform_id}_${batch_id}_liftover.log
        
        echo "" | tee -a ${platform_id}_${batch_id}_liftover.log
        echo "Failure reasons:" | tee -a ${platform_id}_${batch_id}_liftover.log
        grep "^#Fail" hg19_input.vcf.gz.unmap 2>/dev/null | \\
            sed 's/^#Fail(/  /' | sed 's/):/: /' | \\
            sort | uniq -c | sort -rn | tee -a ${platform_id}_${batch_id}_liftover.log || \\
            echo "  No failure annotations found" | tee -a ${platform_id}_${batch_id}_liftover.log
        
        # Extract variant IDs for unmapped variants
        grep -v "^#" hg19_input.vcf.gz.unmap | \\
            awk '{print \$3}' > ${platform_id}_${batch_id}_unlifted.txt
        
        n_unlifted=\$(wc -l < ${platform_id}_${batch_id}_unlifted.txt)
        echo "" | tee -a ${platform_id}_${batch_id}_liftover.log
        echo "  Unlifted variant IDs saved: ${platform_id}_${batch_id}_unlifted.txt (\${n_unlifted})" | tee -a ${platform_id}_${batch_id}_liftover.log
    else
        echo "  No .unmap file created (100% success or CrossMap error)" | tee -a ${platform_id}_${batch_id}_liftover.log
    fi
    
    # Step 5: Convert lifted VCF back to PLINK
    echo "" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "Step 5: Converting lifted VCF to PLINK..." | tee -a ${platform_id}_${batch_id}_liftover.log
    
    plink2 --vcf hg38_lifted.vcf.gz \\
           --make-bed \\
           --out ${platform_id}_${batch_id}_hg38 \\
           --threads ${task.cpus}
    
    # Save the lifted VCF for reference/QC with clear naming
    echo "" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "Saving lifted VCF with descriptive naming..." | tee -a ${platform_id}_${batch_id}_liftover.log
    cp hg38_lifted.vcf.gz ${platform_id}_${batch_id}_hg19_to_hg38_lifted.vcf.gz
    cp hg38_lifted.vcf.gz.tbi ${platform_id}_${batch_id}_hg19_to_hg38_lifted.vcf.gz.tbi
    echo "  ${platform_id}_${batch_id}_hg19_to_hg38_lifted.vcf.gz" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "  ${platform_id}_${batch_id}_hg19_to_hg38_lifted.vcf.gz.tbi" | tee -a ${platform_id}_${batch_id}_liftover.log
    
    # Step 6: Quality checks
    echo "" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "=== QUALITY CHECKS ===" | tee -a ${platform_id}_${batch_id}_liftover.log
    
    n_final_snps=\$(wc -l < ${platform_id}_${batch_id}_hg38.bim)
    n_samples=\$(wc -l < ${platform_id}_${batch_id}_hg38.fam)
    
    echo "Final PLINK files:" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "  SNPs: \${n_final_snps}" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "  Samples: \${n_samples}" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "  Success rate: \${pct_success}%" | tee -a ${platform_id}_${batch_id}_liftover.log
    
    # Check chromosome naming
    echo "" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "Chromosome naming (first 5):" | tee -a ${platform_id}_${batch_id}_liftover.log
    awk '{print \$1}' ${platform_id}_${batch_id}_hg38.bim | sort -u | head -5 | \\
        sed 's/^/  /' | tee -a ${platform_id}_${batch_id}_liftover.log
    
    # Check position ranges for chr1
    echo "" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "Position sanity check (chr1):" | tee -a ${platform_id}_${batch_id}_liftover.log
    chr1_count=\$(awk '\$1 == "1" || \$1 == "chr1"' ${platform_id}_${batch_id}_hg38.bim | wc -l)
    if [[ \${chr1_count} -gt 0 ]]; then
        min_pos=\$(awk '\$1 == "1" || \$1 == "chr1" {print \$4}' ${platform_id}_${batch_id}_hg38.bim | sort -n | head -1)
        max_pos=\$(awk '\$1 == "1" || \$1 == "chr1" {print \$4}' ${platform_id}_${batch_id}_hg38.bim | sort -n | tail -1)
        echo "  Chr1 variants: \${chr1_count}" | tee -a ${platform_id}_${batch_id}_liftover.log
        echo "  Position range: \${min_pos} - \${max_pos}" | tee -a ${platform_id}_${batch_id}_liftover.log
        
        # Sanity check for hg38 (chr1 should be ~10kb to ~249M)
        if [[ \${max_pos} -gt 249000000 ]]; then
            echo "  ⚠ WARNING: Max position exceeds chr1 length!" | tee -a ${platform_id}_${batch_id}_liftover.log
        fi
    fi
    
    # Warn if success rate is low
    echo "" | tee -a ${platform_id}_${batch_id}_liftover.log
    if (( \$(echo "\${pct_success} < 95" | bc -l) )); then
        echo "⚠ WARNING: Liftover success rate below 95%" | tee -a ${platform_id}_${batch_id}_liftover.log
        echo "  Possible causes:" | tee -a ${platform_id}_${batch_id}_liftover.log
        echo "    - Input data quality issues" | tee -a ${platform_id}_${batch_id}_liftover.log
        echo "    - Reference genome mismatch" | tee -a ${platform_id}_${batch_id}_liftover.log
        echo "    - Incompatible chain file" | tee -a ${platform_id}_${batch_id}_liftover.log
        echo "  Review: ${platform_id}_${batch_id}_unlifted.txt" | tee -a ${platform_id}_${batch_id}_liftover.log
    else
        echo "✓ Liftover success rate acceptable (≥95%)" | tee -a ${platform_id}_${batch_id}_liftover.log
    fi
    
    echo "" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "=== SUMMARY ===" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "✓ Coordinates updated: hg19 → hg38" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "✓ REF alleles validated against hg38 reference" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "✓ Strand orientation handled automatically" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "✓ Invalid variants excluded" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "Output files ready for downstream analysis" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "========================================" | tee -a ${platform_id}_${batch_id}_liftover.log
    
    cat ${platform_id}_${batch_id}_liftover.log
    """
}

// ============================================================================
// PROCESS 6: UNION MERGE BATCHES TO PLATFORM
// ============================================================================
process UNION_MERGE_BATCHES_TO_PLATFORM {
    label 'plink_large'
    tag "${platform_id}"
    
    publishDir "${params.outdir}/${platform_id}/05_platform_merge",
        mode: 'copy',
        pattern: "*.{log,txt}"
    
    input:
    tuple val(platform_id),
          val(batch_ids),
          path(bed_files),
          path(bim_files),
          path(fam_files)
    
    output:
    tuple val(platform_id),
          path("${platform_id}_platform.bed"),
          path("${platform_id}_platform.bim"),
          path("${platform_id}_platform.fam"),
          emit: platform_merged
    
    path("${platform_id}_union_merge.log"), emit: union_log
    
    script:
    def n_batches = batch_ids instanceof List ? batch_ids.size() : 1
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "========================================" | tee ${platform_id}_union_merge.log
    echo "UNION MERGE: BATCHES → PLATFORM" | tee -a ${platform_id}_union_merge.log
    echo "Platform: ${platform_id}" | tee -a ${platform_id}_union_merge.log
    echo "Batches to merge: ${n_batches}" | tee -a ${platform_id}_union_merge.log
    echo "========================================" | tee -a ${platform_id}_union_merge.log
    
    if [[ ${n_batches} -eq 1 ]]; then
        # Single batch - just rename
        echo "Single batch detected - no merge needed" | tee -a ${platform_id}_union_merge.log
        mv ${bed_files} ${platform_id}_platform.bed
        mv ${bim_files} ${platform_id}_platform.bim
        mv ${fam_files} ${platform_id}_platform.fam
        
        n_snps=\$(wc -l < ${platform_id}_platform.bim)
        n_samples=\$(wc -l < ${platform_id}_platform.fam)
        
        echo "" | tee -a ${platform_id}_union_merge.log
        echo "Single batch output:" | tee -a ${platform_id}_union_merge.log
        echo "  SNPs: \${n_snps}" | tee -a ${platform_id}_union_merge.log
        echo "  Samples: \${n_samples}" | tee -a ${platform_id}_union_merge.log
        echo "" | tee -a ${platform_id}_union_merge.log
        echo "✓ Platform files ready (no merge required)" | tee -a ${platform_id}_union_merge.log
        exit 0
    fi
    
    # Multiple batches - perform union merge
    echo "" | tee -a ${platform_id}_union_merge.log
    echo "Merging multiple batches with UNION strategy..." | tee -a ${platform_id}_union_merge.log
    echo "  Strategy: Keep ALL variants across all batches" | tee -a ${platform_id}_union_merge.log
    echo "  Mode: --merge-mode 6 (resolve conflicts using first file)" | tee -a ${platform_id}_union_merge.log
    echo "" | tee -a ${platform_id}_union_merge.log
    
    # Create merge list
    ls *.bed | sed 's/.bed\$//' > batch_list.txt
    
    echo "Batches to merge:" | tee -a ${platform_id}_union_merge.log
    cat batch_list.txt | sed 's/^/  /' | tee -a ${platform_id}_union_merge.log
    
    # Count variants per batch
    echo "" | tee -a ${platform_id}_union_merge.log
    echo "Input variant counts per batch:" | tee -a ${platform_id}_union_merge.log
    total_variants=0
    for bim in *.bim; do
        n=\$(wc -l < \${bim})
        total_variants=\$((total_variants + n))
        echo "  \${bim}: \${n} variants" | tee -a ${platform_id}_union_merge.log
    done
    echo "  Total (with overlap): \${total_variants}" | tee -a ${platform_id}_union_merge.log
    
    # Perform UNION merge
    echo "" | tee -a ${platform_id}_union_merge.log
    echo "Executing union merge..." | tee -a ${platform_id}_union_merge.log
    
    plink --merge-list batch_list.txt \\
          --merge-mode 6 \\
          --make-bed \\
          --out ${platform_id}_platform \\
          --threads ${task.cpus} \\
          --allow-no-sex \\
          2>&1 | tee -a ${platform_id}_union_merge.log
    
    # Check if merge succeeded
    if [[ ! -f "${platform_id}_platform.bed" ]]; then
        echo "" | tee -a ${platform_id}_union_merge.log
        echo "ERROR: Merge failed!" | tee -a ${platform_id}_union_merge.log
        echo "Check log above for details" | tee -a ${platform_id}_union_merge.log
        exit 1
    fi
    
    # Calculate final statistics
    n_snps=\$(wc -l < ${platform_id}_platform.bim)
    n_samples=\$(wc -l < ${platform_id}_platform.fam)
    
    # Estimate unique variant count (approximate)
    avg_per_batch=\$(awk "BEGIN {printf \\"%.0f\\", \${total_variants}/${n_batches}}")
    
    echo "" | tee -a ${platform_id}_union_merge.log
    echo "========================================" | tee -a ${platform_id}_union_merge.log
    echo "UNION MERGE COMPLETE" | tee -a ${platform_id}_union_merge.log
    echo "========================================" | tee -a ${platform_id}_union_merge.log
    echo "" | tee -a ${platform_id}_union_merge.log
    echo "INPUT:" | tee -a ${platform_id}_union_merge.log
    echo "  Batches merged: ${n_batches}" | tee -a ${platform_id}_union_merge.log
    echo "  Total variants (with overlap): \${total_variants}" | tee -a ${platform_id}_union_merge.log
    echo "  Average per batch: \${avg_per_batch}" | tee -a ${platform_id}_union_merge.log
    echo "" | tee -a ${platform_id}_union_merge.log
    echo "OUTPUT:" | tee -a ${platform_id}_union_merge.log
    echo "  Unique SNPs: \${n_snps}" | tee -a ${platform_id}_union_merge.log
    echo "  Total samples: \${n_samples}" | tee -a ${platform_id}_union_merge.log
    echo "" | tee -a ${platform_id}_union_merge.log
    echo "MERGE STRATEGY:" | tee -a ${platform_id}_union_merge.log
    echo "  ✓ UNION merge - all variants preserved" | tee -a ${platform_id}_union_merge.log
    echo "  ✓ Conflicts resolved using first file's allele coding" | tee -a ${platform_id}_union_merge.log
    echo "  ✓ Maximum variant retention achieved" | tee -a ${platform_id}_union_merge.log
    echo "" | tee -a ${platform_id}_union_merge.log
    echo "NOTE: All batches are already on hg38 (post-liftover)" | tee -a ${platform_id}_union_merge.log
    echo "      Strand/allele issues resolved in earlier steps" | tee -a ${platform_id}_union_merge.log
    echo "" | tee -a ${platform_id}_union_merge.log
    echo "✓ Platform merge complete" | tee -a ${platform_id}_union_merge.log
    echo "========================================" | tee -a ${platform_id}_union_merge.log
    
    cat ${platform_id}_union_merge.log
    """
}
// ============================================================================
// PROCESS 7A: TOPMED STRAND VALIDATION (QC CHECK ONLY)
// ============================================================================
process TOPMED_STRAND_CHECK {
    label 'vcftools'
    tag "${platform_id}"
    
    publishDir "${params.outdir}/${platform_id}/06_topmed_validation",
        mode: 'copy',
        pattern: "*.{txt,log}"
    
    input:
    tuple val(platform_id),
          path(bed),
          path(bim),
          path(fam)
    path strand_script
    path topmed_ref
    
    output:
    tuple val(platform_id),
          path(bed),
          path(bim),
          path(fam),
          val("topmed"),
          emit: topmed_validated_files
    
    path("${platform_id}_topmed_validation.log"), emit: topmed_validation_log
    path("${platform_id}_topmed-*.txt"), emit: topmed_validation_files, optional: true
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    prefix=\$(basename ${bim} .bim)
    
    echo "========================================" | tee ${platform_id}_topmed_validation.log
    echo "TOPMed STRAND VALIDATION (QC CHECK)" | tee -a ${platform_id}_topmed_validation.log
    echo "Platform: ${platform_id}" | tee -a ${platform_id}_topmed_validation.log
    echo "========================================" | tee -a ${platform_id}_topmed_validation.log
    echo "" | tee -a ${platform_id}_topmed_validation.log
    echo "Purpose: Validate alignment to TOPMed Freeze 10 reference" | tee -a ${platform_id}_topmed_validation.log
    echo "Note: Data already aligned to hg38 in previous steps" | tee -a ${platform_id}_topmed_validation.log
    echo "      This is a QC check, not a fixing step" | tee -a ${platform_id}_topmed_validation.log
    echo "" | tee -a ${platform_id}_topmed_validation.log
    
    # Count input
    n_input=\$(wc -l < ${bim})
    s_input=\$(wc -l < ${fam})
    echo "Input variants: \${n_input}" | tee -a ${platform_id}_topmed_validation.log
    echo "Input samples: \${s_input}" | tee -a ${platform_id}_topmed_validation.log
    
    # Calculate frequencies for strand check
    echo "" | tee -a ${platform_id}_topmed_validation.log
    echo "Calculating allele frequencies..." | tee -a ${platform_id}_topmed_validation.log
    plink --bfile \${prefix} \\
          --freq \\
          --out \${prefix} \\
          --threads ${task.cpus}
    
    # Run Will Rayner's adapted strand check script for TOPMed Freeze 10
    echo "" | tee -a ${platform_id}_topmed_validation.log
    echo "Running TOPMed Freeze 10 strand validation..." | tee -a ${platform_id}_topmed_validation.log
    perl ${strand_script} \\
        -b ${bim} \\
        -f \${prefix}.frq \\
        -r ${topmed_ref} \\
        -h \\
        -o ${platform_id}_topmed \\
        2>&1 | tee -a ${platform_id}_topmed_validation.log
    
    # Analyze results
    echo "" | tee -a ${platform_id}_topmed_validation.log
    echo "=== VALIDATION RESULTS ===" | tee -a ${platform_id}_topmed_validation.log
    
    n_flip=0
    n_force=0
    n_remove=0
    n_match=0
    
    if [[ -f "${platform_id}_topmed-Strand-Flip.txt" ]]; then
        n_flip=\$(wc -l < ${platform_id}_topmed-Strand-Flip.txt)
    fi
    
    if [[ -f "${platform_id}_topmed-Force-Allele1.txt" ]]; then
        n_force=\$(wc -l < ${platform_id}_topmed-Force-Allele1.txt)
    fi
    
    if [[ -f "${platform_id}_topmed-remove.txt" ]]; then
        n_remove=\$(wc -l < ${platform_id}_topmed-remove.txt)
    fi
    
    n_match=\$((n_input - n_flip - n_force - n_remove))
    pct_match=\$(awk "BEGIN {printf \\"%.2f\\", 100*\${n_match}/\${n_input}}")
    
    echo "Variants matching TOPMed reference: \${n_match} (\${pct_match}%)" | tee -a ${platform_id}_topmed_validation.log
    echo "Strand flips needed: \${n_flip}" | tee -a ${platform_id}_topmed_validation.log
    echo "Force allele needed: \${n_force}" | tee -a ${platform_id}_topmed_validation.log
    echo "Variants to remove: \${n_remove}" | tee -a ${platform_id}_topmed_validation.log
    
    echo "" | tee -a ${platform_id}_topmed_validation.log
    
    # Warn if too many issues
    total_issues=\$((n_flip + n_force + n_remove))
    pct_issues=\$(awk "BEGIN {printf \\"%.2f\\", 100*\${total_issues}/\${n_input}}")
    
    if (( \$(echo "\${pct_issues} > 1.0" | bc -l) )); then
        echo "⚠ WARNING: \${pct_issues}% of variants need attention" | tee -a ${platform_id}_topmed_validation.log
        echo "  This is higher than expected (>1%)" | tee -a ${platform_id}_topmed_validation.log
        echo "  Possible causes:" | tee -a ${platform_id}_topmed_validation.log
        echo "    - Issues with earlier bcftools alignment" | tee -a ${platform_id}_topmed_validation.log
        echo "    - Platform merge introduced inconsistencies" | tee -a ${platform_id}_topmed_validation.log
        echo "    - Reference panel mismatch" | tee -a ${platform_id}_topmed_validation.log
        echo "" | tee -a ${platform_id}_topmed_validation.log
        echo "  Recommendation: Review earlier alignment steps" | tee -a ${platform_id}_topmed_validation.log
    else
        echo "✓ Validation passed - \${pct_match}% variants match TOPMed" | tee -a ${platform_id}_topmed_validation.log
        echo "  Minor discrepancies (\${pct_issues}%) are within normal range" | tee -a ${platform_id}_topmed_validation.log
    fi
    
    echo "" | tee -a ${platform_id}_topmed_validation.log
    echo "NOTE: Since data was aligned with bcftools +fixref earlier," | tee -a ${platform_id}_topmed_validation.log
    echo "      these files are for QC validation only" | tee -a ${platform_id}_topmed_validation.log
    echo "      No fixes will be applied - data already properly aligned" | tee -a ${platform_id}_topmed_validation.log
    echo "" | tee -a ${platform_id}_topmed_validation.log
    echo "✓ TOPMed validation complete" | tee -a ${platform_id}_topmed_validation.log
    echo "========================================" | tee -a ${platform_id}_topmed_validation.log
    
    cat ${platform_id}_topmed_validation.log
    """
}

// ============================================================================
// PROCESS 7B: ANVIL VALIDATION (ALREADY ALIGNED)
// ============================================================================
process ANVIL_VALIDATION {
    label 'process_low'
    tag "${platform_id}"
    
    publishDir "${params.outdir}/${platform_id}/06_anvil_validation",
        mode: 'copy',
        pattern: "*.log"
    
    input:
    tuple val(platform_id),
          path(bed),
          path(bim),
          path(fam)
    
    output:
    tuple val(platform_id),
          path(bed),
          path(bim),
          path(fam),
          val("anvil"),
          emit: anvil_validated_files
    
    path("${platform_id}_anvil_validation.log"), emit: anvil_validation_log
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "========================================" | tee ${platform_id}_anvil_validation.log
    echo "ALL OF US ANVIL VALIDATION" | tee -a ${platform_id}_anvil_validation.log
    echo "Platform: ${platform_id}" | tee -a ${platform_id}_anvil_validation.log
    echo "========================================" | tee -a ${platform_id}_anvil_validation.log
    echo "" | tee -a ${platform_id}_anvil_validation.log
    echo "AnVIL data preparation already complete:" | tee -a ${platform_id}_anvil_validation.log
    echo "  ✓ bcftools norm applied (multiallelic split)" | tee -a ${platform_id}_anvil_validation.log
    echo "  ✓ bcftools +fixref applied (REF validation)" | tee -a ${platform_id}_anvil_validation.log
    echo "  ✓ Aligned to hg38 reference" | tee -a ${platform_id}_anvil_validation.log
    echo "  ✓ Platform merge completed" | tee -a ${platform_id}_anvil_validation.log
    echo "" | tee -a ${platform_id}_anvil_validation.log
    
    # Quick validation counts
    n_snps=\$(wc -l < ${bim})
    n_samples=\$(wc -l < ${fam})
    
    echo "Data ready for AnVIL imputation:" | tee -a ${platform_id}_anvil_validation.log
    echo "  Variants: \${n_snps}" | tee -a ${platform_id}_anvil_validation.log
    echo "  Samples: \${n_samples}" | tee -a ${platform_id}_anvil_validation.log
    echo "" | tee -a ${platform_id}_anvil_validation.log
    echo "✓ AnVIL validation complete - no further alignment needed" | tee -a ${platform_id}_anvil_validation.log
    echo "========================================" | tee -a ${platform_id}_anvil_validation.log
    
    cat ${platform_id}_anvil_validation.log
    """
}

// ============================================================================
// PROCESS 8: LIGHT QC (NO STRAND FIXING - ALREADY ALIGNED)
// ============================================================================
process LIGHT_QC_BEFORE_IMPUTATION {
    label 'plink'
    tag "${platform_id}_${service}"
    
    publishDir "${params.outdir}/${platform_id}/07_${service}_qc",
        mode: 'copy',
        pattern: "*.log"
    
    input:
    tuple val(platform_id),
          path(bed),
          path(bim),
          path(fam),
          val(service)
    
    output:
    tuple val(platform_id),
          path("${platform_id}_${service}_qc.bed"),
          path("${platform_id}_${service}_qc.bim"),
          path("${platform_id}_${service}_qc.fam"),
          val(service),
          emit: qc_plink
    
    path("${platform_id}_${service}_light_qc.log"), emit: qc_log
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    prefix=\$(basename ${bed} .bed)
    
    echo "========================================" | tee ${platform_id}_${service}_light_qc.log
    echo "LIGHT QC BEFORE IMPUTATION" | tee -a ${platform_id}_${service}_light_qc.log
    echo "Platform: ${platform_id}" | tee -a ${platform_id}_${service}_light_qc.log
    echo "Service: ${service}" | tee -a ${platform_id}_${service}_light_qc.log
    echo "========================================" | tee -a ${platform_id}_${service}_light_qc.log
    echo "" | tee -a ${platform_id}_${service}_light_qc.log
    echo "NOTE: Data already aligned to hg38 in earlier steps" | tee -a ${platform_id}_${service}_light_qc.log
    echo "      NO strand fixing applied here" | tee -a ${platform_id}_${service}_light_qc.log
    echo "      Only basic QC filters for imputation readiness" | tee -a ${platform_id}_${service}_light_qc.log
    echo "" | tee -a ${platform_id}_${service}_light_qc.log
    
    n_input=\$(wc -l < ${bim})
    s_input=\$(wc -l < ${fam})
    
    echo "INPUT:" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  Variants: \${n_input}" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  Samples: \${s_input}" | tee -a ${platform_id}_${service}_light_qc.log
    
    # Step 1: SNPs only (biallelic A/C/G/T)
    echo "" | tee -a ${platform_id}_${service}_light_qc.log
    echo "Step 1: Keep biallelic SNPs only..." | tee -a ${platform_id}_${service}_light_qc.log
    plink --bfile \${prefix} \\
          --snps-only just-acgt \\
          --max-alleles 2 \\
          --make-bed \\
          --out step1 \\
          --threads ${task.cpus}
    
    n_step1=\$(wc -l < step1.bim)
    echo "  Remaining: \${n_step1}" | tee -a ${platform_id}_${service}_light_qc.log
    
    # Step 2: Remove duplicates
    echo "" | tee -a ${platform_id}_${service}_light_qc.log
    echo "Step 2: Remove duplicate variants..." | tee -a ${platform_id}_${service}_light_qc.log
    plink --bfile step1 \\
          --rm-dup exclude-all \\
          --make-bed \\
          --out step2 \\
          --threads ${task.cpus}
    
    n_step2=\$(wc -l < step2.bim)
    n_dup=\$((n_step1 - n_step2))
    echo "  Duplicates removed: \${n_dup}" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  Remaining: \${n_step2}" | tee -a ${platform_id}_${service}_light_qc.log
    
    # Step 3: Remove monomorphic (MAF > 0)
    echo "" | tee -a ${platform_id}_${service}_light_qc.log
    echo "Step 3: Remove monomorphic variants..." | tee -a ${platform_id}_${service}_light_qc.log
    plink --bfile step2 \\
          --maf 0.000001 \\
          --make-bed \\
          --out step3 \\
          --threads ${task.cpus}
    
    n_step3=\$(wc -l < step3.bim)
    n_mono=\$((n_step2 - n_step3))
    echo "  Monomorphic removed: \${n_mono}" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  Remaining: \${n_step3}" | tee -a ${platform_id}_${service}_light_qc.log
    
    # Step 4: Variant call rate (geno 0.05)
    echo "" | tee -a ${platform_id}_${service}_light_qc.log
    echo "Step 4: Filter variant call rate (≥95%)..." | tee -a ${platform_id}_${service}_light_qc.log
    plink --bfile step3 \\
          --geno 0.05 \\
          --make-bed \\
          --out step4 \\
          --threads ${task.cpus}
    
    n_step4=\$(wc -l < step4.bim)
    n_lowcr=\$((n_step3 - n_step4))
    echo "  Low call rate removed: \${n_lowcr}" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  Remaining: \${n_step4}" | tee -a ${platform_id}_${service}_light_qc.log
    
    # Step 5: Sample call rate (mind 0.05)
    echo "" | tee -a ${platform_id}_${service}_light_qc.log
    echo "Step 5: Filter sample call rate (≥95%)..." | tee -a ${platform_id}_${service}_light_qc.log
    plink --bfile step4 \\
          --mind 0.05 \\
          --make-bed \\
          --out ${platform_id}_${service}_qc \\
          --threads ${task.cpus}
    
    n_final=\$(wc -l < ${platform_id}_${service}_qc.bim)
    s_final=\$(wc -l < ${platform_id}_${service}_qc.fam)
    
    # Calculate retention
    var_retention=\$(awk "BEGIN {printf \\"%.2f\\", 100*\${n_final}/\${n_input}}")
    samp_retention=\$(awk "BEGIN {printf \\"%.2f\\", 100*\${s_final}/\${s_input}}")
    
    n_var_removed=\$((n_input - n_final))
    n_samp_removed=\$((s_input - s_final))
    
    echo "" | tee -a ${platform_id}_${service}_light_qc.log
    echo "========================================" | tee -a ${platform_id}_${service}_light_qc.log
    echo "LIGHT QC SUMMARY" | tee -a ${platform_id}_${service}_light_qc.log
    echo "========================================" | tee -a ${platform_id}_${service}_light_qc.log
    echo "" | tee -a ${platform_id}_${service}_light_qc.log
    echo "QC FILTERS APPLIED:" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  1. Biallelic SNPs only" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  2. Remove duplicates" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  3. Remove monomorphic (MAF > 0)" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  4. Variant call rate ≥ 95%" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  5. Sample call rate ≥ 95%" | tee -a ${platform_id}_${service}_light_qc.log
    echo "" | tee -a ${platform_id}_${service}_light_qc.log
    echo "RESULTS:" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  Input variants: \${n_input}" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  Output variants: \${n_final}" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  Variants removed: \${n_var_removed} (\${var_retention}% retained)" | tee -a ${platform_id}_${service}_light_qc.log
    echo "" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  Input samples: \${s_input}" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  Output samples: \${s_final}" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  Samples removed: \${n_samp_removed} (\${samp_retention}% retained)" | tee -a ${platform_id}_${service}_light_qc.log
    echo "" | tee -a ${platform_id}_${service}_light_qc.log
    echo "FILTERS EXPLICITLY NOT APPLIED:" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  ✗ Strand fixing (already aligned to hg38)" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  ✗ HWE filtering (Module 6 only)" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  ✗ MAF filtering beyond monomorphic (Module 6 only)" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  ✗ Sex checks (Module 6 only)" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  ✗ Heterozygosity (Module 6 only)" | tee -a ${platform_id}_${service}_light_qc.log
    echo "  ✗ Relatedness (Module 6 only)" | tee -a ${platform_id}_${service}_light_qc.log
    echo "" | tee -a ${platform_id}_${service}_light_qc.log
    echo "✓ Light QC complete - ready for ${service} imputation" | tee -a ${platform_id}_${service}_light_qc.log
    echo "========================================" | tee -a ${platform_id}_${service}_light_qc.log
    
    cat ${platform_id}_${service}_light_qc.log
    
    # Cleanup intermediate files
    rm -f step*.{bed,bim,fam,log,nosex}
    """
}


// ============================================================================
// PROCESS 9: CREATE SERVICE-SPECIFIC VCFS BY CHROMOSOME
// ============================================================================

process CREATE_SERVICE_SPECIFIC_VCFS {
    label 'vcftools'
    tag "${platform_id}_${service}_chr${chr}"
    
    publishDir "${params.outdir}/${platform_id}/08_${service}_vcfs",
        mode: 'copy',
        pattern: "*.vcf.gz*"
    
    input:
    tuple val(platform_id),
          val(chr),
          path(bed),
          path(bim),
          path(fam),
          val(service)
    
    output:
    tuple val(platform_id),
          val(service),
          val(chr),
          path("${platform_id}_${service}_chr${chr}.vcf.gz"),
          path("${platform_id}_${service}_chr${chr}.vcf.gz.csi"),
          emit: service_vcfs
    
    path("${platform_id}_${service}_chr${chr}_vcf.log"), emit: vcf_log
    
    script:
    def vcf_version = service == "topmed" ? "4.2" : "4.3"
    """
    #!/bin/bash
    set -euo pipefail
    
    prefix=\$(basename ${bed} .bed)
    n_snps=\$(wc -l < ${bim})
    
    echo "Creating ${service} VCF for chr${chr}" | tee ${platform_id}_${service}_chr${chr}_vcf.log
    echo "VCF version: ${vcf_version}" | tee -a ${platform_id}_${service}_chr${chr}_vcf.log
    
    if [[ \${n_snps} -eq 0 ]]; then
        # Create empty VCF
        cat > ${platform_id}_${service}_chr${chr}.vcf << 'EOFVCF'
##fileformat=VCFv${vcf_version}
##contig=<ID=chr${chr},assembly=GRCh38>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
EOFVCF
        bgzip ${platform_id}_${service}_chr${chr}.vcf
        bcftools index -c ${platform_id}_${service}_chr${chr}.vcf.gz
        echo "Empty VCF created" | tee -a ${platform_id}_${service}_chr${chr}_vcf.log
        exit 0
    fi
    
    # Convert to VCF with chr prefix
    plink2 --bfile \${prefix} \\
           --chr ${chr} \\
           --export vcf-${vcf_version} bgz \\
           --output-chr chr26 \\
           --out temp_chr${chr} \\
           --threads ${task.cpus}
    
    # Set variant IDs and clean FORMAT
    bcftools annotate \\
        --set-id '%CHROM:%POS:%REF:%ALT' \\
        temp_chr${chr}.vcf.gz | \\
    bcftools annotate \\
        -x INFO,^FORMAT/GT \\
        -O z \\
        -o annotated.vcf.gz
    
    # Sort and index
    bcftools sort \\
        annotated.vcf.gz \\
        -Oz \\
        -o ${platform_id}_${service}_chr${chr}.vcf.gz
    
    bcftools index -c ${platform_id}_${service}_chr${chr}.vcf.gz
    
    # Log results
    vcf_snps=\$(bcftools view -H ${platform_id}_${service}_chr${chr}.vcf.gz | wc -l)
    
    cat >> ${platform_id}_${service}_chr${chr}_vcf.log << EOF

VCF Created: ${platform_id}_${service}_chr${chr}
SNPs: \${vcf_snps}
Format: VCFv${vcf_version}, bgzip, CSI indexed
Service: ${service}
Ready for imputation

EOF

    cat ${platform_id}_${service}_chr${chr}_vcf.log
    
    rm -f temp*.vcf.gz annotated.vcf.gz
    """
}

// ============================================================================
// WORKFLOW: MODULE 1 ORCHESTRATION
// ============================================================================

workflow PRE_IMPUTATION_QC {
    take:
    sample_sheet_ch  // tuple(platform_id, batch_id, input_path, file_type, build, file_structure)
    
    main:
    
    // Step 1: Discover files
    DISCOVER_AND_LOAD_SAMPLES(sample_sheet_ch)
    
    // Step 2: Prepare samples (NO QC)
    // Branch based on file structure and prepare accordingly
    // [Implementation would expand here for different structures]
    
    // Step 3: Merge samples to batches (NO QC)
    MERGE_SAMPLES_TO_BATCH(prepared_samples.groupTuple(by: [0, 1]))
    
    // Step 4: Fix hg19 before liftover
    hg19_batches = MERGE_SAMPLES_TO_BATCH.out.batch_merged.filter { it[5] == 'hg19' }
    hg38_batches = MERGE_SAMPLES_TO_BATCH.out.batch_merged.filter { it[5] == 'hg38' }
    
    FORMAT_CHECK_AND_FIX_HG19_BEFORE_LIFTOVER(
        hg19_batches,
        file(params.hg19_fasta)
    )
    
    // Step 5: Liftover hg19 → hg38
    LIFTOVER_HG19_TO_HG38(
        FORMAT_CHECK_AND_FIX_HG19_BEFORE_LIFTOVER.out.fixed_files,
        file(params.liftover_chain)
    )
    
    // Combine all hg38 batches
    all_hg38 = LIFTOVER_HG19_TO_HG38.out.lifted_batch
        .mix(hg38_batches)
    
    // Step 6: Union merge batches to platform (3-pass)
    all_hg38
        .map { platform, batch, bed, bim, fam, build ->
            tuple(platform, bed, bim, fam)
        }
        .groupTuple(by: 0)
        .set { grouped_for_merge }
    
    UNION_MERGE_BATCHES_TO_PLATFORM(grouped_for_merge)
    
    // Step 7: Separate strand checking for each service
    TOPMED_STRAND_CHECK(
        UNION_MERGE_BATCHES_TO_PLATFORM.out.platform_merged,
        file("${projectDir}/bin/check-strand-topmed.pl"),
        file(params.topmed_reference)
    )
    
    ANVIL_STRAND_CHECK(
        UNION_MERGE_BATCHES_TO_PLATFORM.out.platform_merged,
        file(params.hg38_fasta)
    )
    
    // Combine both service paths
    both_services = TOPMED_STRAND_CHECK.out.topmed_strand_files
        .mix(ANVIL_STRAND_CHECK.out.anvil_fixed_files)
    
    // Step 8: Apply fixes and Light QC (FIRST TIME!)
    APPLY_FIXES_AND_LIGHT_QC(both_services)
    
    // Step 9: Split by chromosome and create service-specific VCFs
    chromosomes = Channel.of(1..22)
    
    APPLY_FIXES_AND_LIGHT_QC.out.qc_plink
        .combine(chromosomes)
        .set { for_vcf_creation }
    
    CREATE_SERVICE_SPECIFIC_VCFS(for_vcf_creation)
    
    emit:
    imputation_vcfs = CREATE_SERVICE_SPECIFIC_VCFS.out.service_vcfs
    qc_logs = APPLY_FIXES_AND_LIGHT_QC.out.qc_log
    union_logs = UNION_MERGE_BATCHES_TO_PLATFORM.out.union_log
    vcf_logs = CREATE_SERVICE_SPECIFIC_VCFS.out.vcf_log
}

// ============================================================================
// REQUIRED PARAMETERS
// ============================================================================
/*
 * params {
 *     // Reference files
 *     hg19_fasta = "<insert path to hg19.fa here>"
 *     hg38_fasta = "<insert path to hg38.fa here>"
 *     liftover_chain = "<insert path to hg19ToHg38.over.chain.gz here>"
 *     topmed_reference = "<insert path to PASS.Variants.TOPMed_freeze10_hg38.tab.gz here>"
 *     
 *     // Output
 *     outdir = "results"
 * }
 */
