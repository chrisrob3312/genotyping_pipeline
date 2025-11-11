#!/usr/bin/env nextflow

/*
 * ============================================================================
 * MODULE 1: PRE-IMPUTATION QC AND PREPARATION - OPTIMIZED v5.0
 * ============================================================================
 * 
 * KEY OPTIMIZATIONS IN THIS VERSION:
 * 1. Reduced format conversions (stay in PLINK when possible)
 * 2. Proper bcftools workflow: norm -m-any → norm --check-ref → sort → index
 * 3. PLINK used for union merging (bcftools can't do true union merges)
 * 4. Verbose flags added to Will Rayner script (-v)
 * 5. Better code organization while maintaining flexibility
 * 6. All Apptainer containers properly assigned
 * 
 * MAINTAINED FEATURES:
 * - Flexible input structures (individual_samples, individual_chr_split, 
 *   merged_batch, merged_chr_split)
 * - Mixed PLINK/VCF formats
 * - Per-batch genome build detection
 * - Union merge strategy (preserves max variants)
 * - Dual-pathway outputs (TOPMed: chr1-22 separate, AnVIL: concatenated)
 * 
 * WORKFLOW:
 * 1. Discover files → 2. Prepare samples → 3. Merge to batches
 * 4. Align to reference (hg19 or hg38) → 5. Liftover (if hg19)
 * 6. Union merge to platform → 7. Service validation
 * 8. Light QC → 9. Create service-specific VCFs
 * 
 * ============================================================================
 */

nextflow.enable.dsl = 2

// ============================================================================
// PROCESS 1: DISCOVER AND LOAD SAMPLES
// ============================================================================

process DISCOVER_AND_LOAD_SAMPLES {
    label 'process_low'
    container 'docker://quay.io/biocontainers/python:3.11'
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
    
    echo "=== DISCOVERING INPUT FILES ==="
    echo "Platform: ${platform_id} | Batch: ${batch_id}"
    echo "Structure: ${file_structure} | Type: ${file_type}"
    echo "Path: ${input_path}"
    
    case "${file_structure}" in
        individual_samples)
            if [[ "${file_type}" == "plink" ]]; then
                find ${input_path} -name "*.bed" -type f | sed 's/.bed\$//' > file_list.txt
            else
                find ${input_path} -name "*.vcf.gz" -type f > file_list.txt
            fi
            ;;
        
        individual_chr_split)
            if [[ "${file_type}" == "plink" ]]; then
                find ${input_path} -name "*_chr*.bed" -type f | \
                    sed 's/_chr[0-9]*.bed\$//' | sort -u > file_list.txt
            else
                find ${input_path} -name "*_chr*.vcf.gz" -type f | \
                    sed 's/_chr[0-9]*.vcf.gz\$//' | sort -u > file_list.txt
            fi
            ;;
        
        merged_batch|merged_chr_split)
            echo "${input_path}" > file_list.txt
            ;;
    esac
    
    n_files=\$(wc -l < file_list.txt)
    echo "✓ Found \${n_files} files/samples"
    
    [[ \${n_files} -eq 0 ]] && { echo "ERROR: No files found!"; exit 1; }
    
    head -5 file_list.txt
    """
}

// ============================================================================
// PROCESS 2: PREPARE SAMPLES (CONCATENATE CHROMOSOMES IF SPLIT)
// ============================================================================

process PREPARE_SAMPLES {
    label 'process_medium'
    container 'docker://quay.io/biocontainers/mulled-v2-b0f77c6e0af9c4ab47e2ad0f3e0bfa8042e2d2cb:a1bdb1c3dc6288db80ae3652db3c42ac6166b90b-0'
    // Contains: plink2, bcftools, bgzip, tabix
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
          val(file_structure),
          emit: prepared_samples
    
    path("${sample_id}_prep.log"), emit: prep_log
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "=== SAMPLE PREPARATION (NO QC) ===" | tee ${sample_id}_prep.log
    echo "Sample: ${sample_id}" | tee -a ${sample_id}_prep.log
    echo "Structure: ${file_structure}" | tee -a ${sample_id}_prep.log
    
    # Handle chromosome-split files
    if [[ "${file_structure}" =~ "chr_split" ]]; then
        echo "Concatenating chromosomes..." | tee -a ${sample_id}_prep.log
        
        if [[ "${file_type}" == "vcf" ]]; then
            # VCF: concat chr1-22
            for chr in {1..22}; do
                [[ -f "${sample_id}_chr\${chr}.vcf.gz" ]] && \
                    echo "${sample_id}_chr\${chr}.vcf.gz" >> vcf_list.txt
            done
            
            bcftools concat \\
                --file-list vcf_list.txt \\
                --output temp.vcf.gz \\
                --output-type z \\
                --threads ${task.cpus}
            
            plink2 --vcf temp.vcf.gz \\
                   --make-bed \\
                   --out ${sample_id} \\
                   --threads ${task.cpus}
        
        else
            # PLINK: merge chr1-22
            for chr in {1..22}; do
                [[ -f "${sample_id}_chr\${chr}.bed" ]] && \
                    echo "${sample_id}_chr\${chr}" >> merge_list.txt
            done
            
            plink --merge-list merge_list.txt \\
                  --make-bed \\
                  --out ${sample_id} \\
                  --threads ${task.cpus}
        fi
    
    # Handle single-file inputs
    elif [[ "${file_type}" == "vcf" ]]; then
        plink2 --vcf ${sample_id}.vcf.gz \\
               --make-bed \\
               --out ${sample_id} \\
               --threads ${task.cpus}
    
    else
        # Already PLINK format - files staged by Nextflow
        echo "Files already in PLINK format" | tee -a ${sample_id}_prep.log
    fi
    
    n_vars=\$(wc -l < ${sample_id}.bim)
    n_samp=\$(wc -l < ${sample_id}.fam)
    
    echo "Output: \${n_vars} variants, \${n_samp} samples" | tee -a ${sample_id}_prep.log
    echo "✓ Format conversion complete (NO QC)" | tee -a ${sample_id}_prep.log
    """
}

// ============================================================================
// PROCESS 3: MERGE SAMPLES TO BATCH (UNION MERGE)
// ============================================================================

process MERGE_SAMPLES_TO_BATCH {
    label 'plink_merge'
    container 'docker://quay.io/biocontainers/plink:1.90b6.21--h031d066_5'
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
          val(build),
          val(file_structure)
    
    output:
    tuple val(platform_id),
          val(batch_id),
          path("${platform_id}_${batch_id}.bed"),
          path("${platform_id}_${batch_id}.bim"),
          path("${platform_id}_${batch_id}.fam"),
          val(build),
          emit: batch_merged
    
    path("${platform_id}_${batch_id}_merge.log"), emit: merge_log
    
    when:
    file_structure != 'merged_batch' && file_structure != 'merged_chr_split'
    
    script:
    def n_files = (bed_files instanceof List) ? bed_files.size() : 1
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "=== MERGING SAMPLES TO BATCH (UNION) ===" | tee ${platform_id}_${batch_id}_merge.log
    echo "Platform: ${platform_id} | Batch: ${batch_id}" | tee -a ${platform_id}_${batch_id}_merge.log
    echo "Files: ${n_files}" | tee -a ${platform_id}_${batch_id}_merge.log
    
    if [[ ${n_files} -eq 1 ]]; then
        echo "Single file - no merge needed" | tee -a ${platform_id}_${batch_id}_merge.log
        mv ${bed_files} ${platform_id}_${batch_id}.bed
        mv ${bim_files} ${platform_id}_${batch_id}.bim
        mv ${fam_files} ${platform_id}_${batch_id}.fam
    else
        echo "Performing UNION merge (PLINK --merge-mode 6)..." | tee -a ${platform_id}_${batch_id}_merge.log
        
        ls *.bed | sed 's/.bed\$//' > merge_list.txt
        
        # Union merge keeps ALL variants from all files
        plink --merge-list merge_list.txt \\
              --merge-mode 6 \\
              --make-bed \\
              --out ${platform_id}_${batch_id} \\
              --threads ${task.cpus} \\
              --allow-no-sex
    fi
    
    n_vars=\$(wc -l < ${platform_id}_${batch_id}.bim)
    n_samp=\$(wc -l < ${platform_id}_${batch_id}.fam)
    
    echo "Output: \${n_vars} variants, \${n_samp} samples" | tee -a ${platform_id}_${batch_id}_merge.log
    echo "✓ Batch merge complete" | tee -a ${platform_id}_${batch_id}_merge.log
    """
}

// ============================================================================
// PROCESS 3B: PASSTHROUGH FOR PRE-MERGED BATCHES
// ============================================================================

process PASSTHROUGH_MERGED_BATCH {
    label 'process_low'
    container 'docker://quay.io/biocontainers/python:3.11'
    tag "${platform_id}_${batch_id}"
    
    publishDir "${params.outdir}/${platform_id}/${batch_id}/02_batch_merge",
        mode: 'copy',
        pattern: "*.log"
    
    input:
    tuple val(platform_id),
          val(batch_id),
          path(bed),
          path(bim),
          path(fam),
          val(build),
          val(file_structure)
    
    output:
    tuple val(platform_id),
          val(batch_id),
          path("${platform_id}_${batch_id}.bed"),
          path("${platform_id}_${batch_id}.bim"),
          path("${platform_id}_${batch_id}.fam"),
          val(build),
          emit: batch_merged
    
    path("${platform_id}_${batch_id}_passthrough.log"), emit: passthrough_log
    
    when:
    file_structure == 'merged_batch' || file_structure == 'merged_chr_split'
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "=== PRE-MERGED BATCH - PASSTHROUGH ===" | tee ${platform_id}_${batch_id}_passthrough.log
    echo "Platform: ${platform_id} | Batch: ${batch_id}" | tee -a ${platform_id}_${batch_id}_passthrough.log
    
    cp ${bed} ${platform_id}_${batch_id}.bed
    cp ${bim} ${platform_id}_${batch_id}.bim
    cp ${fam} ${platform_id}_${batch_id}.fam
    
    n_vars=\$(wc -l < ${platform_id}_${batch_id}.bim)
    n_samp=\$(wc -l < ${platform_id}_${batch_id}.fam)
    
    echo "Batch already merged: \${n_vars} variants, \${n_samp} samples" | tee -a ${platform_id}_${batch_id}_passthrough.log
    echo "✓ Passthrough complete" | tee -a ${platform_id}_${batch_id}_passthrough.log
    """
}

// ============================================================================
// PROCESS 4: ALIGN TO REFERENCE (hg19 or hg38)
// ============================================================================

process ALIGN_TO_REFERENCE {
    label 'vcftools'
    container 'docker://quay.io/biocontainers/bcftools:1.18--h8b25389_0'
    tag "${platform_id}_${batch_id}_${build}"
    
    publishDir "${params.outdir}/${platform_id}/${batch_id}/03_${build}_alignment",
        mode: 'copy',
        pattern: "*.{log,txt}"
    
    input:
    tuple val(platform_id),
          val(batch_id),
          path(bed),
          path(bim),
          path(fam),
          val(build)
    path fasta
    path fasta_fai
    
    output:
    tuple val(platform_id),
          val(batch_id),
          path("${platform_id}_${batch_id}_${build}_aligned.bed"),
          path("${platform_id}_${batch_id}_${build}_aligned.bim"),
          path("${platform_id}_${batch_id}_${build}_aligned.fam"),
          val(build),
          emit: aligned_batch
    
    path("${platform_id}_${batch_id}_${build}_alignment.log"), emit: alignment_log
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    prefix=\$(basename ${bed} .bed)
    
    echo "=== ALIGNING TO ${build} REFERENCE ===" | tee ${platform_id}_${batch_id}_${build}_alignment.log
    echo "Platform: ${platform_id} | Batch: ${batch_id}" | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    
    # Convert PLINK → VCF
    echo "Converting PLINK to VCF..." | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    plink2 --bfile \${prefix} \\
           --export vcf bgz \\
           --out temp \\
           --threads ${task.cpus}
    
    n_before=\$(bcftools view -H temp.vcf.gz | wc -l)
    echo "Variants before alignment: \${n_before}" | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    
    # Step 1: Split multiallelics
    echo "" | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    echo "Step 1: Splitting multiallelics..." | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    bcftools norm -m-any temp.vcf.gz \\
        -O z \\
        -o step1.vcf.gz \\
        --threads ${task.cpus} \\
        2>&1 | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    
    # Step 2: Check REF alleles and normalize
    echo "" | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    echo "Step 2: Checking REF alleles and normalizing..." | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    bcftools norm \\
        --check-ref w \\
        --fasta-ref ${fasta} \\
        step1.vcf.gz \\
        -O z \\
        -o step2.vcf.gz \\
        --threads ${task.cpus} \\
        2>&1 | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    
    # Step 3: Fix REF mismatches
    echo "" | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    echo "Step 3: Fixing REF mismatches..." | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    bcftools +fixref step2.vcf.gz \\
        --fasta-ref ${fasta} \\
        -O z \\
        -o step3.vcf.gz \\
        --threads ${task.cpus} \\
        2>&1 | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    
    # Step 4: Sort
    echo "" | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    echo "Step 4: Sorting..." | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    bcftools sort step3.vcf.gz \\
        -O z \\
        -o aligned.vcf.gz \\
        2>&1 | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    
    # Step 5: Index
    bcftools index -c aligned.vcf.gz
    
    n_after=\$(bcftools view -H aligned.vcf.gz | wc -l)
    n_removed=\$((n_before - n_after))
    retention=\$(awk "BEGIN {printf \\"%.2f\\", 100*\${n_after}/\${n_before}}")
    
    echo "" | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    echo "=== ALIGNMENT RESULTS ===" | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    echo "Before: \${n_before} | After: \${n_after} | Removed: \${n_removed} | Retention: \${retention}%" | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    
    # Convert back to PLINK
    echo "" | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    echo "Converting back to PLINK..." | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    plink2 --vcf aligned.vcf.gz \\
           --make-bed \\
           --out ${platform_id}_${batch_id}_${build}_aligned \\
           --threads ${task.cpus}
    
    echo "✓ Alignment to ${build} complete" | tee -a ${platform_id}_${batch_id}_${build}_alignment.log
    """
}

// ============================================================================
// PROCESS 5: LIFTOVER HG19 → HG38
// ============================================================================

process LIFTOVER_HG19_TO_HG38 {
    label 'liftover'
    container 'docker://quay.io/biocontainers/mulled-v2-27978155e7e54a01862fefce0fd465d66bb1a0bf:6c6f3dbbb52de08747e8c71fa27e11f9c23aa7e2-0'
    // Contains: CrossMap, plink2, bcftools, bgzip, tabix
    tag "${platform_id}_${batch_id}"
    
    publishDir "${params.outdir}/${platform_id}/${batch_id}/04_liftover",
        mode: 'copy',
        pattern: "*.{log,vcf.gz,vcf.gz.tbi,txt}"
    
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
    path("${platform_id}_${batch_id}_lifted.vcf.gz"), emit: lifted_vcf
    path("${platform_id}_${batch_id}_lifted.vcf.gz.tbi"), emit: lifted_vcf_index
    path("${platform_id}_${batch_id}_unlifted.txt"), emit: unlifted_variants, optional: true
    
    when:
    build == 'hg19'
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    prefix=\$(basename ${bed} .bed)
    
    echo "=== LIFTOVER: hg19 → hg38 ===" | tee ${platform_id}_${batch_id}_liftover.log
    echo "Platform: ${platform_id} | Batch: ${batch_id}" | tee -a ${platform_id}_${batch_id}_liftover.log
    
    n_input=\$(wc -l < ${bim})
    echo "Input variants (hg19): \${n_input}" | tee -a ${platform_id}_${batch_id}_liftover.log
    
    # Convert PLINK → VCF (hg19)
    echo "" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "Step 1: Converting to VCF..." | tee -a ${platform_id}_${batch_id}_liftover.log
    plink2 --bfile \${prefix} \\
           --export vcf bgz \\
           --output-chr chrM \\
           --out hg19_input \\
           --threads ${task.cpus}
    
    tabix -p vcf hg19_input.vcf.gz
    
    # CrossMap liftover (VCF mode with REF validation)
    echo "" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "Step 2: Running CrossMap liftover (VCF mode)..." | tee -a ${platform_id}_${batch_id}_liftover.log
    CrossMap.py vcf \\
        ${chain_file} \\
        hg19_input.vcf.gz \\
        ${hg38_fasta} \\
        hg38_lifted.vcf \\
        2>&1 | tee -a ${platform_id}_${batch_id}_liftover.log
    
    # Compress, sort, index
    bgzip -c hg38_lifted.vcf > temp.vcf.gz
    
    bcftools sort temp.vcf.gz \\
        -O z \\
        -o ${platform_id}_${batch_id}_lifted.vcf.gz
    
    tabix -p vcf ${platform_id}_${batch_id}_lifted.vcf.gz
    
    n_lifted=\$(bcftools view -H ${platform_id}_${batch_id}_lifted.vcf.gz | wc -l)
    n_failed=\$((n_input - n_lifted))
    success_pct=\$(awk "BEGIN {printf \\"%.2f\\", 100*\${n_lifted}/\${n_input}}")
    
    echo "" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "=== LIFTOVER RESULTS ===" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "Lifted: \${n_lifted} (\${success_pct}%) | Failed: \${n_failed}" | tee -a ${platform_id}_${batch_id}_liftover.log
    
    # Extract unlifted variants
    if [[ -f "hg19_input.vcf.gz.unmap" ]]; then
        grep -v "^#" hg19_input.vcf.gz.unmap | \\
            awk '{print \$3}' > ${platform_id}_${batch_id}_unlifted.txt || true
    fi
    
    # Convert lifted VCF → PLINK
    echo "" | tee -a ${platform_id}_${batch_id}_liftover.log
    echo "Step 3: Converting to PLINK..." | tee -a ${platform_id}_${batch_id}_liftover.log
    plink2 --vcf ${platform_id}_${batch_id}_lifted.vcf.gz \\
           --make-bed \\
           --out ${platform_id}_${batch_id}_hg38 \\
           --threads ${task.cpus}
    
    echo "✓ Liftover complete" | tee -a ${platform_id}_${batch_id}_liftover.log
    """
}

// ============================================================================
// PROCESS 6: UNION MERGE BATCHES → PLATFORM
// ============================================================================

process UNION_MERGE_TO_PLATFORM {
    label 'plink_large'
    container 'docker://quay.io/biocontainers/plink:1.90b6.21--h031d066_5'
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
    
    path("${platform_id}_platform_merge.log"), emit: merge_log
    
    script:
    def n_batches = (batch_ids instanceof List) ? batch_ids.size() : 1
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "=== UNION MERGE: BATCHES → PLATFORM ===" | tee ${platform_id}_platform_merge.log
    echo "Platform: ${platform_id} | Batches: ${n_batches}" | tee -a ${platform_id}_platform_merge.log
    
    if [[ ${n_batches} -eq 1 ]]; then
        echo "Single batch - no merge needed" | tee -a ${platform_id}_platform_merge.log
        mv ${bed_files} ${platform_id}_platform.bed
        mv ${bim_files} ${platform_id}_platform.bim
        mv ${fam_files} ${platform_id}_platform.fam
    else
        echo "Merging ${n_batches} batches with UNION strategy..." | tee -a ${platform_id}_platform_merge.log
        
        ls *.bed | sed 's/.bed\$//' > batch_list.txt
        
        # Log input variant counts
        echo "" | tee -a ${platform_id}_platform_merge.log
        for bim in *.bim; do
            n=\$(wc -l < \${bim})
            echo "  \${bim}: \${n} variants" | tee -a ${platform_id}_platform_merge.log
        done
        
        # Union merge
        plink --merge-list batch_list.txt \\
              --merge-mode 6 \\
              --make-bed \\
              --out ${platform_id}_platform \\
              --threads ${task.cpus} \\
              --allow-no-sex
    fi
    
    n_vars=\$(wc -l < ${platform_id}_platform.bim)
    n_samp=\$(wc -l < ${platform_id}_platform.fam)
    
    echo "" | tee -a ${platform_id}_platform_merge.log
    echo "Platform output: \${n_vars} variants, \${n_samp} samples" | tee -a ${platform_id}_platform_merge.log
    echo "✓ Platform merge complete" | tee -a ${platform_id}_platform_merge.log
    """
}

// ============================================================================
// PROCESS 7A: TOPMED VALIDATION WITH WILL RAYNER (VERBOSE!)
// ============================================================================

process TOPMED_VALIDATION {
    label 'vcftools'
    container 'docker://quay.io/biocontainers/mulled-v2-b0f77c6e0af9c4ab47e2ad0f3e0bfa8042e2d2cb:a1bdb1c3dc6288db80ae3652db3c42ac6166b90b-0'
    // Contains: plink, perl, bcftools
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
          emit: topmed_validated
    
    path("${platform_id}_topmed_validation.log"), emit: validation_log
    path("${platform_id}_topmed-*.txt"), emit: validation_files, optional: true
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    prefix=\$(basename ${bim} .bim)
    
    echo "=== TOPMed VALIDATION (QC CHECK) ===" | tee ${platform_id}_topmed_validation.log
    echo "Platform: ${platform_id}" | tee -a ${platform_id}_topmed_validation.log
    echo "Using: -n flag (keep ALL SNPs) + -v flag (verbose)" | tee -a ${platform_id}_topmed_validation.log
    
    n_input=\$(wc -l < ${bim})
    echo "Input: \${n_input} variants" | tee -a ${platform_id}_topmed_validation.log
    
    # Calculate frequencies
    echo "" | tee -a ${platform_id}_topmed_validation.log
    echo "Calculating allele frequencies..." | tee -a ${platform_id}_topmed_validation.log
    plink --bfile \${prefix} \\
          --freq \\
          --out \${prefix} \\
          --threads ${task.cpus}
    
    # Run Will Rayner with -n (no exclusion) and -v (verbose)
    echo "" | tee -a ${platform_id}_topmed_validation.log
    echo "Running Will Rayner script with verbose output..." | tee -a ${platform_id}_topmed_validation.log
    perl ${strand_script} \\
        -b ${bim} \\
        -f \${prefix}.frq \\
        -r ${topmed_ref} \\
        -h \\
        -n \\
        -v \\
        -o ${platform_id}_topmed \\
        2>&1 | tee -a ${platform_id}_topmed_validation.log
    
    # Summarize results
    echo "" | tee -a ${platform_id}_topmed_validation.log
    echo "=== VALIDATION SUMMARY ===" | tee -a ${platform_id}_topmed_validation.log
    
    [[ -f "${platform_id}_topmed-Strand-Flip.txt" ]] && \\
        echo "Strand flips noted: \$(wc -l < ${platform_id}_topmed-Strand-Flip.txt)" | tee -a ${platform_id}_topmed_validation.log
    
    [[ -f "${platform_id}_topmed-Force-Allele1.txt" ]] && \\
        echo "Force allele noted: \$(wc -l < ${platform_id}_topmed-Force-Allele1.txt)" | tee -a ${platform_id}_topmed_validation.log
    
    [[ -f "${platform_id}_topmed-remove.txt" ]] && \\
        echo "Flagged for removal: \$(wc -l < ${platform_id}_topmed-remove.txt)" | tee -a ${platform_id}_topmed_validation.log
    
    echo "" | tee -a ${platform_id}_topmed_validation.log
    echo "NOTE: With -n flag, these are QC checks only" | tee -a ${platform_id}_topmed_validation.log
    echo "      No variants excluded by AF differences" | tee -a ${platform_id}_topmed_validation.log
    echo "✓ TOPMed validation complete" | tee -a ${platform_id}_topmed_validation.log
    """
}

// ============================================================================
// PROCESS 7B: ANVIL VALIDATION (ALREADY ALIGNED)
// ============================================================================

process ANVIL_VALIDATION {
    label 'process_low'
    container 'docker://quay.io/biocontainers/python:3.11'
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
          emit: anvil_validated
    
    path("${platform_id}_anvil_validation.log"), emit: validation_log
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "=== ALL OF US ANVIL VALIDATION ===" | tee ${platform_id}_anvil_validation.log
    echo "Platform: ${platform_id}" | tee -a ${platform_id}_anvil_validation.log
    echo "Data already aligned to hg38 - no additional validation needed" | tee -a ${platform_id}_anvil_validation.log
    
    n_vars=\$(wc -l < ${bim})
    n_samp=\$(wc -l < ${fam})
    
    echo "Ready for AnVIL: \${n_vars} variants, \${n_samp} samples" | tee -a ${platform_id}_anvil_validation.log
    echo "✓ AnVIL validation complete" | tee -a ${platform_id}_anvil_validation.log
    """
}

// ============================================================================
// PROCESS 8: LIGHT QC (FIRST TIME QC IS APPLIED)
// ============================================================================

process LIGHT_QC {
    label 'plink'
    container 'docker://quay.io/biocontainers/plink:1.90b6.21--h031d066_5'
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
    
    path("${platform_id}_${service}_qc.log"), emit: qc_log
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    prefix=\$(basename ${bed} .bed)
    
    echo "=== LIGHT QC (FIRST QC APPLICATION) ===" | tee ${platform_id}_${service}_qc.log
    echo "Platform: ${platform_id} | Service: ${service}" | tee -a ${platform_id}_${service}_qc.log
    
    n_in_v=\$(wc -l < ${bim})
    n_in_s=\$(wc -l < ${fam})
    echo "Input: \${n_in_v} variants, \${n_in_s} samples" | tee -a ${platform_id}_${service}_qc.log
    
    # QC pipeline
    echo "" | tee -a ${platform_id}_${service}_qc.log
    
    # Step 1: Biallelic SNPs only
    echo "1. Keep biallelic SNPs (A/C/G/T)..." | tee -a ${platform_id}_${service}_qc.log
    plink --bfile \${prefix} \\
          --snps-only just-acgt \\
          --max-alleles 2 \\
          --make-bed --out step1 --threads ${task.cpus}
    
    # Step 2: Remove duplicates
    echo "2. Remove duplicates..." | tee -a ${platform_id}_${service}_qc.log
    plink --bfile step1 \\
          --rm-dup exclude-all \\
          --make-bed --out step2 --threads ${task.cpus}
    
    # Step 3: Remove monomorphic
    echo "3. Remove monomorphic (MAF > 0)..." | tee -a ${platform_id}_${service}_qc.log
    plink --bfile step2 \\
          --maf 0.000001 \\
          --make-bed --out step3 --threads ${task.cpus}
    
    # Step 4: Variant call rate ≥95%
    echo "4. Variant call rate ≥95%..." | tee -a ${platform_id}_${service}_qc.log
    plink --bfile step3 \\
          --geno 0.05 \\
          --make-bed --out step4 --threads ${task.cpus}
    
    # Step 5: Sample call rate ≥95%
    echo "5. Sample call rate ≥95%..." | tee -a ${platform_id}_${service}_qc.log
    plink --bfile step4 \\
          --mind 0.05 \\
          --make-bed --out ${platform_id}_${service}_qc --threads ${task.cpus}
    
    n_out_v=\$(wc -l < ${platform_id}_${service}_qc.bim)
    n_out_s=\$(wc -l < ${platform_id}_${service}_qc.fam)
    
    v_ret=\$(awk "BEGIN {printf \\"%.2f\\", 100*\${n_out_v}/\${n_in_v}}")
    s_ret=\$(awk "BEGIN {printf \\"%.2f\\", 100*\${n_out_s}/\${n_in_s}}")
    
    echo "" | tee -a ${platform_id}_${service}_qc.log
    echo "=== QC RESULTS ===" | tee -a ${platform_id}_${service}_qc.log
    echo "Variants: \${n_out_v} (\${v_ret}% retained)" | tee -a ${platform_id}_${service}_qc.log
    echo "Samples: \${n_out_s} (\${s_ret}% retained)" | tee -a ${platform_id}_${service}_qc.log
    echo "✓ Light QC complete - ready for ${service}" | tee -a ${platform_id}_${service}_qc.log
    
    rm -f step*.{bed,bim,fam,log,nosex}
    """
}

// ============================================================================
// PROCESS 9A: CREATE TOPMED VCFS (SEPARATE chr1-22)
// ============================================================================

process CREATE_TOPMED_VCFS {
    label 'vcftools'
    container 'docker://quay.io/biocontainers/bcftools:1.18--h8b25389_0'
    tag "${platform_id}_topmed_chr${chr}"
    
    publishDir "${params.outdir}/${platform_id}/08_topmed_vcfs",
        mode: 'copy',
        pattern: "*.vcf.gz*"
    
    input:
    tuple val(platform_id),
          val(chr),
          path(bed),
          path(bim),
          path(fam)
    
    output:
    tuple val(platform_id),
          val("topmed"),
          val(chr),
          path("${platform_id}_topmed_chr${chr}.vcf.gz"),
          path("${platform_id}_topmed_chr${chr}.vcf.gz.csi"),
          emit: topmed_vcfs
    
    path("${platform_id}_topmed_chr${chr}.log"), emit: vcf_log
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    prefix=\$(basename ${bed} .bed)
    
    echo "=== Creating TOPMed VCF: chr${chr} ===" | tee ${platform_id}_topmed_chr${chr}.log
    
    # Convert to VCF (VCF 4.2, chr prefix)
    plink2 --bfile \${prefix} \\
           --chr ${chr} \\
           --export vcf-4.2 bgz \\
           --output-chr chr26 \\
           --out temp_chr${chr} \\
           --threads ${task.cpus}
    
    # Set variant IDs and clean
    bcftools annotate \\
        --set-id '%CHROM:%POS:%REF:%ALT' \\
        temp_chr${chr}.vcf.gz | \\
    bcftools annotate \\
        -x INFO,^FORMAT/GT \\
        -O z \\
        -o annotated.vcf.gz
    
    # Sort and index
    bcftools sort annotated.vcf.gz \\
        -O z \\
        -o ${platform_id}_topmed_chr${chr}.vcf.gz
    
    bcftools index -c ${platform_id}_topmed_chr${chr}.vcf.gz
    
    n=\$(bcftools view -H ${platform_id}_topmed_chr${chr}.vcf.gz | wc -l)
    echo "VCF created: \${n} variants" | tee -a ${platform_id}_topmed_chr${chr}.log
    echo "Format: VCFv4.2, bgzip, CSI indexed" | tee -a ${platform_id}_topmed_chr${chr}.log
    
    rm -f temp*.vcf.gz annotated.vcf.gz
    """
}

// ============================================================================
// PROCESS 9B: CREATE ANVIL VCF (ALL chr1-22 CONCATENATED)
// ============================================================================

process CREATE_ANVIL_VCF {
    label 'vcftools'
    container 'docker://quay.io/biocontainers/bcftools:1.18--h8b25389_0'
    tag "${platform_id}_anvil"
    
    publishDir "${params.outdir}/${platform_id}/08_anvil_vcfs",
        mode: 'copy',
        pattern: "*.vcf.gz*"
    
    input:
    tuple val(platform_id),
          path(bed),
          path(bim),
          path(fam)
    
    output:
    tuple val(platform_id),
          val("anvil"),
          path("${platform_id}_anvil_all_autosomes.vcf.gz"),
          path("${platform_id}_anvil_all_autosomes.vcf.gz.csi"),
          emit: anvil_vcf
    
    path("${platform_id}_anvil.log"), emit: vcf_log
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    prefix=\$(basename ${bed} .bed)
    
    echo "=== Creating AnVIL VCF (ALL chr1-22) ===" | tee ${platform_id}_anvil.log
    echo "Platform: ${platform_id}" | tee -a ${platform_id}_anvil.log
    echo "CRITICAL: All autosomes in ONE file (saves quota!)" | tee -a ${platform_id}_anvil.log
    
    # Create per-chromosome VCFs
    echo "" | tee -a ${platform_id}_anvil.log
    echo "Creating per-chromosome VCFs..." | tee -a ${platform_id}_anvil.log
    
    for chr in {1..22}; do
        echo "  Processing chr\${chr}..." | tee -a ${platform_id}_anvil.log
        
        plink2 --bfile \${prefix} \\
               --chr \${chr} \\
               --export vcf-4.3 bgz \\
               --output-chr chr26 \\
               --out temp_chr\${chr} \\
               --threads ${task.cpus}
        
        bcftools annotate \\
            --set-id '%CHROM:%POS:%REF:%ALT' \\
            temp_chr\${chr}.vcf.gz | \\
        bcftools annotate \\
            -x INFO,^FORMAT/GT \\
            -O z \\
            -o chr\${chr}.vcf.gz
        
        bcftools index -c chr\${chr}.vcf.gz
        echo "chr\${chr}.vcf.gz" >> vcf_list.txt
    done
    
    # Concatenate all chromosomes
    echo "" | tee -a ${platform_id}_anvil.log
    echo "Concatenating chr1-22 into single VCF..." | tee -a ${platform_id}_anvil.log
    
    bcftools concat \\
        --file-list vcf_list.txt \\
        --output temp_concat.vcf.gz \\
        --output-type z \\
        --threads ${task.cpus}
    
    # Sort and index
    bcftools sort temp_concat.vcf.gz \\
        -O z \\
        -o ${platform_id}_anvil_all_autosomes.vcf.gz
    
    bcftools index -c ${platform_id}_anvil_all_autosomes.vcf.gz
    
    n=\$(bcftools view -H ${platform_id}_anvil_all_autosomes.vcf.gz | wc -l)
    
    echo "" | tee -a ${platform_id}_anvil.log
    echo "=== VCF CREATED ===" | tee -a ${platform_id}_anvil.log
    echo "File: ${platform_id}_anvil_all_autosomes.vcf.gz" | tee -a ${platform_id}_anvil.log
    echo "Total variants: \${n}" | tee -a ${platform_id}_anvil.log
    echo "Chromosomes: chr1-chr22 (all in one file)" | tee -a ${platform_id}_anvil.log
    echo "Format: VCFv4.3, bgzip, CSI indexed" | tee -a ${platform_id}_anvil.log
    echo "✓ Ready for AnVIL (single job = optimal quota)" | tee -a ${platform_id}_anvil.log
    
    rm -f temp*.vcf.gz chr[0-9]*.vcf.gz chr[0-9]*.vcf.gz.csi
    """
}

// ============================================================================
// WORKFLOW: MODULE 1 ORCHESTRATION
// ============================================================================

workflow PRE_IMPUTATION_QC {
    take:
    sample_sheet_ch
    
    main:
    
    // 1. Discover files
    DISCOVER_AND_LOAD_SAMPLES(sample_sheet_ch)
    
    // 2. Expand to individual samples and prepare
    DISCOVER_AND_LOAD_SAMPLES.out.discovered_files
        .flatMap { platform_id, batch_id, file_list, file_type, build, file_structure ->
            def samples = file(file_list).readLines()
            samples.collect { sample_path ->
                def sample_id = file(sample_path).getBaseName()
                tuple(platform_id, batch_id, sample_id, file_type, build, file_structure)
            }
        }
        .set { samples_to_prepare }
    
    PREPARE_SAMPLES(samples_to_prepare)
    
    // 3. Group samples by batch and route
    PREPARE_SAMPLES.out.prepared_samples
        .groupTuple(by: [0, 1])
        .branch {
            merge: (it[6] instanceof List ? it[6][0] : it[6]) !in ['merged_batch', 'merged_chr_split']
            passthrough: true
        }
        .set { batch_routing }
    
    MERGE_SAMPLES_TO_BATCH(batch_routing.merge)
    PASSTHROUGH_MERGED_BATCH(batch_routing.passthrough)
    
    // Combine all batches
    all_batches = MERGE_SAMPLES_TO_BATCH.out.batch_merged
        .mix(PASSTHROUGH_MERGED_BATCH.out.batch_merged)
    
    // 4. Split by genome build and align
    hg19_batches = all_batches.filter { it[5] == 'hg19' }
    hg38_batches = all_batches.filter { it[5] == 'hg38' }
    
    ALIGN_TO_REFERENCE(
        hg19_batches,
        file(params.hg19_fasta),
        file("${params.hg19_fasta}.fai")
    )
    
    ALIGN_TO_REFERENCE(
        hg38_batches,
        file(params.hg38_fasta),
        file("${params.hg38_fasta}.fai")
    )
    
    // 5. Liftover hg19 → hg38
    hg19_aligned = ALIGN_TO_REFERENCE.out.aligned_batch.filter { it[5] == 'hg19' }
    
    LIFTOVER_HG19_TO_HG38(
        hg19_aligned,
        file(params.liftover_chain),
        file(params.hg38_fasta),
        file("${params.hg38_fasta}.fai")
    )
    
    // 6. Combine all hg38 batches and merge to platform
    hg38_aligned = ALIGN_TO_REFERENCE.out.aligned_batch.filter { it[5] == 'hg38' }
    
    all_hg38 = LIFTOVER_HG19_TO_HG38.out.lifted_batch
        .mix(hg38_aligned)
        .map { platform_id, batch_id, bed, bim, fam, build ->
            tuple(platform_id, batch_id, bed, bim, fam)
        }
        .groupTuple(by: 0)
    
    UNION_MERGE_TO_PLATFORM(all_hg38)
    
    // 7. Service-specific validation
    TOPMED_VALIDATION(
        UNION_MERGE_TO_PLATFORM.out.platform_merged,
        file(params.strand_check_script),
        file(params.topmed_reference)
    )
    
    ANVIL_VALIDATION(UNION_MERGE_TO_PLATFORM.out.platform_merged)
    
    // 8. Light QC for both services
    both_services = TOPMED_VALIDATION.out.topmed_validated
        .mix(ANVIL_VALIDATION.out.anvil_validated)
    
    LIGHT_QC(both_services)
    
    // 9. Create service-specific VCFs
    
    // TOPMed: Separate chr1-22
    topmed_data = LIGHT_QC.out.qc_plink.filter { it[4] == "topmed" }
    chromosomes = Channel.of(1..22)
    
    topmed_data
        .combine(chromosomes)
        .map { platform_id, bed, bim, fam, service, chr ->
            tuple(platform_id, chr, bed, bim, fam)
        }
        .set { topmed_for_vcf }
    
    CREATE_TOPMED_VCFS(topmed_for_vcf)
    
    // AnVIL: Concatenated file
    anvil_data = LIGHT_QC.out.qc_plink
        .filter { it[4] == "anvil" }
        .map { platform_id, bed, bim, fam, service ->
            tuple(platform_id, bed, bim, fam)
        }
    
    CREATE_ANVIL_VCF(anvil_data)
    
    emit:
    topmed_vcfs = CREATE_TOPMED_VCFS.out.topmed_vcfs
    anvil_vcf = CREATE_ANVIL_VCF.out.anvil_vcf
    qc_logs = LIGHT_QC.out.qc_log
}
