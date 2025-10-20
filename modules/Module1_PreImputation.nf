#!/usr/bin/env nextflow

/*
================================================================================
MODULE 1: PRE-IMPUTATION QC AND PREPARATION (UPDATED)
================================================================================
UPDATES:
- REMOVED liftover process (servers handle it automatically)
- Added build tracking throughout pipeline
- Dual Rayner references (hg19 and hg38)
- Pass build info to Module 2 for server submission
================================================================================
*/

nextflow.enable.dsl = 2

// ============================================================================
// PROCESS 1.1: CONVERT INPUT TO PLINK FORMAT
// ============================================================================
process CONVERT_TO_PLINK {
    label 'preqc'
    tag "${platform_id}"
    conda "${projectDir}/envs/plink.yml"
    
    publishDir "${params.outdir}/${platform_id}/01_converted",
        mode: 'copy',
        pattern: "*.{bed,bim,fam,log}"
    
    input:
    tuple val(platform_id),
          val(file_path),
          val(file_type),
          val(build),       // Track build from input
          val(batch)
    
    output:
    tuple val(platform_id),
          path("${platform_id}.bed"),
          path("${platform_id}.bim"),
          path("${platform_id}.fam"),
          val(build),       // Pass build forward
          val(batch), emit: plink_files
    
    path("${platform_id}_convert.log"), emit: convert_log
    
    script:
    if (file_type == 'plink') {
        """
        # Already PLINK format - just copy/symlink
        ln -s ${file_path}.bed ${platform_id}.bed
        ln -s ${file_path}.bim ${platform_id}.bim
        ln -s ${file_path}.fam ${platform_id}.fam
        
        echo "Platform: ${platform_id}" > ${platform_id}_convert.log
        echo "Input format: PLINK" >> ${platform_id}_convert.log
        echo "Build: ${build}" >> ${platform_id}_convert.log
        echo "No conversion needed" >> ${platform_id}_convert.log
        """
    } else if (file_type == 'vcf') {
        """
        # Convert VCF to PLINK
        plink2 \\
            --vcf ${file_path} \\
            --make-bed \\
            --out ${platform_id} \\
            --threads ${task.cpus} \\
            --memory ${task.memory.toMega()}
        
        echo "Platform: ${platform_id}" > ${platform_id}_convert.log
        echo "Input format: VCF" >> ${platform_id}_convert.log
        echo "Build: ${build}" >> ${platform_id}_convert.log
        echo "Converted to PLINK format" >> ${platform_id}_convert.log
        """
    }
}

// ============================================================================
// PROCESS 1.2: BASIC QC WITH PLINK
// ============================================================================
process BASIC_QC {
    label 'preqc'
    tag "${platform_id}"
    conda "${projectDir}/envs/plink.yml"
    
    publishDir "${params.outdir}/${platform_id}/02_qc",
        mode: 'copy',
        pattern: "*.{txt,log}"
    
    input:
    tuple val(platform_id),
          path(bed),
          path(bim),
          path(fam),
          val(build),
          val(batch)
    
    output:
    tuple val(platform_id),
          path("${platform_id}_qc.bed"),
          path("${platform_id}_qc.bim"),
          path("${platform_id}_qc.fam"),
          val(build),       // Keep passing build
          val(batch), emit: qc_plink
    
    path("${platform_id}_qc_report.txt"), emit: qc_report
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Get input prefix
    input_prefix=\$(basename ${bed} .bed)
    
    echo "Starting Basic QC for ${platform_id}"
    echo "Build: ${build}"
    
    # Count input
    snps_input=\$(wc -l < ${bim})
    samples_input=\$(wc -l < ${fam})
    
    # Step 1: Keep SNPs only, remove duplicates
    plink \\
        --bfile \${input_prefix} \\
        --snps-only just-acgt \\
        --rm-dup exclude-all \\
        --make-bed \\
        --out ${platform_id}_step1 \\
        --threads ${task.cpus}
    
    snps_step1=\$(wc -l < ${platform_id}_step1.bim)
    removed_step1=\$((snps_input - snps_step1))
    
    # Step 2: Remove monomorphic SNPs
    plink \\
        --bfile ${platform_id}_step1 \\
        --maf 0.0 \\
        --make-bed \\
        --out ${platform_id}_step2 \\
        --threads ${task.cpus}
    
    snps_step2=\$(wc -l < ${platform_id}_step2.bim)
    removed_step2=\$((snps_step1 - snps_step2))
    
    # Step 3: Remove multi-allelic variants
    awk '\$5 != "0" && \$6 != "0"' ${platform_id}_step2.bim > snps_biallelic.txt
    
    plink \\
        --bfile ${platform_id}_step2 \\
        --extract snps_biallelic.txt \\
        --make-bed \\
        --out ${platform_id}_step3 \\
        --threads ${task.cpus}
    
    snps_step3=\$(wc -l < ${platform_id}_step3.bim)
    removed_step3=\$((snps_step2 - snps_step3))
    
    # Step 4: Missingness filtering
    plink \\
        --bfile ${platform_id}_step3 \\
        --geno ${params.preqc_geno} \\
        --mind ${params.preqc_mind} \\
        --make-bed \\
        --out ${platform_id}_qc \\
        --threads ${task.cpus}
    
    snps_final=\$(wc -l < ${platform_id}_qc.bim)
    samples_final=\$(wc -l < ${platform_id}_qc.fam)
    removed_snps_miss=\$((snps_step3 - snps_final))
    removed_samples=\$((samples_input - samples_final))
    
    # Generate QC report
    cat > ${platform_id}_qc_report.txt << EOF
===============================================
BASIC QC REPORT: ${platform_id}
===============================================

INPUT:
  Samples: \${samples_input}
  SNPs: \${snps_input}
  Build: ${build}
  Batch: ${batch}

FILTERING STEPS:
  Step 1 (SNPs only, duplicates):
    Removed: \${removed_step1} variants
    Remaining: \${snps_step1} variants
  
  Step 2 (Monomorphic):
    Removed: \${removed_step2} variants
    Remaining: \${snps_step2} variants
  
  Step 3 (Multi-allelic):
    Removed: \${removed_step3} variants
    Remaining: \${snps_step3} variants
  
  Step 4 (Missingness):
    SNP threshold: ${params.preqc_geno}
    Sample threshold: ${params.preqc_mind}
    Removed SNPs: \${removed_snps_miss}
    Removed samples: \${removed_samples}

FINAL RESULTS:
  Samples: \${samples_final}
  SNPs: \${snps_final}
  Retention rate (SNPs): \$(echo "scale=2; 100*\${snps_final}/\${snps_input}" | bc)%
  Retention rate (samples): \$(echo "scale=2; 100*\${samples_final}/\${samples_input}" | bc)%

BUILD INFO:
  Current build: ${build}
  Note: Imputation servers will convert to hg38 if needed

===============================================
EOF

    cat ${platform_id}_qc_report.txt
    
    # Cleanup
    rm -f ${platform_id}_step*.{bed,bim,fam,log,nosex}
    
    echo "Basic QC complete for ${platform_id}"
    """
}

// ============================================================================
// PROCESS 1.3: RAYNER STRAND CHECK (BUILD-AWARE)
// ============================================================================
process RAYNER_STRAND_CHECK {
    label 'preqc'
    tag "${platform_id}"
    conda "${projectDir}/envs/plink.yml"
    
    publishDir "${params.outdir}/${platform_id}/03_strand_check",
        mode: 'copy',
        pattern: "*.{txt,log,sh}"
    
    input:
    tuple val(platform_id),
          path(bed),
          path(bim),
          path(fam),
          val(build),
          val(batch)
    path(rayner_script)
    path(rayner_ref_hg19)
    path(rayner_ref_hg38)
    
    output:
    tuple val(platform_id),
          path("${platform_id}-updated*"),
          path("${platform_id}_qc*"),
          val(build),
          val(batch), emit: strand_files
    
    path("${platform_id}_rayner.log"), emit: rayner_log
    
    script:
    // Select appropriate reference based on build
    def ref_file = build == 'hg19' ? rayner_ref_hg19 : rayner_ref_hg38
    
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "Running Rayner strand check for ${platform_id}"
    echo "Build: ${build}"
    echo "Reference: ${ref_file}"
    
    # Copy input files
    cp ${bed} ${platform_id}_qc.bed
    cp ${bim} ${platform_id}_qc.bim
    cp ${fam} ${platform_id}_qc.fam
    
    # Run Rayner check with appropriate reference
    perl ${rayner_script} \\
        -b ${platform_id}_qc.bim \\
        -f ${ref_file} \\
        -h \\
        -o ${platform_id}
    
    # Create log
    cat > ${platform_id}_rayner.log << EOF
===============================================
RAYNER STRAND CHECK: ${platform_id}
===============================================

Build: ${build}
Reference: ${ref_file}

Generated files:
\$(ls -1 ${platform_id}-* 2>/dev/null || echo "No update files generated")

Strand flip script:
\$(ls -1 Run-plink.sh 2>/dev/null || echo "Run-plink.sh not found")

Note: Next step will apply these corrections

===============================================
EOF

    cat ${platform_id}_rayner.log
    
    echo "Rayner strand check complete"
    """
}

// ============================================================================
// PROCESS 1.4: APPLY STRAND FIXES
// ============================================================================
process APPLY_STRAND_FIXES {
    label 'preqc'
    tag "${platform_id}"
    conda "${projectDir}/envs/plink.yml"
    
    publishDir "${params.outdir}/${platform_id}/04_strand_fixed",
        mode: 'copy',
        pattern: "*.{bed,bim,fam,log}"
    
    input:
    tuple val(platform_id),
          path(update_files),
          path(input_files),
          val(build),
          val(batch)
    
    output:
    tuple val(platform_id),
          path("${platform_id}_fixed.bed"),
          path("${platform_id}_fixed.bim"),
          path("${platform_id}_fixed.fam"),
          val(build),
          val(batch), emit: fixed_plink
    
    path("${platform_id}_fixes.log"), emit: fixes_log
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "Applying strand corrections for ${platform_id}"
    
    # Check for Run-plink.sh script
    if [ -f "Run-plink.sh" ]; then
        echo "Found Run-plink.sh, executing corrections..."
        bash Run-plink.sh
        
        # Find the final output
        FINAL=\$(ls -t ${platform_id}-updated*.bed 2>/dev/null | head -1 | sed 's/.bed//')
        
        if [ -n "\$FINAL" ]; then
            mv \${FINAL}.bed ${platform_id}_fixed.bed
            mv \${FINAL}.bim ${platform_id}_fixed.bim
            mv \${FINAL}.fam ${platform_id}_fixed.fam
            echo "Applied corrections from Run-plink.sh" > ${platform_id}_fixes.log
        else
            echo "ERROR: No output from Run-plink.sh" > ${platform_id}_fixes.log
            exit 1
        fi
    else
        echo "No Run-plink.sh found, copying input files..."
        cp ${platform_id}_qc.bed ${platform_id}_fixed.bed
        cp ${platform_id}_qc.bim ${platform_id}_fixed.bim
        cp ${platform_id}_qc.fam ${platform_id}_fixed.fam
        echo "No strand corrections needed" > ${platform_id}_fixes.log
    fi
    
    # Count variants
    snps_fixed=\$(wc -l < ${platform_id}_fixed.bim)
    
    cat >> ${platform_id}_fixes.log << EOF

Platform: ${platform_id}
Build: ${build}
Final SNPs: \${snps_fixed}

EOF

    cat ${platform_id}_fixes.log
    """
}

// ============================================================================
// PROCESS 1.5: SPLIT BY CHROMOSOME
// ============================================================================
process SPLIT_BY_CHROMOSOME {
    label 'preqc'
    tag "${platform_id}_chr${chr}"
    conda "${projectDir}/envs/plink.yml"
    
    publishDir "${params.outdir}/${platform_id}/05_split_chr",
        mode: 'copy',
        pattern: "*.{bed,bim,fam}"
    
    input:
    tuple val(platform_id),
          path(bed),
          path(bim),
          path(fam),
          val(build),
          val(batch),
          val(chr)
    
    output:
    tuple val(platform_id),
          val(chr),
          path("${platform_id}_chr${chr}.bed"),
          path("${platform_id}_chr${chr}.bim"),
          path("${platform_id}_chr${chr}.fam"),
          val(build),
          val(batch), emit: chr_plink
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    input_prefix=\$(basename ${bed} .bed)
    
    plink \\
        --bfile \${input_prefix} \\
        --chr ${chr} \\
        --make-bed \\
        --out ${platform_id}_chr${chr} \\
        --threads ${task.cpus}
    
    echo "Split chromosome ${chr} for ${platform_id} (${build})"
    """
}

// ============================================================================
// PROCESS 1.6: CREATE VCFs FOR IMPUTATION (BUILD-AWARE)
// ============================================================================
process CREATE_IMPUTATION_VCFS {
    label 'bcftools'
    tag "${platform_id}_chr${chr}"
    conda "${projectDir}/envs/plink.yml"
    
    publishDir "${params.outdir}/${platform_id}/06_vcf_ready",
        mode: 'copy',
        pattern: "*.{vcf.gz,tbi,txt}"
    
    input:
    tuple val(platform_id),
          val(chr),
          path(bed),
          path(bim),
          path(fam),
          val(build),
          val(batch)
    
    output:
    tuple val(platform_id),
          val(chr),
          path("${platform_id}_chr${chr}.vcf.gz"),
          path("${platform_id}_chr${chr}.vcf.gz.tbi"),
          val(build),       // IMPORTANT: Pass build to Module 2
          val(batch), emit: imputation_vcfs
    
    path("${platform_id}_chr${chr}_metadata.txt"), emit: metadata
    
    script:
    // Chromosome naming based on build
    // hg38/GRCh38 typically uses "chr" prefix, hg19/GRCh37 may not
    // TOPMed expects chr prefix for hg38
    def chr_name = build == 'hg38' ? "chr${chr}" : chr
    
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "Creating VCF for ${platform_id} chromosome ${chr} (${build})"
    
    input_prefix=\$(basename ${bed} .bed)
    
    # Convert PLINK to VCF with appropriate chromosome naming
    plink2 \\
        --bfile \${input_prefix} \\
        --recode vcf-iid bgz \\
        --output-chr ${chr_name} \\
        --out ${platform_id}_chr${chr}_temp \\
        --threads ${task.cpus}
    
    # Sort and normalize
    bcftools sort \\
        ${platform_id}_chr${chr}_temp.vcf.gz \\
        -O z \\
        -o ${platform_id}_chr${chr}.vcf.gz
    
    # Index
    bcftools index -t ${platform_id}_chr${chr}.vcf.gz
    
    # Count variants
    variants=\$(bcftools view -H ${platform_id}_chr${chr}.vcf.gz | wc -l)
    samples=\$(bcftools query -l ${platform_id}_chr${chr}.vcf.gz | wc -l)
    
    # Create metadata file
    cat > ${platform_id}_chr${chr}_metadata.txt << EOF
===============================================
VCF METADATA: ${platform_id} Chr${chr}
===============================================

Platform: ${platform_id}
Chromosome: ${chr}
Chromosome name in VCF: ${chr_name}
Build: ${build}
Batch: ${batch}

Variants: \${variants}
Samples: \${samples}

File: ${platform_id}_chr${chr}.vcf.gz
Index: ${platform_id}_chr${chr}.vcf.gz.tbi

IMPORTANT FOR IMPUTATION:
- This VCF is in ${build} coordinates
- Imputation servers will automatically convert to hg38 if needed
- Specify build="${build}" in server submission

===============================================
EOF

    cat ${platform_id}_chr${chr}_metadata.txt
    
    # Cleanup
    rm -f ${platform_id}_chr${chr}_temp.vcf.gz
    
    echo "VCF creation complete"
    """
}

// ============================================================================
// MODULE 1 WORKFLOW (UPDATED - NO LIFTOVER)
// ============================================================================
workflow MODULE1_PREIMPUTATION {
    take:
    sample_sheet_ch  // tuple(platform_id, file_path, file_type, build, batch)
    
    main:
    
    // Process 1.1: Convert to PLINK if needed
    CONVERT_TO_PLINK(sample_sheet_ch)
    
    // Process 1.2: Basic QC
    BASIC_QC(CONVERT_TO_PLINK.out.plink_files)
    
    // ========================================================================
    // LIFTOVER REMOVED - Servers handle it!
    // ========================================================================
    // We keep track of build and pass it through the pipeline
    // The imputation servers will automatically convert hg19 → hg38
    
    // Process 1.3: Rayner strand check (build-aware)
    RAYNER_STRAND_CHECK(
        BASIC_QC.out.qc_plink,
        file(params.rayner_script),
        file(params.rayner_ref_hg19),
        file(params.rayner_ref_hg38)
    )
    
    // Process 1.4: Apply strand fixes
    APPLY_STRAND_FIXES(RAYNER_STRAND_CHECK.out.strand_files)
    
    // Process 1.5: Split by chromosome
    chromosomes = Channel.of(1..22)
    
    // Combine each platform with all chromosomes
    chr_input = APPLY_STRAND_FIXES.out.fixed_plink
        .combine(chromosomes)
    
    SPLIT_BY_CHROMOSOME(chr_input)
    
    // Process 1.6: Create VCFs for imputation
    CREATE_IMPUTATION_VCFS(SPLIT_BY_CHROMOSOME.out.chr_plink)
    
    emit:
    // Output VCFs WITH BUILD INFO for Module 2
    imputation_vcfs = CREATE_IMPUTATION_VCFS.out.imputation_vcfs
    
    // QC reports
    qc_reports = BASIC_QC.out.qc_report
    strand_reports = RAYNER_STRAND_CHECK.out.rayner_log
    fixes_reports = APPLY_STRAND_FIXES.out.fixes_log
    vcf_metadata = CREATE_IMPUTATION_VCFS.out.metadata
}

/*
================================================================================
USAGE NOTES:
================================================================================

1. Input sample sheet format (TSV):
   platform_id  file_path  file_type  build  batch
   Platform1    /path/to   plink      hg19   batch1
   Platform2    /path/to   plink      hg38   batch1

2. Build tracking:
   - Build information flows through entire pipeline
   - No liftover performed in Module 1
   - Build passed to Module 2 for server submission

3. Rayner references:
   - Need both hg19 and hg38 TOPMed references
   - Script automatically selects correct reference based on build

4. Output:
   - VCFs in original build (hg19 or hg38)
   - Build info passed to Module 2
   - Servers handle liftover automatically

================================================================================
*/
