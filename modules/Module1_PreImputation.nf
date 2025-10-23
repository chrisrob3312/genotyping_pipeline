/*
 * ===================================================================
 * MODULE 1: PRE-IMPUTATION QC AND PREPARATION (UPDATED v2.1)
 * ===================================================================
 * 
 * FEATURES:
 * - MERGE_PLATFORM_SAMPLES process for combining multiple files per platform
 * - Supports multiple PLINK or VCF files per platform (different batches/samples)
 * - Uses Apptainer/Singularity containers (no conda)
 * - CrossMap for genome liftover (hg19 → hg38)
 * - check-strand-topmed.pl from bin/ directory
 * - TOPMed Freeze 10 PASS variants for strand checking
 * - Strict VCF formatting for All of Us AnVIL + TOPMed
 * - Reference files from resources/references/ directory
 * 
 * PURPOSE:
 * Prepares genotyping data for imputation by:
 * 0. Merges multiple sample files per platform (if needed)
 * 1. Basic QC with PLINK
 * 2. Genome build conversion using CrossMap (hg19 → hg38)
 * 3. Strand checking with TOPMed Freeze 10 references
 * 4. Allele switch and strand flip corrections
 * 5. Strict VCF formatting for imputation servers
 * 
 * INPUT:
 * - Sample sheet (CSV) with platform, file paths (can be multiple), build, batch
 * - PLINK files (.bed/.bim/.fam) or VCF files
 * - Multiple files per platform will be automatically merged
 * 
 * OUTPUT:
 * - Per-chromosome VCF files (chr1-22)
 * - VCFv4.2, bgzip compressed, CSI indexed
 * - Chromosome naming: "chr1", "chr2", etc. (with "chr" prefix)
 * - Ready for TOPMed and All of Us AnVIL imputation
 * 
 * PARALLELIZATION:
 * - Each platform processed independently and in parallel
 * - Each chromosome split and processed in parallel
 */

// ===================================================================
// NEXTFLOW DSL2
// ===================================================================
nextflow.enable.dsl = 2

// ===================================================================
// PROCESS 1.0: MERGE PLATFORM SAMPLES (NEW!)
// ===================================================================
/*
 * WHAT IT DOES:
 * - Combines multiple PLINK or VCF files per platform
 * - Handles scenarios where samples are split across multiple files
 * - Merges by samples (not by SNPs) - keeps intersection of SNPs
 * 
 * WHEN IT RUNS:
 * - Only when platform has multiple input files
 * - For single files, passes through directly
 * 
 * EXAMPLE USE CASES:
 * - Platform1 has 3 batches: batch1.bed, batch2.bed, batch3.bed
 * - Platform2 has samples split by cohort
 * - Platform3 has incremental data additions
 * 
 * Container: PLINK2 for efficient merging
 */

process MERGE_PLATFORM_SAMPLES {
    container 'docker://quay.io/biocontainers/plink2:2.00a3.7--h4ac6f70_0'
    
    label 'merge_samples'
    tag "${platform_id}"
    
    publishDir "${params.outdir}/${platform_id}/00_sample_merge",
        mode: 'copy',
        pattern: "*.{log,txt}"
    
    input:
    tuple val(platform_id),
          val(file_paths),      // List of file paths (can be 1 or more)
          val(file_type),
          val(build),
          val(batch)
    
    output:
    tuple val(platform_id),
          path("${platform_id}_merged.bed"),
          path("${platform_id}_merged.bim"),
          path("${platform_id}_merged.fam"),
          val(build),
          val(batch), emit: merged_files
    
    path("${platform_id}_merge.log"), emit: merge_log
    
    script:
    def file_list = file_paths instanceof List ? file_paths : [file_paths]
    def n_files = file_list.size()
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "=========================================="
    echo "Merging samples for platform: ${platform_id}"
    echo "Number of input files: ${n_files}"
    echo "File type: ${file_type}"
    echo "Build: ${build}"
    echo "=========================================="
    
    # If only one file, convert to PLINK and pass through
    if [[ ${n_files} -eq 1 ]]; then
        echo "Single file detected - no merge needed"
        
        if [[ "${file_type}" == "plink" ]]; then
            # Link PLINK files
            file_path="${file_list[0]}"
            ln -s \${file_path}.bed ${platform_id}_merged.bed
            ln -s \${file_path}.bim ${platform_id}_merged.bim
            ln -s \${file_path}.fam ${platform_id}_merged.fam
            
        elif [[ "${file_type}" =~ ^vcf ]]; then
            # Convert VCF to PLINK
            file_path="${file_list[0]}"
            [[ ! -f "\${file_path}" ]] && file_path="${file_list[0]}.vcf.gz"
            [[ ! -f "\${file_path}" ]] && file_path="${file_list[0]}.vcf"
            
            plink2 \\
                --vcf \${file_path} \\
                --make-bed \\
                --out ${platform_id}_merged \\
                --threads ${task.cpus} \\
                --memory ${task.memory.toMega()}
        fi
        
        n_snps=\$(wc -l < ${platform_id}_merged.bim)
        n_samples=\$(wc -l < ${platform_id}_merged.fam)
        
        cat > ${platform_id}_merge.log << EOF
Platform: ${platform_id}
Input files: 1 (no merge needed)
Output SNPs: \${n_snps}
Output Samples: \${n_samples}
EOF
        
        echo "✓ Single file processed"
        exit 0
    fi
    
    echo "Multiple files detected - merging ${n_files} files"
    
    # CASE 1: Multiple PLINK files
    if [[ "${file_type}" == "plink" ]]; then
        echo "Converting all PLINK files and preparing for merge..."
        
        # Create list of files to merge
        first_file="${file_list[0]}"
        echo "\${first_file}" > merge_list.txt
        
        # Add remaining files to merge list
        for i in \$(seq 1 \$((${n_files} - 1))); do
            file_path="${file_list[\$i]}"
            echo "\${file_path}" >> merge_list.txt
        done
        
        # Check all files exist
        while read -r file_path; do
            if [[ ! -f "\${file_path}.bed" ]]; then
                echo "ERROR: File not found: \${file_path}.bed"
                exit 1
            fi
        done < merge_list.txt
        
        # Perform merge using PLINK2
        # --pmerge-list: merge by samples (not SNPs)
        # Keeps intersection of SNPs across all files
        plink2 \\
            --pfile-merge-list merge_list.txt \\
            --make-bed \\
            --out ${platform_id}_merged \\
            --threads ${task.cpus} \\
            --memory ${task.memory.toMega()}
    
    # CASE 2: Multiple VCF files
    elif [[ "${file_type}" =~ ^vcf ]]; then
        echo "Converting VCF files to PLINK and merging..."
        
        # Convert each VCF to PLINK first
        for i in \$(seq 0 \$((${n_files} - 1))); do
            file_path="${file_list[\$i]}"
            
            # Find VCF file (try different extensions)
            vcf_file="\${file_path}"
            [[ ! -f "\${vcf_file}" ]] && vcf_file="\${file_path}.vcf.gz"
            [[ ! -f "\${vcf_file}" ]] && vcf_file="\${file_path}.vcf"
            
            if [[ ! -f "\${vcf_file}" ]]; then
                echo "ERROR: VCF file not found: \${file_path}"
                exit 1
            fi
            
            echo "Converting VCF file \$i: \${vcf_file}"
            
            plink2 \\
                --vcf \${vcf_file} \\
                --make-bed \\
                --out temp_\${i} \\
                --threads ${task.cpus}
            
            echo "temp_\${i}" >> merge_list.txt
        done
        
        # Merge all converted PLINK files
        plink2 \\
            --pfile-merge-list merge_list.txt \\
            --make-bed \\
            --out ${platform_id}_merged \\
            --threads ${task.cpus} \\
            --memory ${task.memory.toMega()}
        
        # Cleanup temp files
        rm -f temp_*.{bed,bim,fam,log}
    
    else
        echo "ERROR: Unknown file type: ${file_type}"
        exit 1
    fi
    
    # Verify merge was successful
    if [[ ! -f "${platform_id}_merged.bed" ]]; then
        echo "ERROR: Merge failed - output files not created"
        exit 1
    fi
    
    # Calculate statistics
    n_snps=\$(wc -l < ${platform_id}_merged.bim)
    n_samples=\$(wc -l < ${platform_id}_merged.fam)
    
    # Check for duplicate samples
    n_dup_samples=\$(cut -f2 ${platform_id}_merged.fam | sort | uniq -d | wc -l)
    
    if [[ \${n_dup_samples} -gt 0 ]]; then
        echo "WARNING: Found \${n_dup_samples} duplicate sample IDs after merge"
        cut -f2 ${platform_id}_merged.fam | sort | uniq -d > duplicate_samples.txt
        cat duplicate_samples.txt
    fi
    
    # Generate merge report
    cat > ${platform_id}_merge.log << EOF
===============================================
SAMPLE MERGE REPORT: ${platform_id}
===============================================

INPUT:
  Number of files: ${n_files}
  File type: ${file_type}
  Build: ${build}
  
MERGE STRATEGY:
  Method: Sample-based merge
  SNP handling: Intersection (common SNPs only)
  
OUTPUT:
  Total samples: \${n_samples}
  Total SNPs: \${n_snps}
  Duplicate sample IDs: \${n_dup_samples}
  
FILES:
  ${platform_id}_merged.{bed,bim,fam}

MERGE STATUS: ✓ SUCCESS

===============================================
EOF

    cat ${platform_id}_merge.log
    
    echo "✓ Platform merge complete"
    echo "  Samples: \${n_samples}"
    echo "  SNPs: \${n_snps}"
    """
}



// ===================================================================
// PROCESS 1.1: BASIC QC
// ===================================================================

process BASIC_QC {
    container 'docker://quay.io/biocontainers/plink:1.90b6.21--h031d066_5'
    
    label 'preqc'
    tag "${platform_id}"
    
    publishDir "${params.outdir}/${platform_id}/01_basic_qc", 
        mode: 'copy',
        pattern: "*.{log,txt}"
    
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
          val(build),
          val(batch), emit: qc_plink
    
    path("${platform_id}_qc_report.txt"), emit: qc_report
    
    script:
    def geno = params.preqc_geno ?: '0.05'
    def mind = params.preqc_mind ?: '0.05'
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "=========================================="
    echo "Basic QC: ${platform_id}"
    echo "Build: ${build}"
    echo "=========================================="
    
    input_prefix=\$(basename ${bed} .bed)
    snps_input=\$(wc -l < ${bim})
    samples_input=\$(wc -l < ${fam})
    
    # Step 1: Biallelic SNPs only
    echo "Step 1: Filtering to biallelic SNPs..."
    plink \\
        --bfile \${input_prefix} \\
        --snps-only just-acgt \\
        --rm-dup exclude-all \\
        --make-bed \\
        --out temp1 \\
        --threads ${task.cpus}
    
    snps_step1=\$(wc -l < temp1.bim)
    
    # Step 2: Remove monomorphic
    echo "Step 2: Removing monomorphic SNPs..."
    plink \\
        --bfile temp1 \\
        --maf 0.000001 \\
        --make-bed \\
        --out temp2 \\
        --threads ${task.cpus}
    
    snps_step2=\$(wc -l < temp2.bim)
    
    # Step 3: Missingness filters
    echo "Step 3: Applying missingness filters..."
    plink \\
        --bfile temp2 \\
        --geno ${geno} \\
        --mind ${mind} \\
        --make-bed \\
        --out ${platform_id}_qc \\
        --threads ${task.cpus}
    
    snps_final=\$(wc -l < ${platform_id}_qc.bim)
    samples_final=\$(wc -l < ${platform_id}_qc.fam)
    
    snp_retention=\$(awk "BEGIN {printf \\"%.2f\\", 100*\${snps_final}/\${snps_input}}")
    sample_retention=\$(awk "BEGIN {printf \\"%.2f\\", 100*\${samples_final}/\${samples_input}}")
    
    # Report
    cat > ${platform_id}_qc_report.txt << EOF
===============================================
BASIC QC REPORT: ${platform_id}
===============================================

INPUT:
  Samples: \${samples_input}
  SNPs: \${snps_input}
  Build: ${build}

FILTERING:
  1. Biallelic SNPs: \${snps_step1}
  2. Remove monomorphic: \${snps_step2}
  3. Missingness (geno=${geno}, mind=${mind}): \${snps_final} SNPs, \${samples_final} samples

OUTPUT:
  Samples: \${samples_final} (\${sample_retention}% retained)
  SNPs: \${snps_final} (\${snp_retention}% retained)
  
QC PASS: ✓

===============================================
EOF

    cat ${platform_id}_qc_report.txt
    rm -f temp*.{bed,bim,fam,log,nosex}
    """
}

// ===================================================================
// PROCESS 1.2: LIFTOVER WITH CROSSMAP
// ===================================================================

process LIFTOVER_CROSSMAP {
    container 'docker://quay.io/biocontainers/crossmap:0.6.4--pyh7cba7a3_0'
    
    label 'liftover'
    tag "${platform_id}"
    
    publishDir "${params.outdir}/${platform_id}/02_liftover",
        mode: 'copy',
        pattern: "*.{log,txt}"
    
    when:
    build == 'hg19'
    
    input:
    tuple val(platform_id),
          path(bed),
          path(bim),
          path(fam),
          val(build),
          val(batch)
    path chain_file
    path hg38_fasta
    
    output:
    tuple val(platform_id),
          path("${platform_id}_hg38.bed"),
          path("${platform_id}_hg38.bim"),
          path("${platform_id}_hg38.fam"),
          val("hg38"),
          val(batch), emit: lifted_plink
    
    path("${platform_id}_liftover.log"), emit: liftover_log
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "=========================================="
    echo "CrossMap Liftover: ${platform_id}"
    echo "hg19 → hg38"
    echo "=========================================="
    
    input_prefix=\$(basename ${bed} .bed)
    snps_input=\$(wc -l < ${bim})
    
    # Convert .bim to BED format
    awk 'BEGIN {OFS="\\t"} {
        chr = (\$1 ~ /^chr/) ? \$1 : "chr"\$1
        print chr, \$4-1, \$4, \$2, \$5, \$6
    }' ${bim} > hg19_positions.bed
    
    # Run CrossMap
    CrossMap.py bed ${chain_file} hg19_positions.bed hg38_positions.bed || true
    
    # Process results
    awk '{print \$4}' hg38_positions.bed > mapped_snps.txt
    snps_mapped=\$(wc -l < mapped_snps.txt)
    pct_mapped=\$(awk "BEGIN {printf \\"%.2f\\", 100*\${snps_mapped}/\${snps_input}}")
    
    # Create mapping
    awk 'BEGIN {OFS="\\t"} {
        chr = \$1; gsub(/^chr/, "", chr)
        print \$4, chr, \$3
    }' hg38_positions.bed > snp_to_hg38_pos.txt
    
    # Update .bim
    awk 'BEGIN {OFS="\\t"} 
    NR==FNR {map[\$1] = \$2"\\t"\$3; next}
    \$2 in map {
        split(map[\$2], pos, "\\t")
        print pos[1], \$2, \$3, pos[2], \$5, \$6
    }' snp_to_hg38_pos.txt ${bim} > ${platform_id}_hg38.bim
    
    # Extract mapped SNPs
    plink \\
        --bfile \${input_prefix} \\
        --extract mapped_snps.txt \\
        --make-bed \\
        --out ${platform_id}_hg38_temp \\
        --threads ${task.cpus}
    
    mv ${platform_id}_hg38_temp.bed ${platform_id}_hg38.bed
    mv ${platform_id}_hg38_temp.fam ${platform_id}_hg38.fam
    
    # Report
    cat > ${platform_id}_liftover.log << EOF
===============================================
CROSSMAP LIFTOVER: ${platform_id}
hg19 → hg38
===============================================

INPUT (hg19): \${snps_input} SNPs
MAPPED (hg38): \${snps_mapped} SNPs (\${pct_mapped}%)
UNMAPPED: \$((snps_input - snps_mapped)) SNPs

LIFTOVER PASS: ✓
===============================================
EOF

    cat ${platform_id}_liftover.log
    """
}

// ===================================================================
// PROCESS 1.3: STRAND CHECKING WITH TOPMED FREEZE 10
// ===================================================================

process STRAND_CHECK_TOPMED {
    container 'docker://perl:5.34'
    
    label 'strand_check'
    tag "${platform_id}"
    
    publishDir "${params.outdir}/${platform_id}/03_strand_check",
        mode: 'copy',
        pattern: "*.{txt,log}"
    
    input:
    tuple val(platform_id),
          path(bed),
          path(bim),
          path(fam),
          val(build),
          val(batch)
    path strand_script
    path topmed_ref
    
    output:
    tuple val(platform_id),
          path(bed),
          path(bim),
          path(fam),
          val(build),
          val(batch),
          path("${platform_id}-*.txt"), emit: strand_files
    
    path("${platform_id}_strand_check.log"), emit: strand_log
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "=========================================="
    echo "Strand Check: ${platform_id}"
    echo "Reference: TOPMed Freeze 10"
    echo "=========================================="
    
    input_prefix=\$(basename ${bim} .bim)
    
    # Calculate frequencies
    plink \\
        --bfile \${input_prefix} \\
        --freq \\
        --out \${input_prefix} \\
        --threads ${task.cpus}
    
    # Run strand check
    perl ${strand_script} \\
        -b ${bim} \\
        -f \${input_prefix}.frq \\
        -r ${topmed_ref} \\
        -h \\
        -o ${platform_id}
    
    # Count results
    n_flip=0; [[ -f "${platform_id}-Strand-Flip.txt" ]] && n_flip=\$(wc -l < ${platform_id}-Strand-Flip.txt)
    n_force=0; [[ -f "${platform_id}-Force-Allele1.txt" ]] && n_force=\$(wc -l < ${platform_id}-Force-Allele1.txt)
    n_remove=0; [[ -f "${platform_id}-remove.txt" ]] && n_remove=\$(wc -l < ${platform_id}-remove.txt)
    
    cat > ${platform_id}_strand_check.log << EOF
===============================================
STRAND CHECK: ${platform_id}
Reference: TOPMed Freeze 10
===============================================

RESULTS:
  Flip: \${n_flip}
  Force: \${n_force}
  Remove: \${n_remove}
  
STRAND CHECK PASS: ✓
===============================================
EOF

    cat ${platform_id}_strand_check.log
    """
}

// ===================================================================
// PROCESS 1.4: APPLY STRAND FIXES
// ===================================================================

process APPLY_STRAND_FIXES {
    container 'docker://quay.io/biocontainers/plink:1.90b6.21--h031d066_5'
    
    label 'preqc'
    tag "${platform_id}"
    
    publishDir "${params.outdir}/${platform_id}/04_strand_fixed",
        mode: 'copy',
        pattern: "*.log"
    
    input:
    tuple val(platform_id),
          path(bed),
          path(bim),
          path(fam),
          val(build),
          val(batch),
          path(strand_files)
    
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
    
    input_prefix=\$(basename ${bed} .bed)
    
    # Step 1: Remove
    if [[ -f "${platform_id}-remove.txt" && -s "${platform_id}-remove.txt" ]]; then
        plink --bfile \${input_prefix} --exclude ${platform_id}-remove.txt --make-bed --out temp1 --threads ${task.cpus}
    else
        ln -s ${bed} temp1.bed; ln -s ${bim} temp1.bim; ln -s ${fam} temp1.fam
    fi
    
    # Step 2: Flip
    if [[ -f "${platform_id}-Strand-Flip.txt" && -s "${platform_id}-Strand-Flip.txt" ]]; then
        plink --bfile temp1 --flip ${platform_id}-Strand-Flip.txt --make-bed --out temp2 --threads ${task.cpus}
    else
        ln -s temp1.bed temp2.bed; ln -s temp1.bim temp2.bim; ln -s temp1.fam temp2.fam
    fi
    
    # Step 3: Force allele
    if [[ -f "${platform_id}-Force-Allele1.txt" && -s "${platform_id}-Force-Allele1.txt" ]]; then
        plink --bfile temp2 --a1-allele ${platform_id}-Force-Allele1.txt --make-bed --out ${platform_id}_fixed --threads ${task.cpus}
    else
        mv temp2.bed ${platform_id}_fixed.bed; mv temp2.bim ${platform_id}_fixed.bim; mv temp2.fam ${platform_id}_fixed.fam
    fi
    
    echo "Strand fixes applied: ${platform_id}" > ${platform_id}_fixes.log
    rm -f temp*.{bed,bim,fam,log,nosex}
    """
}

// ===================================================================
// PROCESS 1.5: SPLIT BY CHROMOSOME
// ===================================================================

process SPLIT_BY_CHROMOSOME {
    container 'docker://quay.io/biocontainers/plink2:2.00a3.7--h4ac6f70_0'
    
    label 'preqc'
    tag "${platform_id}_chr${chr}"
    
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
    
    plink2 \\
        --bfile \${input_prefix} \\
        --chr ${chr} \\
        --make-bed \\
        --out ${platform_id}_chr${chr} \\
        --threads ${task.cpus}
    """
}

// ===================================================================
// PROCESS 1.6: CREATE IMPUTATION-READY VCFS
// ===================================================================

process CREATE_IMPUTATION_VCFS {
    container 'docker://quay.io/biocontainers/bcftools:1.18--h8b25389_0'
    
    label 'preqc'
    tag "${platform_id}_chr${chr}"
    
    publishDir "${params.outdir}/${platform_id}/05_imputation_ready",
        mode: 'copy',
        pattern: "*.vcf.gz*"
    
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
          path("${platform_id}_chr${chr}.vcf.gz.csi"),
          val(build),
          val(batch), emit: imputation_vcfs
    
    path("${platform_id}_chr${chr}_vcf.log"), emit: vcf_log
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    input_prefix=\$(basename ${bed} .bed)
    n_snps=\$(wc -l < ${bim})
    
    if [[ \${n_snps} -eq 0 ]]; then
        # Empty VCF
        cat > ${platform_id}_chr${chr}.vcf << 'EOFVCF'
##fileformat=VCFv4.2
##contig=<ID=chr${chr},assembly=GRCh38>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
EOFVCF
        bgzip ${platform_id}_chr${chr}.vcf
        bcftools index -c ${platform_id}_chr${chr}.vcf.gz
        echo "Empty VCF" > ${platform_id}_chr${chr}_vcf.log
        exit 0
    fi
    
    # Convert to VCF with chr prefix
    plink2 \\
        --bfile \${input_prefix} \\
        --export vcf-4.2 bgz \\
        --output-chr chr26 \\
        --out temp \\
        --threads ${task.cpus}
    
    # Clean headers and keep only GT
    bcftools annotate \\
        --set-id '%CHROM:%POS:%REF:%ALT' \\
        temp.vcf.gz | \\
    bcftools annotate \\
        -x INFO,^FORMAT/GT \\
        -O z \\
        -o annotated.vcf.gz
    
    # Sort
    bcftools sort \\
        annotated.vcf.gz \\
        -Oz \\
        -o ${platform_id}_chr${chr}.vcf.gz
    
    # Index with CSI
    bcftools index -c ${platform_id}_chr${chr}.vcf.gz
    
    # Log
    vcf_snps=\$(bcftools view -H ${platform_id}_chr${chr}.vcf.gz | wc -l)
    cat > ${platform_id}_chr${chr}_vcf.log << EOF
VCF Created: ${platform_id}_chr${chr}
SNPs: \${vcf_snps}
Format: VCFv4.2, bgzip, CSI indexed
Ready for: TOPMed & All of Us AnVIL
EOF
    
    rm -f temp.vcf.gz annotated.vcf.gz
    """
}

// ===================================================================
// MODULE 1 WORKFLOW
// ===================================================================

workflow MODULE1_PREQC {
    take:
    sample_sheet_ch  // Can contain multiple files per platform
    
    main:
    // Process 1.0: Merge multiple sample files per platform (NEW!)
    MERGE_PLATFORM_SAMPLES(sample_sheet_ch)
    
    // Process 1.1: Convert to standard PLINK format
    CONVERT_TO_PLINK(MERGE_PLATFORM_SAMPLES.out.merged_files)
    
    // Process 1.2: Basic QC
    BASIC_QC(CONVERT_TO_PLINK.out.plink_files)
    
    // Process 1.3: Liftover for hg19
    hg19_ch = BASIC_QC.out.qc_plink.filter { it[4] == 'hg19' }
    hg38_ch = BASIC_QC.out.qc_plink.filter { it[4] == 'hg38' }
    
    LIFTOVER_CROSSMAP(
        hg19_ch,
        file(params.crossmap_chain),
        file(params.hg38_fasta)
    )
    
    all_hg38_ch = LIFTOVER_CROSSMAP.out.lifted_plink.mix(hg38_ch)
    
    // Process 1.4: Strand checking
    STRAND_CHECK_TOPMED(
        all_hg38_ch,
        file("${projectDir}/bin/check-strand-topmed.pl"),
        file(params.topmed_freeze10_ref)
    )
    
    // Process 1.5: Apply strand fixes
    APPLY_STRAND_FIXES(STRAND_CHECK_TOPMED.out.strand_files)
    
    // Process 1.6: Split by chromosome
    chromosomes = Channel.of(1..22)
    chr_input = APPLY_STRAND_FIXES.out.fixed_plink.combine(chromosomes)
    SPLIT_BY_CHROMOSOME(chr_input)
    
    // Process 1.7: Create VCFs
    CREATE_IMPUTATION_VCFS(SPLIT_BY_CHROMOSOME.out.chr_plink)
    
    emit:
    imputation_vcfs = CREATE_IMPUTATION_VCFS.out.imputation_vcfs
    merge_reports = MERGE_PLATFORM_SAMPLES.out.merge_log
    qc_reports = BASIC_QC.out.qc_report
    liftover_reports = LIFTOVER_CROSSMAP.out.liftover_log
    strand_reports = STRAND_CHECK_TOPMED.out.strand_log
    fixes_reports = APPLY_STRAND_FIXES.out.fixes_log
    vcf_logs = CREATE_IMPUTATION_VCFS.out.vcf_log
}

/*
 * ===================================================================
 * SAMPLE SHEET FORMATS:
 * ===================================================================
 * 
 * FORMAT 1: Single file per platform
 * platform_id,file_path,file_type,build,batch
 * Platform1,/data/plat1,plink,hg19,batch1
 * Platform2,/data/plat2,plink,hg38,batch1
 * 
 * FORMAT 2: Multiple files per platform (NEW!)
 * platform_id,file_paths,file_type,build,batch
 * Platform1,"/data/p1_batch1;/data/p1_batch2;/data/p1_batch3",plink,hg19,batch1
 * Platform2,"/data/p2_cohort1;/data/p2_cohort2",plink,hg38,batch1
 * 
 * Note: Multiple files are separated by semicolons
 * 
 * ===================================================================
 * REQUIRED PARAMETERS (nextflow.config):
 * ===================================================================
 * params {
 *     crossmap_chain = "${projectDir}/resources/references/hg19ToHg38.over.chain.gz"
 *     hg38_fasta = "${projectDir}/resources/references/hg38.fa"
 *     topmed_freeze10_ref = "${projectDir}/resources/references/PASS.Variants.TOPMed_freeze10_hg38.tab.gz"
 *     preqc_geno = 0.05
 *     preqc_mind = 0.05
 *     outdir = "results"
 * }
 * ===================================================================
 */
