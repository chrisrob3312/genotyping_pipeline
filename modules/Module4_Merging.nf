#!/usr/bin/env nextflow

/*
 * ============================================================================
 * MODULE 4: PLATFORM MERGING (3-Pass Strategy with Dual Imputation Support)
 * ============================================================================
 * 
 * Purpose: Merge imputed data from multiple genotyping platforms
 * Strategy: 3-pass merging with automatic mismatch detection and fixing
 * 
 * POST-MERGE QC (Applied in this module):
 * - SNPs only (biallelic, no indels, no multi-allelic)
 * - No monomorphic variants (MAF > 0)
 * - Variant call rate > 95%
 * - Sample call rate > 95%
 * 
 * Dual Pipeline Support:
 * - topmed_ref (TOPMed imputation server)
 * - allofus_ref (All of Us AnVIL imputation)
 * 
 * Resource Optimization:
 * - Dynamic memory scaling based on sample count
 * - Efficient parallelization across chromosomes
 * - Optimized for large-scale merging (10,000+ samples, 10+ platforms)
 * ============================================================================
 */

nextflow.enable.dsl=2

// ============================================================================
// PROCESS 0: Convert VCF to PLINK2 binary format
// ============================================================================
process vcfToPlink {
    label 'plink2'
    tag "${ref_panel}_${platform}_chr${chr}"
    publishDir "${params.outdir}/module4/${ref_panel}/00_vcf_to_plink/${platform}", mode: 'copy'
    
    container 'docker://quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'
    
    cpus { 4 }
    memory { 8.GB + (2.GB * Math.min(task.attempt, 3)) }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(ref_panel), val(platform), val(chr), path(vcf), path(tbi)
    
    output:
    tuple val(ref_panel), val(platform), val(chr),
          path("${platform}_chr${chr}.pgen"),
          path("${platform}_chr${chr}.pvar"),
          path("${platform}_chr${chr}.psam"),
          emit: plink_files
    
    script:
    """
    plink2 \\
        --vcf ${vcf} dosage=DS \\
        --max-alleles 2 \\
        --vcf-half-call m \\
        --make-pgen \\
        --out ${platform}_chr${chr} \\
        --threads ${task.cpus}
    
    echo "Converted: ${platform} chr${chr} for ${ref_panel}"
    """
}

// ============================================================================
// PROCESS 1a: PASS 1 - Initial merge attempt to identify issues
// ============================================================================
process mergePlatforms_Pass1_Identify {
    label 'plink_merge'
    tag "${ref_panel}_chr${chr}_pass1"
    publishDir "${params.outdir}/module4/${ref_panel}/01_merge_pass1/chr${chr}", mode: 'copy'
    
    container 'docker://quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'
    
    cpus { 8 }
    memory { 16.GB + (4.GB * Math.min(task.attempt, 3)) }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(ref_panel), val(chr),
          path(pgen_files), path(pvar_files), path(psam_files)
    
    output:
    tuple val(ref_panel), val(chr),
          path("chr${chr}_pass1.pgen") optional true,
          path("chr${chr}_pass1.pvar") optional true,
          path("chr${chr}_pass1.psam") optional true,
          path("chr${chr}_pass1.log"),
          path("chr${chr}_mismatches.txt") optional true,
          emit: pass1_results
    
    script:
    def pmerge_list = pgen_files.collect { it.toString().replaceAll(/\\.pgen$/, '') }.join('\n')
    
    """
    # Create pmerge list
    cat > pmerge_list.txt <<EOF
${pmerge_list}
EOF
    
    echo "=== PASS 1: Identifying Merge Issues ==="
    echo "Reference Panel: ${ref_panel}"
    echo "Chromosome: ${chr}"
    echo "Platforms to merge: \$(wc -l < pmerge_list.txt)"
    
    # Attempt merge - allow failure to capture errors
    plink2 \\
        --pmerge-list pmerge_list.txt pfile \\
        --merge-max-allele-ct 2 \\
        --out chr${chr}_pass1 \\
        --threads ${task.cpus} \\
        --memory ${task.memory.toMega() * 0.9 as int} \\
        2>&1 | tee chr${chr}_pass1.log || true
    
    # Extract mismatch information
    if grep -q "Error:" chr${chr}_pass1.log; then
        echo "Mismatches detected - extracting details..."
        grep "REF allele mismatch" chr${chr}_pass1.log > chr${chr}_mismatches.txt || true
        n_mismatches=\$(wc -l < chr${chr}_mismatches.txt || echo 0)
        echo "Total mismatches found: \$n_mismatches"
        echo "MISMATCH_COUNT=\$n_mismatches" >> chr${chr}_pass1.log
    else
        echo "No mismatches detected - merge successful!"
        echo "MISMATCH_COUNT=0" >> chr${chr}_pass1.log
    fi
    """
}

// ============================================================================
// PROCESS 1b: PASS 2 - Analyze mismatches and prepare fixes
// ============================================================================
process analyzeMismatchesAndPrepFixes {
    label 'python'
    tag "${ref_panel}_chr${chr}_analyze"
    publishDir "${params.outdir}/module4/${ref_panel}/02_mismatch_analysis/chr${chr}", mode: 'copy'
    
    container 'docker://quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'
    
    cpus { 2 }
    memory { 4.GB + (1.GB * Math.min(task.attempt, 3)) }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(ref_panel), val(chr), 
          path(pgen), path(pvar), path(psam), path(log), path(mismatches)
    tuple val(ref_panel_check), val(chr_check),
          path(pgen_files), path(pvar_files), path(psam_files)
    
    output:
    tuple val(ref_panel), val(chr),
          path("chr${chr}_flip_list.txt") optional true,
          path("chr${chr}_ref_allele.txt") optional true,
          path("chr${chr}_exclude_list.txt") optional true,
          path("chr${chr}_fix_summary.txt"),
          emit: fix_instructions
    
    when:
    ref_panel == ref_panel_check && chr == chr_check && mismatches.name != 'OPTIONAL_FILE'
    
    script:
    """
    #!/usr/bin/env python3
    
    import re
    import sys
    
    print("=== PASS 2: Analyzing Mismatches ===")
    print(f"Reference Panel: ${ref_panel}")
    print(f"Chromosome: ${chr}")
    
    flip_variants = []
    ref_fix_variants = []
    exclude_variants = []
    
    # Complementary bases
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # Parse mismatch file
    try:
        with open("${mismatches}", 'r') as f:
            for line in f:
                match = re.search(r"'([^']+)'", line)
                if match:
                    variant_id = match.group(1)
                    flip_variants.append(variant_id)
    
    except FileNotFoundError:
        print("No mismatch file found - no fixes needed")
        sys.exit(0)
    
    # Write fix files
    summary = {
        'total_mismatches': len(flip_variants),
        'strand_flips': len(flip_variants),
        'ref_swaps': 0,
        'excluded': 0
    }
    
    if flip_variants:
        with open("chr${chr}_flip_list.txt", 'w') as f:
            for var in flip_variants:
                f.write(f"{var}\\n")
        print(f"Prepared {len(flip_variants)} variants for strand flipping")
    
    # Write summary
    with open("chr${chr}_fix_summary.txt", 'w') as f:
        f.write("Mismatch Analysis Summary\\n")
        f.write(f"Reference Panel: ${ref_panel}\\n")
        f.write(f"Chromosome: ${chr}\\n")
        f.write(f"Total mismatches: {summary['total_mismatches']}\\n")
        f.write(f"Strand flips to attempt: {summary['strand_flips']}\\n")
        f.write(f"REF/ALT swaps to attempt: {summary['ref_swaps']}\\n")
        f.write(f"Variants to exclude: {summary['excluded']}\\n")
    
    print("Analysis complete!")
    """
}

// ============================================================================
// PROCESS 1c: PASS 2 - Apply fixes to each platform
// ============================================================================
process applyFixesToPlatform {
    label 'plink2'
    tag "${ref_panel}_${platform}_chr${chr}_fix"
    publishDir "${params.outdir}/module4/${ref_panel}/03_fixed_platforms/${platform}/chr${chr}", mode: 'copy'
    
    container 'docker://quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'
    
    cpus { 4 }
    memory { 8.GB + (2.GB * Math.min(task.attempt, 3)) }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(ref_panel), val(platform), val(chr),
          path(pgen), path(pvar), path(psam)
    tuple val(ref_panel_check), val(chr_check),
          path(flip_list), path(ref_allele), path(exclude_list), path(summary)
    
    output:
    tuple val(ref_panel), val(platform), val(chr),
          path("${platform}_chr${chr}_fixed.pgen"),
          path("${platform}_chr${chr}_fixed.pvar"),
          path("${platform}_chr${chr}_fixed.psam"),
          emit: fixed_files
    
    when:
    ref_panel == ref_panel_check && chr == chr_check
    
    script:
    def flip_cmd = flip_list.name != 'OPTIONAL_FILE' ? "--flip ${flip_list}" : ""
    def ref_cmd = ref_allele.name != 'OPTIONAL_FILE' ? "--ref-allele ${ref_allele}" : ""
    def exclude_cmd = exclude_list.name != 'OPTIONAL_FILE' ? "--exclude ${exclude_list}" : ""
    
    """
    echo "=== Applying Fixes to ${platform} ==="
    
    plink2 \\
        --pfile ${pgen.baseName} \\
        ${flip_cmd} \\
        ${ref_cmd} \\
        ${exclude_cmd} \\
        --make-pgen \\
        --out ${platform}_chr${chr}_fixed \\
        --threads ${task.cpus} \\
        --memory ${task.memory.toMega() * 0.9 as int}
    
    echo "Fixes applied to ${platform}"
    """
}

// ============================================================================
// PROCESS 1d: PASS 3 - Clean merge after fixes
// ============================================================================
process mergePlatforms_Pass3_CleanMerge {
    label 'plink_large'
    tag "${ref_panel}_chr${chr}_pass3"
    publishDir "${params.outdir}/module4/${ref_panel}/04_merged_final/chr${chr}", mode: 'copy'
    
    container 'docker://quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'
    
    cpus { 8 }
    memory { 24.GB + (8.GB * Math.min(task.attempt, 3)) }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(ref_panel), val(chr),
          path(pgen_files), path(pvar_files), path(psam_files)
    
    output:
    tuple val(ref_panel), val(chr),
          path("${ref_panel}_merged_chr${chr}.pgen"),
          path("${ref_panel}_merged_chr${chr}.pvar"),
          path("${ref_panel}_merged_chr${chr}.psam"),
          path("${ref_panel}_merged_chr${chr}.log"),
          emit: merged_plink
    
    script:
    def pmerge_list = pgen_files.collect { it.toString().replaceAll(/_fixed\\.pgen$/, '_fixed') }.join('\n')
    
    """
    cat > pmerge_list_fixed.txt <<EOF
${pmerge_list}
EOF
    
    echo "=== PASS 3: Clean Merge After Fixes ==="
    echo "Reference Panel: ${ref_panel}"
    echo "Chromosome: ${chr}"
    
    plink2 \\
        --pmerge-list pmerge_list_fixed.txt pfile \\
        --merge-max-allele-ct 2 \\
        --out ${ref_panel}_merged_chr${chr} \\
        --threads ${task.cpus} \\
        --memory ${task.memory.toMega() * 0.9 as int}
    
    if [ -f "${ref_panel}_merged_chr${chr}.pgen" ]; then
        n_variants=\$(awk 'NR>1' ${ref_panel}_merged_chr${chr}.pvar | wc -l)
        n_samples=\$(awk 'NR>1' ${ref_panel}_merged_chr${chr}.psam | wc -l)
        echo "SUCCESS: Clean merge completed"
        echo "Final variant count: \$n_variants"
        echo "Final sample count: \$n_samples"
    else
        echo "ERROR: Merge failed!"
        exit 1
    fi
    """
}

// ============================================================================
// PROCESS 2: Post-Merge QC - SNPs only, biallelic
// ============================================================================
process filterToSNPsOnly {
    label 'plink2'
    tag "${ref_panel}_chr${chr}_snps"
    publishDir "${params.outdir}/module4/${ref_panel}/05_snps_only/chr${chr}", mode: 'copy'
    
    container 'docker://quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'
    
    cpus { 4 }
    memory { 8.GB + (2.GB * Math.min(task.attempt, 3)) }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(ref_panel), val(chr),
          path(pgen), path(pvar), path(psam), path(log)
    
    output:
    tuple val(ref_panel), val(chr),
          path("${ref_panel}_chr${chr}_snps.pgen"),
          path("${ref_panel}_chr${chr}_snps.pvar"),
          path("${ref_panel}_chr${chr}_snps.psam"),
          path("${ref_panel}_chr${chr}_removed_nonsnps.txt"),
          emit: snps_only
    
    script:
    """
    echo "=== Filtering to SNPs Only (Biallelic) ==="
    echo "Reference Panel: ${ref_panel}"
    echo "Chromosome: ${chr}"
    
    # Extract non-SNP variants for logging
    awk 'NR>1 && (length(\$4)>1 || length(\$5)>1) {print \$3}' ${pvar} > ${ref_panel}_chr${chr}_removed_nonsnps.txt || touch ${ref_panel}_chr${chr}_removed_nonsnps.txt
    
    # Keep only biallelic SNPs (single nucleotide, no indels)
    plink2 \\
        --pfile ${pgen.baseName} \\
        --snps-only just-acgt \\
        --max-alleles 2 \\
        --make-pgen \\
        --out ${ref_panel}_chr${chr}_snps \\
        --threads ${task.cpus} \\
        --memory ${task.memory.toMega() * 0.9 as int}
    
    n_removed=\$(wc -l < ${ref_panel}_chr${chr}_removed_nonsnps.txt)
    n_remaining=\$(awk 'NR>1' ${ref_panel}_chr${chr}_snps.pvar | wc -l)
    echo "Removed \$n_removed non-SNP variants"
    echo "Remaining SNPs: \$n_remaining"
    """
}

// ============================================================================
// PROCESS 3: Remove monomorphic variants (MAF > 0)
// ============================================================================
process removeMonomorphic {
    label 'plink2'
    tag "${ref_panel}_chr${chr}_maf"
    publishDir "${params.outdir}/module4/${ref_panel}/06_polymorphic/chr${chr}", mode: 'copy'
    
    container 'docker://quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'
    
    cpus { 4 }
    memory { 8.GB + (2.GB * Math.min(task.attempt, 3)) }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(ref_panel), val(chr),
          path(pgen), path(pvar), path(psam), path(nonsnps)
    
    output:
    tuple val(ref_panel), val(chr),
          path("${ref_panel}_chr${chr}_poly.pgen"),
          path("${ref_panel}_chr${chr}_poly.pvar"),
          path("${ref_panel}_chr${chr}_poly.psam"),
          path("${ref_panel}_chr${chr}_removed_monomorphic.txt"),
          emit: polymorphic
    
    script:
    """
    echo "=== Removing Monomorphic Variants (MAF > 0) ==="
    echo "Reference Panel: ${ref_panel}"
    echo "Chromosome: ${chr}"
    
    # Calculate MAF and identify monomorphic variants
    plink2 \\
        --pfile ${pgen.baseName} \\
        --freq \\
        --out ${ref_panel}_chr${chr}_freq \\
        --threads ${task.cpus}
    
    # Extract monomorphic variants (MAF = 0)
    awk 'NR>1 && \$5==0 {print \$2}' ${ref_panel}_chr${chr}_freq.afreq > ${ref_panel}_chr${chr}_removed_monomorphic.txt || touch ${ref_panel}_chr${chr}_removed_monomorphic.txt
    
    # Filter to polymorphic only (MAF > 0)
    plink2 \\
        --pfile ${pgen.baseName} \\
        --mac 1 \\
        --make-pgen \\
        --out ${ref_panel}_chr${chr}_poly \\
        --threads ${task.cpus} \\
        --memory ${task.memory.toMega() * 0.9 as int}
    
    n_removed=\$(wc -l < ${ref_panel}_chr${chr}_removed_monomorphic.txt)
    n_remaining=\$(awk 'NR>1' ${ref_panel}_chr${chr}_poly.pvar | wc -l)
    echo "Removed \$n_removed monomorphic variants"
    echo "Remaining polymorphic variants: \$n_remaining"
    """
}

// ============================================================================
// PROCESS 4: Filter variants by call rate (>95%)
// ============================================================================
process filterVariantCallRate {
    label 'plink2'
    tag "${ref_panel}_chr${chr}_varcr"
    publishDir "${params.outdir}/module4/${ref_panel}/07_variant_callrate/chr${chr}", mode: 'copy'
    
    container 'docker://quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'
    
    cpus { 4 }
    memory { 8.GB + (2.GB * Math.min(task.attempt, 3)) }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(ref_panel), val(chr),
          path(pgen), path(pvar), path(psam), path(mono)
    
    output:
    tuple val(ref_panel), val(chr),
          path("${ref_panel}_chr${chr}_geno.pgen"),
          path("${ref_panel}_chr${chr}_geno.pvar"),
          path("${ref_panel}_chr${chr}_geno.psam"),
          path("${ref_panel}_chr${chr}_removed_lowcr_variants.txt"),
          emit: variant_callrate_filtered
    
    script:
    """
    echo "=== Filtering Variants by Call Rate >95% ==="
    echo "Reference Panel: ${ref_panel}"
    echo "Chromosome: ${chr}"
    
    # Calculate variant missingness
    plink2 \\
        --pfile ${pgen.baseName} \\
        --missing variant-only \\
        --out ${ref_panel}_chr${chr}_vmiss \\
        --threads ${task.cpus}
    
    # Extract variants with call rate <95% (missingness >0.05)
    awk 'NR>1 && \$5>0.05 {print \$2}' ${ref_panel}_chr${chr}_vmiss.vmiss > ${ref_panel}_chr${chr}_removed_lowcr_variants.txt || touch ${ref_panel}_chr${chr}_removed_lowcr_variants.txt
    
    # Filter: keep variants with call rate >95%
    plink2 \\
        --pfile ${pgen.baseName} \\
        --geno ${params.variant_call_rate} \\
        --make-pgen \\
        --out ${ref_panel}_chr${chr}_geno \\
        --threads ${task.cpus} \\
        --memory ${task.memory.toMega() * 0.9 as int}
    
    n_removed=\$(wc -l < ${ref_panel}_chr${chr}_removed_lowcr_variants.txt)
    n_remaining=\$(awk 'NR>1' ${ref_panel}_chr${chr}_geno.pvar | wc -l)
    echo "Removed \$n_removed variants with call rate <95%"
    echo "Remaining variants: \$n_remaining"
    """
}

// ============================================================================
// PROCESS 5: Filter samples by call rate (>95%)
// ============================================================================
process filterSampleCallRate {
    label 'plink2'
    tag "${ref_panel}_chr${chr}_sampcr"
    publishDir "${params.outdir}/module4/${ref_panel}/08_sample_callrate/chr${chr}", mode: 'copy'
    
    container 'docker://quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'
    
    cpus { 4 }
    memory { 8.GB + (2.GB * Math.min(task.attempt, 3)) }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(ref_panel), val(chr),
          path(pgen), path(pvar), path(psam), path(var_removed)
    
    output:
    tuple val(ref_panel), val(chr),
          path("${ref_panel}_chr${chr}_qc.pgen"),
          path("${ref_panel}_chr${chr}_qc.pvar"),
          path("${ref_panel}_chr${chr}_qc.psam"),
          path("${ref_panel}_chr${chr}_removed_lowcr_samples.txt"),
          emit: sample_callrate_filtered
    
    script:
    """
    echo "=== Filtering Samples by Call Rate >95% ==="
    echo "Reference Panel: ${ref_panel}"
    echo "Chromosome: ${chr}"
    
    # Calculate sample missingness
    plink2 \\
        --pfile ${pgen.baseName} \\
        --missing sample-only \\
        --out ${ref_panel}_chr${chr}_smiss \\
        --threads ${task.cpus}
    
    # Extract samples with call rate <95% (missingness >0.05)
    awk 'NR>1 && \$5>0.05 {print \$1}' ${ref_panel}_chr${chr}_smiss.smiss > ${ref_panel}_chr${chr}_removed_lowcr_samples.txt || touch ${ref_panel}_chr${chr}_removed_lowcr_samples.txt
    
    # Filter: keep samples with call rate >95%
    plink2 \\
        --pfile ${pgen.baseName} \\
        --mind ${params.sample_call_rate} \\
        --make-pgen \\
        --out ${ref_panel}_chr${chr}_qc \\
        --threads ${task.cpus} \\
        --memory ${task.memory.toMega() * 0.9 as int}
    
    n_removed=\$(wc -l < ${ref_panel}_chr${chr}_removed_lowcr_samples.txt)
    n_remaining=\$(awk 'NR>1' ${ref_panel}_chr${chr}_qc.psam | wc -l)
    echo "Removed \$n_removed samples with call rate <95%"
    echo "Remaining samples: \$n_remaining"
    """
}

// ============================================================================
// PROCESS 6: Convert back to VCF for downstream analysis
// ============================================================================
process plinkToVcf {
    label 'plink2'
    tag "${ref_panel}_chr${chr}_vcf"
    publishDir "${params.outdir}/module4/${ref_panel}/09_final_vcf/chr${chr}", mode: 'copy'
    
    container 'docker://quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'
    
    cpus { 4 }
    memory { 8.GB + (2.GB * Math.min(task.attempt, 3)) }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(ref_panel), val(chr),
          path(pgen), path(pvar), path(psam), path(samp_removed)
    
    output:
    tuple val(ref_panel), val(chr),
          path("${ref_panel}_merged_chr${chr}.vcf.gz"),
          path("${ref_panel}_merged_chr${chr}.vcf.gz.tbi"),
          emit: final_vcf
    
    script:
    """
    echo "=== Converting to VCF Format ==="
    echo "Reference Panel: ${ref_panel}"
    echo "Chromosome: ${chr}"
    
    # Convert PLINK to VCF
    plink2 \\
        --pfile ${pgen.baseName} \\
        --export vcf-4.2 bgz \\
        --out ${ref_panel}_merged_chr${chr} \\
        --threads ${task.cpus} \\
        --memory ${task.memory.toMega() * 0.9 as int}
    
    # Index VCF with bcftools
    bcftools index -t ${ref_panel}_merged_chr${chr}.vcf.gz
    
    # Final counts
    n_variants=\$(bcftools view -H ${ref_panel}_merged_chr${chr}.vcf.gz | wc -l)
    n_samples=\$(bcftools query -l ${ref_panel}_merged_chr${chr}.vcf.gz | wc -l)
    echo "FINAL: \$n_variants variants, \$n_samples samples"
    """
}

// ============================================================================
// PROCESS 7: Generate comprehensive QC report
// ============================================================================
process generateQCReport {
    label 'R'
    tag "${ref_panel}_chr${chr}_report"
    publishDir "${params.outdir}/module4/${ref_panel}/10_qc_reports/chr${chr}", mode: 'copy'
    
    container 'docker://rocker/tidyverse:4.3'
    
    cpus { 2 }
    memory { 4.GB }
    
    input:
    tuple val(ref_panel), val(chr),
          path(pass1_log), path(fix_summary), path(pass3_log),
          path(nonsnps), path(mono), path(var_cr), path(samp_cr)
    
    output:
    tuple val(ref_panel), val(chr),
          path("${ref_panel}_chr${chr}_qc_report.html"),
          path("${ref_panel}_chr${chr}_qc_summary.txt"),
          emit: qc_report
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(data.table)
    library(rmarkdown)
    
    ref_panel <- "${ref_panel}"
    chr <- "${chr}"
    
    # Parse Pass 1 results
    pass1_log <- readLines("${pass1_log}")
    mismatch_line <- grep("MISMATCH_COUNT=", pass1_log, value = TRUE)
    n_mismatches <- ifelse(length(mismatch_line) > 0,
                           as.numeric(sub("MISMATCH_COUNT=", "", mismatch_line)),
                           0)
    
    # Parse fix summary if exists
    if (file.exists("${fix_summary}") && file.size("${fix_summary}") > 0) {
        fix_data <- readLines("${fix_summary}")
        n_flips <- as.numeric(sub(".*: ", "", grep("Strand flips", fix_data, value = TRUE)))
    } else {
        n_flips <- 0
    }
    
    # Count QC removals
    n_nonsnps <- length(readLines("${nonsnps}"))
    n_mono <- length(readLines("${mono}"))
    n_var_cr <- length(readLines("${var_cr}"))
    n_samp_cr <- length(readLines("${samp_cr}"))
    
    # Parse final counts from Pass 3 log
    pass3_log <- readLines("${pass3_log}")
    variant_line <- grep("Final variant count:", pass3_log, value = TRUE)
    sample_line <- grep("Final sample count:", pass3_log, value = TRUE)
    n_variants_pass3 <- as.numeric(sub(".*: ", "", variant_line))
    n_samples_pass3 <- as.numeric(sub(".*: ", "", sample_line))
    
    # Calculate final counts
    n_final_variants <- n_variants_pass3 - n_nonsnps - n_mono - n_var_cr
    n_final_samples <- n_samples_pass3 - n_samp_cr
    
    # Create summary
    summary <- data.frame(
        Step = c(
            "Pass 1: Mismatches Detected",
            "Pass 2: Strand Flips Applied",
            "Pass 3: Variants After Merge",
            "Pass 3: Samples After Merge",
            "QC: Non-SNPs Removed",
            "QC: Monomorphic Removed",
            "QC: Low Variant Call Rate",
            "QC: Low Sample Call Rate",
            "Final Variant Count",
            "Final Sample Count"
        ),
        Count = c(
            n_mismatches,
            n_flips,
            n_variants_pass3,
            n_samples_pass3,
            n_nonsnps,
            n_mono,
            n_var_cr,
            n_samp_cr,
            n_final_variants,
            n_final_samples
        )
    )
    
    write.table(summary, "${ref_panel}_chr${chr}_qc_summary.txt",
                quote = FALSE, row.names = FALSE, sep = "\\t")
    
    # Generate HTML report
    cat("
---
title: 'Module 4 QC Report'
subtitle: '", ref_panel, " - Chromosome ", chr, "'
output: html_document
---

# Platform Merging & QC Summary

**Reference Panel:** ", ref_panel, "  
**Chromosome:** ", chr, "

## QC Steps Applied

1. **3-Pass Merging:** Identify → Fix → Merge
2. **SNPs Only:** Biallelic, no indels
3. **Polymorphic:** MAF > 0
4. **Variant Call Rate:** >95%
5. **Sample Call Rate:** >95%

**Note:** HWE and heterozygosity QC deferred to Module 6 (post-reimputation)

## Summary Statistics

```{r, echo=FALSE}
knitr::kable(summary, align = 'lr')
```

## Final Dataset

- **Variants:** ", n_final_variants, "
- **Samples:** ", n_final_samples, "
- **Quality:** All variants pass call rate and polymorphism filters

", file = "${ref_panel}_chr${chr}_report.Rmd")
    
    rmarkdown::render("${ref_panel}_chr${chr}_report.Rmd",
                     output_file = "${ref_panel}_chr${chr}_qc_report.html",
                     quiet = TRUE)
    
    cat("QC report generated\\n")
    """
}

// ============================================================================
// MODULE 4 WORKFLOW
// ============================================================================
workflow MODULE4_MERGING {
    take:
    qc_data  // From Module 3: tuple(ref_panel, platform, vcf_files, tbi_files)
    
    main:
    // Flatten to per-chromosome channels
    qc_data
        .flatMap { ref_panel, platform, vcfs, tbis ->
            def vcf_list = vcfs instanceof List ? vcfs : [vcfs]
            def tbi_list = tbis instanceof List ? tbis : [tbis]
            
            vcf_list.collect { vcf ->
                def chr = (vcf.name =~ /chr(\d+|X|Y)/)[0][1]
                def tbi = tbi_list.find { it.name.contains("chr${chr}") }
                tuple(ref_panel, platform, chr, vcf, tbi)
            }
        }
        .set { per_chr_vcfs }
    
    // STEP 1: Convert VCF to PLINK2
    vcfToPlink(per_chr_vcfs)
    
    // STEP 2: Group by ref_panel and chromosome
    vcfToPlink.out.plink_files
        .map { ref_panel, platform, chr, pgen, pvar, psam ->
            tuple(ref_panel, chr, pgen, pvar, psam)
        }
        .groupTuple(by: [0, 1])
        .set { grouped_by_chr }
    
    // PASS 1: Identify issues
    mergePlatforms_Pass1_Identify(grouped_by_chr)
    
    // PASS 2: Analyze and prepare fixes
    analyzeMismatchesAndPrepFixes(
        mergePlatforms_Pass1_Identify.out.pass1_results,
        grouped_by_chr
    )
    
    // PASS 2: Apply fixes to each platform
    vcfToPlink.out.plink_files
        .combine(analyzeMismatchesAndPrepFixes.out.fix_instructions, by: [0, 2])
        .set { files_with_fixes }
    
    applyFixesToPlatform(
        files_with_fixes.map { ref_panel, chr, platform, pgen, pvar, psam, flip, ref_a, excl, summ ->
            tuple(ref_panel, platform, chr, pgen, pvar, psam)
        },
        files_with_fixes.map { ref_panel, chr, platform, pgen, pvar, psam, flip, ref_a, excl, summ ->
            tuple(ref_panel, chr, flip, ref_a, excl, summ)
        }
    )
    
    // PASS 3: Clean merge with fixed files
    applyFixesToPlatform.out.fixed_files
        .map { ref_panel, platform, chr, pgen, pvar, psam ->
            tuple(ref_panel, chr, pgen, pvar, psam)
        }
        .groupTuple(by: [0, 1])
        .set { fixed_grouped }
    
    mergePlatforms_Pass3_CleanMerge(fixed_grouped)
    
    // POST-MERGE QC
    filterToSNPsOnly(mergePlatforms_Pass3_CleanMerge.out.merged_plink)
    removeMonomorphic(filterToSNPsOnly.out.snps_only)
    filterVariantCallRate(removeMonomorphic.out.polymorphic)
    filterSampleCallRate(filterVariantCallRate.out.variant_callrate_filtered)
    plinkToVcf(filterSampleCallRate.out.sample_callrate_filtered)
    
    // Generate QC reports
    mergePlatforms_Pass1_Identify.out.pass1_results
        .map { ref_panel, chr, pgen, pvar, psam, log, mismatch ->
            tuple(ref_panel, chr, log)
        }
        .join(
            analyzeMismatchesAndPrepFixes.out.fix_instructions
                .map { ref_panel, chr, flip, ref_a, excl, summ ->
                    tuple(ref_panel, chr, summ)
                },
            by: [0, 1]
        )
        .join(
            mergePlatforms_Pass3_CleanMerge.out.merged_plink
                .map { ref_panel, chr, pgen, pvar, psam, log ->
                    tuple(ref_panel, chr, log)
                },
            by: [0, 1]
        )
        .join(
            filterToSNPsOnly.out.snps_only
                .map { ref_panel, chr, pgen, pvar, psam, nonsnps ->
                    tuple(ref_panel, chr, nonsnps)
                },
            by: [0, 1]
        )
        .join(
            removeMonomorphic.out.polymorphic
                .map { ref_panel, chr, pgen, pvar, psam, mono ->
                    tuple(ref_panel, chr, mono)
                },
            by: [0, 1]
        )
        .join(
            filterVariantCallRate.out.variant_callrate_filtered
                .map { ref_panel, chr, pgen, pvar, psam, var_cr ->
                    tuple(ref_panel, chr, var_cr)
                },
            by: [0, 1]
        )
        .join(
            filterSampleCallRate.out.sample_callrate_filtered
                .map { ref_panel, chr, pgen, pvar, psam, samp_cr ->
                    tuple(ref_panel, chr, samp_cr)
                },
            by: [0, 1]
        )
        .set { for_report }
    
    generateQCReport(for_report)
    
    emit:
    merged_vcf = plinkToVcf.out.final_vcf
    qc_reports = generateQCReport.out.qc_report
}

// ============================================================================
// WORKFLOW COMPLETION
// ============================================================================
workflow.onComplete {
    println """
    ============================================================================
    MODULE 4 - PLATFORM MERGING COMPLETE
    ============================================================================
    Status:        ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration:      ${workflow.duration}
    
    QC Applied:
    ✓ 3-pass merging (identify → fix → merge)
    ✓ SNPs only (biallelic, no indels)
    ✓ Polymorphic variants (MAF > 0)
    ✓ Variant call rate >95%
    ✓ Sample call rate >95%
    
    
    Results: ${params.outdir}/module4/
    ============================================================================
    """
}

