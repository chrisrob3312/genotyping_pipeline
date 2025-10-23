#!/usr/bin/env nextflow

/*
================================================================================
    Complete Genotyping Imputation and Ancestry Pipeline
================================================================================
    Seven-module pipeline for:
    Module 1: Pre-imputation QC and file preparation
    Module 2: Imputation (TOPMed & All of Us) - parallelized
    Module 3: Post-imputation QC with MagicalRsq
    Module 4: Platform merging with QC
    Module 5: Re-imputation of merged data
    Module 6: Post-merge QC (2nd pass) with GENESIS/PLINK, HWE, Het
    Module 7: Ancestry estimation (GRAF-pop, ADMIXTURE, RFMix)
================================================================================
*/

nextflow.enable.dsl=2

// ============================================================================
// INCLUDE ALL MODULES
// ============================================================================

include { MODULE1_PREIMPUTATION } from './modules/Module1_PreImputation.nf'
include { MODULE2_IMPUTATION } from './modules/Module2_Imputation.nf'
include { MODULE3_POSTQC } from './modules/Module3_PostQC.nf'
include { MODULE4_MERGING } from './modules/Module4_Merging.nf'
include { MODULE5_REIMPUTATION } from './modules/Modules567_and_Main.nf'
include { MODULE6_POSTMERGE_QC2 } from './modules/Modules567_and_Main.nf'
include { MODULE7_ANCESTRY } from './modules/Modules567_and_Main.nf'

// ============================================================================
// HELP MESSAGE
// ============================================================================

def helpMessage() {
    log.info"""
    ================================================================================
    Complete Genotyping Imputation & Ancestry Pipeline
    ================================================================================
    
    A comprehensive 7-module pipeline for genotyping data QC, imputation,
    merging, and ancestry estimation across multiple array platforms.
    
    Usage:
        nextflow run main.nf --sample_sheet samples.tsv [options]
    
    Required Arguments:
        --sample_sheet PATH       TSV with: platform,bed,bim,fam
    
    Optional Arguments - General:
        --input_build STR         Genome build [default: hg19]
        --outdir PATH             Output directory [default: results]
        --skip_modules LIST       Comma-separated modules to skip (e.g., "5,7")
    
    Optional Arguments - Module 1:
        --reference_fasta PATH    Human reference genome
        --liftover_chain PATH     hg19->hg38 chain file
    
    Optional Arguments - Modules 2 & 5 (Imputation):
        --run_topmed             Run TOPMed imputation [default: true]
        --run_anvil              Run All of Us AnVIL [default: false]
        --topmed_api_token STR   TOPMed API token (required if run_topmed)
        --topmed_password STR    TOPMed download password
        --anvil_workspace STR    AnVIL Terra workspace
        --anvil_project STR      AnVIL Google project
    
    Optional Arguments - Module 6:
        --use_genesis            Use GENESIS PCRelate [default: true]
                                 (false = use PLINK relatedness)
    
    Optional Arguments - Module 7:
        --reference_panel PATH   Ancestry reference panel
        --genetic_map PATH       Genetic map directory for RFMix
        --rfmix_script PATH      RFMix execution script
        --admixture_k LIST       K values for ADMIXTURE [default: "5,6,7"]
    
    Execution Profiles:
        -profile standard        Local execution
        -profile docker          Use Docker containers
        -profile singularity     Use Singularity containers
        -profile slurm           SLURM cluster
        -profile awsbatch        AWS Batch
        -profile google          Google Cloud
    
    Module Flow:
        Module 1 → Pre-imputation QC (parallelized per platform)
        Module 2 → Imputation (parallelized per platform × service)
        Module 3 → Post-imputation QC (parallelized per platform)
        Module 4 → Merge platforms
        Module 5 → Re-impute merged data
        Module 6 → Post-merge QC (2nd pass)
        Module 7 → Ancestry estimation (parallelized: GRAF/ADMIXTURE/RFMix)
    
    Examples:
        # Complete pipeline with TOPMed only
        nextflow run main.nf \\
            --sample_sheet platforms.tsv \\
            --run_topmed \\
            --topmed_api_token \$TOKEN \\
            -profile docker
        
        # Skip re-imputation and ancestry
        nextflow run main.nf \\
            --sample_sheet platforms.tsv \\
            --skip_modules "5,7" \\
            -profile slurm,singularity
        
        # Full pipeline with both services
        nextflow run main.nf \\
            --sample_sheet platforms.tsv \\
            --run_topmed --run_anvil \\
            --topmed_api_token \$TOKEN \\
            --anvil_workspace my-workspace \\
            -profile docker
    
    Output Structure:
        results/
        ├── module1/  # Pre-imputation QC
        ├── module2/  # Imputation results
        ├── module3/  # Post-imputation QC
        ├── module4/  # Merged platforms
        ├── module5/  # Re-imputation
        ├── module6/  # Final QC
        └── module7/  # Ancestry results
    
    For detailed documentation:
        See README.md and docs/ directory
    ================================================================================
    """
}

// ============================================================================
// PARAMETER VALIDATION
// ============================================================================

if (params.help) {
    helpMessage()
    exit 0
}

// Check required parameters
if (!params.sample_sheet) {
    log.error "ERROR: --sample_sheet required"
    helpMessage()
    exit 1
}

if (!file(params.sample_sheet).exists()) {
    log.error "ERROR: Sample sheet not found: ${params.sample_sheet}"
    exit 1
}

// Parse skip_modules
def skip_modules = params.skip_modules ? 
                   params.skip_modules.split(',').collect { it.toInteger() } : []

// Validate imputation settings
if (params.run_topmed && !params.topmed_api_token) {
    log.error "ERROR: --topmed_api_token required when run_topmed=true"
    exit 1
}

// ============================================================================
// PRINT PIPELINE INFO
// ============================================================================

log.info """
================================================================================
Complete Genotyping Imputation & Ancestry Pipeline
================================================================================
Sample sheet      : ${params.sample_sheet}
Input build       : ${params.input_build}
Output directory  : ${params.outdir}
Profile           : ${workflow.profile}

Imputation:
  TOPMed          : ${params.run_topmed}
  All of Us AnVIL : ${params.run_anvil}

Modules to skip   : ${skip_modules ?: 'None'}

Modules enabled:
  Module 1 (Pre-QC)       : ${!skip_modules.contains(1)}
  Module 2 (Imputation)   : ${!skip_modules.contains(2)}
  Module 3 (Post-QC)      : ${!skip_modules.contains(3)}
  Module 4 (Merging)      : ${!skip_modules.contains(4)}
  Module 5 (Re-imputation): ${!skip_modules.contains(5)}
  Module 6 (Final QC)     : ${!skip_modules.contains(6)}
  Module 7 (Ancestry)     : ${!skip_modules.contains(7)}
================================================================================
"""

// ============================================================================
// PARSE INPUT
// ============================================================================

Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header: false, sep: '\t')
    .map { row -> 
        def platform = row[0]
        def bed = file(row[1])
        def bim = file(row[2])
        def fam = file(row[3])
        
        // Validate files exist
        if (!bed.exists() || !bim.exists() || !fam.exists()) {
            log.error "ERROR: Files not found for ${platform}"
            exit 1
        }
        
        return tuple(platform, bed, bim, fam)
    }
    .set { input_plink_ch }

// ============================================================================
// MAIN WORKFLOW
// ============================================================================

workflow {
    log.info "Starting Complete Imputation & Ancestry Pipeline"
    log.info "Detected platforms: ${input_plink_ch.count().val}"
    
    // MODULE 1: Pre-imputation QC and file preparation
    if (!skip_modules.contains(1)) {
        log.info "=== Running Module 1: Pre-Imputation QC ==="
        MODULE1_PREIMPUTATION(input_plink_ch)
        module1_out = MODULE1_PREIMPUTATION.out
    } else {
        log.info "Skipping Module 1"
        // Would need to provide input from previous run
        exit 1, "Cannot skip Module 1 without existing outputs"
    }
    
    // MODULE 2: Imputation submission and retrieval
    if (!skip_modules.contains(2)) {
        log.info "=== Running Module 2: Imputation ==="
        MODULE2_IMPUTATION(
            module1_out.vcf_files,
            module1_out.manifests
        )
        module2_out = MODULE2_IMPUTATION.out
    } else {
        log.info "Skipping Module 2"
        exit 1, "Cannot skip Module 2 without existing outputs"
    }
    
    // MODULE 3: Post-imputation QC
    if (!skip_modules.contains(3)) {
        log.info "=== Running Module 3: Post-Imputation QC ==="
        MODULE3_POSTQC(module2_out.imputed_data)
        module3_out = MODULE3_POSTQC.out
    } else {
        log.info "Skipping Module 3"
        exit 1, "Cannot skip Module 3 without existing outputs"
    }
    
    // MODULE 4: Platform merging
    if (!skip_modules.contains(4)) {
        log.info "=== Running Module 4: Platform Merging ==="
        MODULE4_MERGING(module3_out.qc_data)
        module4_out = MODULE4_MERGING.out
    } else {
        log.info "Skipping Module 4"
        exit 1, "Cannot skip Module 4 without existing outputs"
    }
    
    // MODULE 5: Re-imputation of merged data
    if (!skip_modules.contains(5)) {
        log.info "=== Running Module 5: Re-Imputation ==="
        MODULE5_REIMPUTATION(module4_out.merged_data)
        module5_out = MODULE5_REIMPUTATION.out
    } else {
        log.info "Skipping Module 5 - Using Module 4 output"
        // Skip to Module 6 with Module 4 data
        module5_out = [reimputed_data: module4_out.merged_data]
    }
    
    // MODULE 6: Post-merge QC (2nd pass)
    if (!skip_modules.contains(6)) {
        log.info "=== Running Module 6: Final QC ==="
        MODULE6_POSTMERGE_QC2(module5_out.reimputed_data)
        module6_out = MODULE6_POSTMERGE_QC2.out
    } else {
        log.info "Skipping Module 6"
        module6_out = module5_out
    }
    
    // MODULE 7: Ancestry estimation
    if (!skip_modules.contains(7)) {
        log.info "=== Running Module 7: Ancestry Estimation ==="
        
        // Parse ADMIXTURE K values
        def k_values = params.admixture_k.split(',').collect { it.toInteger() }
        
        MODULE7_ANCESTRY(module6_out.final_data_maf)
        module7_out = MODULE7_ANCESTRY.out
    } else {
        log.info "Skipping Module 7"
    }
    
    log.info "=== Pipeline Complete ==="
}

// ============================================================================
// WORKFLOW COMPLETION
// ============================================================================

workflow.onComplete {
    def duration = workflow.duration
    def success_msg = workflow.success ? "SUCCESS" : "FAILED"
    
    log.info """
    ================================================================================
    Pipeline Execution Summary
    ================================================================================
    Status        : ${success_msg}
    Duration      : ${duration}
    Completed at  : ${workflow.complete}
    Work directory: ${workflow.workDir}
    Exit status   : ${workflow.exitStatus}
    ================================================================================
    """
    
    if (workflow.success) {
        log.info """
        Pipeline completed successfully!
        
        Output structure:
        ${params.outdir}/
        ├── module1/  Platform-specific pre-QC
        ├── module2/  Imputation results per platform
        ├── module3/  Post-imputation QC per platform
        ├── module4/  Merged platforms
        ├── module5/  Re-imputed merged data
        ├── module6/  Final QC'd datasets
        │   ├── final_maf01.*    (MAF >= 0.01 filtered)
        │   └── final_nomaf.*    (No MAF filter)
        └── module7/  Ancestry results
            ├── 01_grafpop/
            ├── 02_admixture/
            └── 03_rfmix/
        
        Key outputs:
        - Final datasets: ${params.outdir}/module6/06_final_datasets/
        - Ancestry: ${params.outdir}/module7/
        - QC reports: ${params.outdir}/module*/reports/
        
        Next steps:
        1. Review QC reports in each module
        2. Check ancestry results in module7/
        3. Use final datasets for downstream analysis
        """
    } else {
        log.error """
        Pipeline failed!
        
        Troubleshooting:
        1. Check error logs: ${workflow.workDir}
        2. Review reports: ${params.outdir}/reports/
        3. Check specific module logs
        
        Common issues:
        - API token errors: Check TOPMed/AnVIL credentials
        - Memory errors: Increase resources in nextflow.config
        - Missing files: Verify all input paths in sample sheet
        """
    }
}

workflow.onError {
    log.error "Pipeline stopped with error: ${workflow.errorMessage}"
    log.error "Failed at: ${workflow.errorReport}"
}
