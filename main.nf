#!/usr/bin/env nextflow

/*
================================================================================
    Complete Genotyping Imputation and Ancestry Pipeline
================================================================================
    Seven-module pipeline for:
    Module 1: Pre-imputation QC and file preparation
    Module 2: Imputation (TOPMed & All of Us) - parallelized
    Module 3: Post-imputation QC with MagicalRsq-X
    Module 4: Platform merging with QC
    Module 5: Re-imputation of merged data
    Module 6: Post-merge QC (2nd pass) with GENESIS/PLINK, HWE, Het
    Module 7: Ancestry estimation (GRAF-pop, ADMIXTURE, RFMix)
    
    Author: Your Name
    Version: 1.0.0
    Updated: October 2025
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
include { MODULE5_REIMPUTATION } from './modules/Module5_Reimputation.nf'
include { MODULE6_POSTMERGE_QC2 } from './modules/Module6_PostMergeQC.nf'
include { MODULE7_ANCESTRY } from './modules/Module7_Ancestry.nf'

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
    
    USAGE:
        nextflow run main.nf --sample_sheet <file> [options]
    
    REQUIRED ARGUMENTS:
        --sample_sheet PATH      Sample sheet with platform info (TSV format)
                                 Columns: platform_id, bed, bim, fam
    
    OPTIONAL ARGUMENTS:
        --outdir PATH            Output directory [default: results]
        --input_build STR        Genome build: hg19 or hg38 [default: hg19]
        
    MODULE CONTROL:
        --skip_module2           Skip imputation (Module 2)
        --skip_module3           Skip post-imputation QC (Module 3)
        --skip_module4           Skip merging (Module 4)
        --skip_module5           Skip re-imputation (Module 5)
        --skip_module6           Skip final QC (Module 6)
        --skip_module7           Skip ancestry estimation (Module 7)
        
    IMPUTATION OPTIONS (Module 2 & 5):
        --run_topmed             Run TOPMed imputation
        --topmed_api_token STR   TOPMed API authentication token
        --topmed_password STR    Password for encrypted results
        
        --run_anvil              Run All of Us AnVIL imputation
        --anvil_workspace STR    Terra workspace name
        --anvil_billing_project  Google Cloud billing project
    
    QC OPTIONS (Module 3 & 6):
        --magicalrsq_threshold   MagicalRsq-X R² threshold [default: 0.3]
        --use_genesis            Use GENESIS for relatedness [default: true]
        --callrate_threshold     Call rate threshold [default: 0.95]
        --hwe_threshold          Hardy-Weinberg p-value [default: 1e-6]
        
    ANCESTRY OPTIONS (Module 7):
        --admixture_k STR        ADMIXTURE K values [default: "5,6,7"]
        --run_tractor_lai        Run TractorWorkflow LAI tools [default: false]
    
    REFERENCE FILES:
        --rayner_reference PATH  Will Rayner reference for strand checking
        --liftover_chain PATH    LiftOver chain file (hg19ToHg38)
        --grafpop_reference PATH GRAF-pop reference panel
        --rfmix_reference_dir    RFMix reference directory
    
    EXECUTION PROFILES:
        -profile local           Run locally with conda
        -profile slurm           Run on SLURM cluster
        -profile pbs             Run on PBS cluster
        -profile docker          Use Docker containers
        -profile singularity     Use Singularity containers
        -profile apptainer       Use Apptainer containers
        -profile awsbatch        Run on AWS Batch
        -profile google          Run on Google Cloud
    
    EXAMPLES:
    
    1. Full pipeline with TOPMed imputation:
       nextflow run main.nf \\
           --sample_sheet platforms.tsv \\
           --run_topmed \\
           --topmed_api_token \$TOKEN \\
           -profile slurm,apptainer
    
    2. QC and ancestry only (skip imputation):
       nextflow run main.nf \\
           --sample_sheet platforms.tsv \\
           --skip_module2 --skip_module5 \\
           -profile local
    
    3. Test with single module:
       nextflow run main.nf \\
           --sample_sheet platforms.tsv \\
           --skip_module2 --skip_module3 \\
           --skip_module4 --skip_module5 \\
           --skip_module6 --skip_module7 \\
           -profile local
    
    For more information:
        - README.md for detailed documentation
        - https://github.com/your-repo/pipeline
    
    ================================================================================
    """.stripIndent()
}

// ============================================================================
// PARAMETER VALIDATION
// ============================================================================

// Show help if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Check required parameters
if (!params.sample_sheet) {
    log.error "ERROR: --sample_sheet is required"
    helpMessage()
    exit 1
}

// Validate imputation parameters
if (params.run_topmed && !params.topmed_api_token) {
    log.error "ERROR: --topmed_api_token required when --run_topmed is set"
    exit 1
}

if (params.run_anvil && (!params.anvil_workspace || !params.anvil_billing_project)) {
    log.error "ERROR: --anvil_workspace and --anvil_billing_project required when --run_anvil is set"
    exit 1
}

// At least one imputation service must be enabled (unless skipping Module 2)
if (!params.skip_module2 && !params.run_topmed && !params.run_anvil) {
    log.error "ERROR: Must enable at least one imputation service (--run_topmed or --run_anvil)"
    exit 1
}

// ============================================================================
// PARSE SAMPLE SHEET
// ============================================================================

def parseSampleSheet(csv_file) {
    def samples = []
    
    csv_file.splitCsv(header: true, sep: '\t')
        .each { row ->
            samples << tuple(
                row.platform_id,
                file(row.file_path + ".bed"),
                file(row.file_path + ".bim"),
                file(row.file_path + ".fam")
            )
        }
    
    return samples
}

// ============================================================================
// MAIN WORKFLOW
// ============================================================================

workflow {
    
    // Print pipeline info
    log.info """
    ================================================================================
    Complete Genotyping Imputation & Ancestry Pipeline
    ================================================================================
    Sample sheet    : ${params.sample_sheet}
    Output directory: ${params.outdir}
    Input build     : ${params.input_build}
    
    Modules enabled:
      Module 1: Pre-imputation QC     - ENABLED
      Module 2: Imputation            - ${params.skip_module2 ? 'SKIPPED' : 'ENABLED'}
      Module 3: Post-imputation QC    - ${params.skip_module3 ? 'SKIPPED' : 'ENABLED'}
      Module 4: Platform merging      - ${params.skip_module4 ? 'SKIPPED' : 'ENABLED'}
      Module 5: Re-imputation         - ${params.skip_module5 ? 'SKIPPED' : 'ENABLED'}
      Module 6: Final QC              - ${params.skip_module6 ? 'SKIPPED' : 'ENABLED'}
      Module 7: Ancestry estimation   - ${params.skip_module7 ? 'SKIPPED' : 'ENABLED'}
    
    Imputation services:
      TOPMed: ${params.run_topmed ? 'ENABLED' : 'DISABLED'}
      AnVIL:  ${params.run_anvil ? 'ENABLED' : 'DISABLED'}
    
    ================================================================================
    """.stripIndent()
    
    // Parse sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.platform_id,
                file(row.file_path + ".bed"),
                file(row.file_path + ".bim"),
                file(row.file_path + ".fam")
            )
        }
        .set { sample_ch }
    
    // Track which modules to skip
    def skip_modules = []
    if (params.skip_module2) skip_modules << 2
    if (params.skip_module3) skip_modules << 3
    if (params.skip_module4) skip_modules << 4
    if (params.skip_module5) skip_modules << 5
    if (params.skip_module6) skip_modules << 6
    if (params.skip_module7) skip_modules << 7
    
    // ========================================================================
    // MODULE 1: Pre-Imputation QC (ALWAYS RUN)
    // ========================================================================
    
    MODULE1_PREIMPUTATION(sample_ch)
    module1_out = MODULE1_PREIMPUTATION.out
    
    // ========================================================================
    // MODULE 2: Imputation
    // ========================================================================
    
    if (!skip_modules.contains(2)) {
        MODULE2_IMPUTATION(module1_out.imputation_ready)
        module2_out = MODULE2_IMPUTATION.out
    } else {
        module2_out = [
            topmed_data: Channel.empty(),
            anvil_data: Channel.empty()
        ]
    }
    
    // ========================================================================
    // MODULE 3: Post-Imputation QC
    // ========================================================================
    
    if (!skip_modules.contains(3)) {
        // Prepare input channel with reference panel tracking
        if (params.run_topmed && params.run_anvil) {
            // Both services: create separate channels
            module2_out.topmed_data
                .map { platform, vcfs, infos -> 
                    tuple("topmed_ref", platform, vcfs, infos) 
                }
                .set { topmed_for_qc }
            
            module2_out.anvil_data
                .map { platform, vcfs, infos -> 
                    tuple("allofus_ref", platform, vcfs, infos) 
                }
                .set { anvil_for_qc }
            
            // Combine for parallel processing
            topmed_for_qc.mix(anvil_for_qc).set { qc3_input }
            
        } else if (params.run_topmed) {
            module2_out.topmed_data
                .map { platform, vcfs, infos -> 
                    tuple("topmed_ref", platform, vcfs, infos) 
                }
                .set { qc3_input }
        } else {
            module2_out.anvil_data
                .map { platform, vcfs, infos -> 
                    tuple("allofus_ref", platform, vcfs, infos) 
                }
                .set { qc3_input }
        }
        
        MODULE3_POSTQC(qc3_input)
        module3_out = MODULE3_POSTQC.out
    } else {
        log.warn "WARNING: Cannot skip Module 3 if Module 2 is enabled"
        exit 1
    }
    
    // ========================================================================
    // MODULE 4: Platform Merging
    // ========================================================================
    
    if (!skip_modules.contains(4)) {
        MODULE4_MERGING(module3_out.qc_data)
        module4_out = MODULE4_MERGING.out
    } else {
        log.warn "WARNING: Cannot skip Module 4 if Module 3 is enabled"
        exit 1
    }
    
    // ========================================================================
    // MODULE 5: Re-Imputation
    // ========================================================================
    
    if (!skip_modules.contains(5)) {
        MODULE5_REIMPUTATION(module4_out.merged_data)
        module5_out = MODULE5_REIMPUTATION.out
    } else {
        // Use merged data directly
        module5_out = [reimputed_data: module4_out.merged_data]
    }
    
    // ========================================================================
    // MODULE 6: Final QC (2nd Pass)
    // ========================================================================
    
    if (!skip_modules.contains(6)) {
        MODULE6_POSTMERGE_QC2(module5_out.reimputed_data)
        module6_out = MODULE6_POSTMERGE_QC2.out
    } else {
        module6_out = module5_out
    }
    
    // ========================================================================
    // MODULE 7: Ancestry Estimation
    // ========================================================================
    
    if (!skip_modules.contains(7)) {
        MODULE7_ANCESTRY(module6_out.ancestry_formatted)
        module7_out = MODULE7_ANCESTRY.out
    }
}

// ============================================================================
// WORKFLOW COMPLETION HANDLER
// ============================================================================

workflow.onComplete {
    def duration = workflow.duration
    def success = workflow.success
    
    log.info """
    ================================================================================
    Pipeline Execution Complete
    ================================================================================
    Status        : ${success ? 'SUCCESS ✓' : 'FAILED ✗'}
    Duration      : ${duration}
    Output        : ${params.outdir}
    
    ${params.run_topmed && params.run_anvil ? 
        "Both TOPMed and AnVIL results available in separate directories" : ""}
    
    Key Output Files:
    ${!params.skip_module6 ? 
        "- Final QC'd datasets: ${params.outdir}/module6/{ref_panel}/04_final_datasets/" : ""}
    ${!params.skip_module7 ? 
        "- Ancestry results: ${params.outdir}/module7/{ref_panel}/" : ""}
    
    Execution Reports:
    - Timeline: ${params.outdir}/reports/timeline.html
    - Report: ${params.outdir}/reports/execution_report.html
    - Trace: ${params.outdir}/reports/trace.txt
    
    ${success ? 'Pipeline completed successfully!' : 'Pipeline failed. Check logs for errors.'}
    ================================================================================
    """.stripIndent()
    
    if (!success) {
        log.error "Error: Pipeline execution failed"
        log.error "Check the error messages above and log files in ${params.outdir}/logs/"
    }
}

workflow.onError {
    log.error """
    ================================================================================
    Pipeline Error
    ================================================================================
    Error message: ${workflow.errorMessage}
    Error report : ${workflow.errorReport}
    
    Work directory: ${workflow.workDir}
    
    To resume from where the pipeline stopped:
    nextflow run main.nf -resume [same parameters]
    ================================================================================
    """.stripIndent()
}

// ============================================================================
// MANIFEST
// ============================================================================

manifest {
    name            = 'Genotyping Imputation & Ancestry Pipeline'
    author          = 'Your Name'
    homePage        = 'https://github.com/your-repo/pipeline'
    description     = 'Complete pipeline for genotyping QC, imputation, and ancestry estimation'
    mainScript      = 'main.nf'
    nextflowVersion = '>=22.04.0'
    version         = '1.0.0'
}
