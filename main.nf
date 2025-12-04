#!/usr/bin/env nextflow

/*
================================================================================
    Complete Genotyping Imputation and Ancestry Pipeline
================================================================================
    Seven-module pipeline for:
    Module 1: Pre-imputation QC and file preparation  
    Module 2: Imputation (TOPMed & All of Us AnVIL) - parallelized by platform
    Module 3: Post-imputation QC with MagicalRsq-X
    Module 4: Platform merging with QC
    Module 5: Re-imputation of merged data
    Module 6: Post-merge QC (2nd pass) with GENESIS/PLINK
    Module 7: Ancestry estimation (global and local)
================================================================================
*/

nextflow.enable.dsl=2

// ============================================================================
// INCLUDE MODULES
// ============================================================================

include { PRE_IMPUTATION_QC } from './modules/Module1_PreImputation'
include { MODULE2_IMPUTATION as IMPUTATION } from './modules/Module2_Imputation'
include { MODULE3_POSTQC as POST_IMPUTATION_QC } from './modules/Module3_PostQC'
include { MODULE4_MERGING as PLATFORM_MERGING } from './modules/Module4_Merging'
include { MODULE5_REIMPUTATION as RE_IMPUTATION } from './modules/Module5_ReImputation'
include { MODULE6_POSTMERGE_QC as POST_MERGE_QC } from './modules/Module6_PostMergeQC'
include { MODULE7_ANCESTRY as ANCESTRY_ESTIMATION } from './modules/Module7_Ancestry'

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
        --sample_sheet PATH       TSV file with columns: platform,bed,bim,fam
    
    Optional Arguments - General:
        --input_build STR         Genome build [default: hg19]
        --outdir PATH             Output directory [default: results]
        --skip_modules LIST       Comma-separated modules to skip (e.g., "5,7")
    
    Optional Arguments - Module 1:
        --reference_fasta PATH    Human reference genome (hg38)
        --liftover_chain PATH     hg19->hg38 chain file
        --topmed_reference PATH   TOPMed reference for strand checking
    
    Optional Arguments - Modules 2 & 5 (Imputation):
        --run_topmed              Run TOPMed imputation [default: true]
        --run_anvil               Run All of Us AnVIL [default: false]
        --topmed_api_token STR    TOPMed API token (required if run_topmed)
        --topmed_password STR     TOPMed download password
        --anvil_workspace STR     AnVIL Terra workspace
        --anvil_project STR       AnVIL Google project
    
    Optional Arguments - Modules 3 & 6 (QC):
        --magicalrsq_threshold NUM    MagicalRsq-X R² threshold [default: 0.3]
        --sample_call_rate NUM        Sample call rate [default: 0.95]
        --variant_call_rate NUM       Variant call rate [default: 0.95]
        --hwe_pvalue NUM              HWE p-value [default: 1e-6]
        --use_genesis                 Use GENESIS PCRelate [default: true]
    
    Optional Arguments - Module 7:
        --reference_panel PATH    Ancestry reference panel
        --genetic_map_dir PATH    Genetic map directory for RFMix
        --admixture_k LIST        K values for ADMIXTURE [default: "5,6,7"]
        --lai_methods LIST        LAI methods (rfmix2,rfmix1,flare,gnomix)
        --run_lai                 Enable local ancestry inference
    
    Execution Profiles:
        -profile standard         Local execution
        -profile apptainer        Use Apptainer containers (recommended)
        -profile singularity      Use Singularity containers
        -profile docker           Use Docker containers (testing only)
        -profile slurm            SLURM cluster execution
        -profile pbs              PBS/Torque cluster execution
        -profile google           Google Cloud (All of Us AnVIL)
        -profile aws              AWS Batch execution
    
    Examples:
        # Full pipeline with TOPMed on SLURM cluster
        nextflow run main.nf \\
            --sample_sheet platforms.tsv \\
            --run_topmed \\
            --topmed_api_token \$TOKEN \\
            -profile slurm,apptainer \\
            -resume
        
        # Skip re-imputation and ancestry for faster testing
        nextflow run main.nf \\
            --sample_sheet platforms.tsv \\
            --skip_modules "5,7" \\
            -profile apptainer
        
        # All of Us AnVIL imputation on Google Cloud
        nextflow run main.nf \\
            --sample_sheet platforms.tsv \\
            --run_anvil \\
            --anvil_workspace my-workspace \\
            --anvil_project my-project \\
            -profile google,apptainer
    
    Output Directory Structure:
        ${params.outdir}/
        ├── module1/          Pre-imputation QC per platform
        ├── module2/          Imputation results per platform
        ├── module3/          Post-imputation QC per platform
        ├── module4/          Merged platforms
        ├── module5/          Re-imputed merged data
        ├── module6/          Final QC'd datasets ★ ANALYSIS-READY ★
        │   ├── final_maf01/     MAF >= 0.01 filtered
        │   └── final_nomaf/     No MAF filter
        └── module7/          Ancestry results
            ├── grafanc/         Global ancestry
            ├── admixture/       Population structure
            └── lai/             Local ancestry (if enabled)
    
    Documentation:
        See README.md and docs/ directory for detailed information.
    ================================================================================
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// ============================================================================
// PARAMETER VALIDATION
// ============================================================================

// Check required parameters
if (!params.sample_sheet) {
    log.error "ERROR: --sample_sheet is required!"
    helpMessage()
    exit 1
}

// Check sample sheet exists
if (!file(params.sample_sheet).exists()) {
    log.error "ERROR: Sample sheet not found: ${params.sample_sheet}"
    exit 1
}

// Parse modules to skip
def skip_modules = params.skip_modules ? 
    params.skip_modules.split(',').collect { it.toInteger() } : []

// Validate imputation parameters
if (params.run_topmed && !params.topmed_api_token) {
    log.warn "WARNING: TOPMed imputation enabled but no API token provided!"
    log.warn "Set --topmed_api_token or disable with --run_topmed false"
}

if (params.run_anvil && (!params.anvil_workspace || !params.anvil_project)) {
    log.warn "WARNING: AnVIL imputation enabled but workspace/project not set!"
    log.warn "Set --anvil_workspace and --anvil_project or disable with --run_anvil false"
}

// ============================================================================
// PRINT CONFIGURATION
// ============================================================================

log.info """
================================================================================
Pipeline Configuration
================================================================================
Sample sheet      : ${params.sample_sheet}
Output directory  : ${params.outdir}
Input build       : ${params.input_build}
Imputation        : TOPMed=${params.run_topmed}, AnVIL=${params.run_anvil}
Use GENESIS       : ${params.use_genesis}
Modules to skip   : ${skip_modules.size() > 0 ? skip_modules : 'None'}

Modules enabled:
  Module 1 (Pre-QC)        : ${!skip_modules.contains(1)}
  Module 2 (Imputation)    : ${!skip_modules.contains(2)}
  Module 3 (Post-QC)       : ${!skip_modules.contains(3)}
  Module 4 (Merging)       : ${!skip_modules.contains(4)}
  Module 5 (Re-imputation) : ${!skip_modules.contains(5)}
  Module 6 (Final QC)      : ${!skip_modules.contains(6)}
  Module 7 (Ancestry)      : ${!skip_modules.contains(7)}
================================================================================
""".stripIndent()

// ============================================================================
// PARSE INPUT SAMPLE SHEET
// ============================================================================

/*
 * Expected format (TSV):
 * platform	bed	bim	fam
 * platform1	path/to/data1.bed	path/to/data1.bim	path/to/data1.fam
 * platform2	path/to/data2.bed	path/to/data2.bim	path/to/data2.fam
 */

Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header: true, sep: '\t')
    .map { row -> 
        def platform = row.platform
        def bed = file(row.bed)
        def bim = file(row.bim)
        def fam = file(row.fam)
        
        // Validate files exist
        if (!bed.exists()) error "BED file not found: ${bed}"
        if (!bim.exists()) error "BIM file not found: ${bim}"
        if (!fam.exists()) error "FAM file not found: ${fam}"
        
        return tuple(platform, bed, bim, fam)
    }
    .set { input_plink_ch }

// ============================================================================
// MAIN WORKFLOW
// ============================================================================

workflow {
    log.info "Starting pipeline execution..."
    
    // ------------------------------------------------------------------------
    // MODULE 1: Pre-imputation QC and file preparation
    // ------------------------------------------------------------------------
    if (!skip_modules.contains(1)) {
        log.info "=== Running Module 1: Pre-Imputation QC ==="
        
        PRE_IMPUTATION_QC(
            input_plink_ch,
            file(params.reference_fasta),
            file(params.liftover_chain),
            file(params.topmed_reference),
            file(params.strand_check_script)
        )
        
        // Outputs: per-platform VCFs (bgzipped, indexed, by chromosome)
        module1_vcfs = PRE_IMPUTATION_QC.out.vcf_files
        module1_manifests = PRE_IMPUTATION_QC.out.manifests
    } else {
        error "Cannot skip Module 1 - it's required for all downstream modules"
    }
    
    // ------------------------------------------------------------------------
    // MODULE 2: Imputation via TOPMed and/or All of Us AnVIL
    // ------------------------------------------------------------------------
    if (!skip_modules.contains(2)) {
        log.info "=== Running Module 2: Imputation ==="
        
        IMPUTATION(
            module1_vcfs,
            module1_manifests,
            params.run_topmed,
            params.run_anvil,
            params.topmed_api_token,
            params.topmed_password,
            params.anvil_workspace,
            params.anvil_project
        )
        
        // Outputs: imputed VCFs per platform (bgzipped, indexed, by chromosome)
        module2_imputed = IMPUTATION.out.imputed_vcfs
    } else {
        error "Cannot skip Module 2 - imputation is required"
    }
    
    // ------------------------------------------------------------------------
    // MODULE 3: Post-imputation QC with MagicalRsq-X
    // ------------------------------------------------------------------------
    if (!skip_modules.contains(3)) {
        log.info "=== Running Module 3: Post-Imputation QC ==="
        
        POST_IMPUTATION_QC(
            module2_imputed,
            file(params.magicalrsq_filter_script),
            file(params.magicalrsq_models_dir),
            params.magicalrsq_threshold,
            params.sample_call_rate,
            params.variant_call_rate
        )
        
        // Outputs: QC'd VCFs per platform (bgzipped, indexed, by chromosome)
        module3_qc_data = POST_IMPUTATION_QC.out.qc_vcfs
    } else {
        log.info "Skipping Module 3 - using Module 2 output directly"
        module3_qc_data = module2_imputed
    }
    
    // ------------------------------------------------------------------------
    // MODULE 4: Platform merging
    // ------------------------------------------------------------------------
    if (!skip_modules.contains(4)) {
        log.info "=== Running Module 4: Platform Merging ==="
        
        // Collect all platforms for merging
        PLATFORM_MERGING(
            module3_qc_data.collect()
        )
        
        // Outputs: single merged PLINK file (bed/bim/fam)
        module4_merged = PLATFORM_MERGING.out.merged_plink
    } else {
        error "Cannot skip Module 4 - merging is required for downstream analysis"
    }
    
    // ------------------------------------------------------------------------
    // MODULE 5: Re-imputation of merged data
    // ------------------------------------------------------------------------
    if (!skip_modules.contains(5)) {
        log.info "=== Running Module 5: Re-Imputation ==="
        
        RE_IMPUTATION(
            module4_merged,
            params.run_topmed,
            params.run_anvil,
            params.topmed_api_token,
            params.topmed_password,
            params.anvil_workspace,
            params.anvil_project
        )
        
        // Outputs: re-imputed VCFs (bgzipped, indexed, by chromosome)
        module5_reimputed = RE_IMPUTATION.out.reimputed_vcfs
    } else {
        log.info "Skipping Module 5 - using Module 4 merged data for Module 6"
        // Convert Module 4 PLINK to VCF format for Module 6
        module5_reimputed = module4_merged
    }
    
    // ------------------------------------------------------------------------
    // MODULE 6: Post-merge QC (2nd pass)
    // ------------------------------------------------------------------------
    if (!skip_modules.contains(6)) {
        log.info "=== Running Module 6: Post-Merge QC (Final) ==="
        
        POST_MERGE_QC(
            module5_reimputed,
            file(params.magicalrsq_filter_script),
            file(params.magicalrsq_models_dir),
            params.use_genesis,
            params.magicalrsq_threshold,
            params.hwe_pvalue,
            params.het_sd_threshold,
            params.kinship_threshold
        )
        
        // Outputs: final analysis-ready datasets (MAF filtered and unfiltered)
        module6_final_maf = POST_MERGE_QC.out.final_data_maf
        module6_final_nomaf = POST_MERGE_QC.out.final_data_nomaf
        module6_pcs = POST_MERGE_QC.out.pca_results
    } else {
        log.info "Skipping Module 6 - using Module 5 output directly"
        module6_final_maf = module5_reimputed
        module6_final_nomaf = module5_reimputed
    }
    
    // ------------------------------------------------------------------------
    // MODULE 7: Ancestry estimation
    // ------------------------------------------------------------------------
    if (!skip_modules.contains(7)) {
        log.info "=== Running Module 7: Ancestry Estimation ==="
        
        // Parse ADMIXTURE K values
        def k_values = params.admixture_k.split(',').collect { it.toInteger() }
        
        // Parse LAI methods
        def lai_methods = params.lai_methods.split(',')
        
        ANCESTRY_ESTIMATION(
            module6_final_nomaf,
            file(params.reference_panel),
            file(params.genetic_map_dir),
            k_values,
            params.run_lai,
            lai_methods
        )
        
        // Outputs: ancestry estimation results
        module7_grafanc = ANCESTRY_ESTIMATION.out.grafanc_results
        module7_admixture = ANCESTRY_ESTIMATION.out.admixture_results
        module7_lai = ANCESTRY_ESTIMATION.out.lai_results
    } else {
        log.info "Skipping Module 7 - Ancestry estimation disabled"
    }
    
    log.info "=== Pipeline workflow defined successfully ==="
}

// ============================================================================
// WORKFLOW COMPLETION HANDLERS
// ============================================================================

workflow.onComplete {
    def duration = workflow.duration
    def success_status = workflow.success ? "SUCCESS ✓" : "FAILED ✗"
    
    log.info """
    ================================================================================
    Pipeline Execution Summary
    ================================================================================
    Status        : ${success_status}
    Duration      : ${duration}
    Completed at  : ${workflow.complete}
    Work directory: ${workflow.workDir}
    Exit status   : ${workflow.exitStatus}
    Error report  : ${workflow.errorReport ?: 'None'}
    ================================================================================
    """.stripIndent()
    
    if (workflow.success) {
        log.info """
        Pipeline completed successfully!
        
        ✓ Analysis-ready datasets available in:
          - ${params.outdir}/module6/final_datasets/
        
        ✓ Ancestry results available in:
          - ${params.outdir}/module7/
        
        ✓ QC reports and logs in each module directory
        
        Next steps:
        1. Review QC reports in ${params.outdir}/module*/reports/
        2. Check ancestry results in ${params.outdir}/module7/
        3. Use final datasets in ${params.outdir}/module6/final_datasets/
        4. See execution report: ${params.tracedir}/execution_report.html
        ================================================================================
        """.stripIndent()
    } else {
        log.error """
        Pipeline failed!
        
        Troubleshooting:
        1. Check error logs in: ${workflow.workDir}
        2. Review trace file: ${params.tracedir}/execution_trace.txt
        3. Check module-specific logs in ${params.outdir}/
        
        Common issues:
        - API token errors: Verify TOPMed/AnVIL credentials
        - Memory errors: Increase resources in nextflow.config
        - Missing files: Verify all paths in sample sheet
        - Container issues: Ensure Apptainer/Singularity is properly installed
        
        For help, see documentation in README.md
        ================================================================================
        """.stripIndent()
    }
}

workflow.onError {
    log.error "Pipeline execution stopped with error:"
    log.error "  ${workflow.errorMessage}"
    log.error "  ${workflow.errorReport}"
}
