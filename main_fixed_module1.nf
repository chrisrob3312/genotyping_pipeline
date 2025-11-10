#!/usr/bin/env nextflow

/*
================================================================================
    Complete Genotyping Imputation and Ancestry Pipeline
================================================================================
    Seven-module pipeline for:
    Module 1: Pre-imputation QC and preparation  
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
include { IMPUTATION } from './modules/Module2_Imputation'
include { POST_IMPUTATION_QC } from './modules/Module3_PostQC'
include { PLATFORM_MERGING } from './modules/Module4_Merging'
include { RE_IMPUTATION } from './modules/Module5_ReImputation'
include { POST_MERGE_QC } from './modules/Module6_PostMergeQC'
include { ANCESTRY_ESTIMATION } from './modules/Module7_Ancestry'

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
        nextflow run main.nf --sample_sheet samples.csv [options]
    
    Required Arguments:
        --sample_sheet PATH       CSV file with platform/batch information
    
    Sample Sheet Format (CSV):
        platform_id,batch_id,input_path,file_type,build,file_structure
        GSAv1,batch1,/data/gsa1/batch1,plink,hg19,individual_samples
        GSAv2,batch2,/data/gsa2/batch2,plink,hg38,individual_chr_split
        
    File Structure Options:
        - individual_samples     : One file per sample, all chromosomes together
        - individual_chr_split   : One file per sample per chromosome
        - merged_batch          : Pre-merged batch file (all samples, all chr)
        - merged_chr_split      : Pre-merged batch, split by chromosome
    
    Optional Arguments - General:
        --outdir PATH             Output directory [default: results]
        --skip_modules LIST       Comma-separated modules to skip (e.g., "5,7")
    
    Optional Arguments - Module 1 (Pre-imputation QC):
        --hg19_fasta PATH         hg19 reference genome FASTA
        --hg38_fasta PATH         hg38 reference genome FASTA
        --liftover_chain PATH     hg19->hg38 chain file
        --topmed_reference PATH   TOPMed Freeze 10 reference for strand checking
    
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
    
    Optional Arguments - Module 7 (Ancestry):
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
            --sample_sheet platforms.csv \\
            --run_topmed \\
            --topmed_api_token \$TOKEN \\
            -profile slurm,apptainer \\
            -resume
        
        # Skip re-imputation and ancestry for faster testing
        nextflow run main.nf \\
            --sample_sheet platforms.csv \\
            --skip_modules "5,7" \\
            -profile apptainer
        
        # All of Us AnVIL imputation on Google Cloud
        nextflow run main.nf \\
            --sample_sheet platforms.csv \\
            --run_anvil \\
            --anvil_workspace my-workspace \\
            --anvil_project my-project \\
            -profile google,apptainer
    
    Output Directory Structure:
        \${params.outdir}/
        ├── <platform_id>/
        │   ├── <batch_id>/
        │   │   ├── 01_sample_prep/       Sample format conversion
        │   │   ├── 02_batch_merge/       Batch-level merging
        │   │   ├── 03_hg19_fixes/        Pre-liftover alignment (if hg19)
        │   │   ├── 04_liftover/          hg19→hg38 conversion (if needed)
        │   │   └── ...
        │   ├── 05_platform_merge/        Platform-level union merge
        │   ├── 06_topmed_validation/     TOPMed strand check (QC)
        │   ├── 06_anvil_validation/      AnVIL validation (QC)
        │   ├── 07_topmed_qc/            Light QC for TOPMed
        │   ├── 07_anvil_qc/             Light QC for AnVIL
        │   └── 08_topmed_vcfs/          Service-specific VCFs by chr
        │       └── 08_anvil_vcfs/
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
        Input setup: See input_file_setup_v2.md for sample sheet examples.
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

// Validate Module 1 reference files
if (!file(params.hg19_fasta).exists()) {
    log.error "ERROR: hg19 reference not found: ${params.hg19_fasta}"
    exit 1
}

if (!file(params.hg38_fasta).exists()) {
    log.error "ERROR: hg38 reference not found: ${params.hg38_fasta}"
    exit 1
}

if (!file(params.liftover_chain).exists()) {
    log.error "ERROR: Liftover chain file not found: ${params.liftover_chain}"
    exit 1
}

if (!file(params.topmed_reference).exists()) {
    log.error "ERROR: TOPMed reference not found: ${params.topmed_reference}"
    exit 1
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
Imputation        : TOPMed=${params.run_topmed}, AnVIL=${params.run_anvil}
Use GENESIS       : ${params.use_genesis}
Modules to skip   : ${skip_modules.size() > 0 ? skip_modules : 'None'}

Reference Files:
  hg19 FASTA      : ${params.hg19_fasta}
  hg38 FASTA      : ${params.hg38_fasta}
  Liftover chain  : ${params.liftover_chain}
  TOPMed reference: ${params.topmed_reference}

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
 * Expected format (CSV):
 * platform_id,batch_id,input_path,file_type,build,file_structure
 * GSAv1,batch1,/data/gsa1/batch1,plink,hg19,individual_samples
 * GSAv2,batch2,/data/gsa2/batch2,plink,hg38,individual_chr_split
 * MEGA,batch3,/data/mega/batch3_merged,plink,hg38,merged_batch
 * 
 * File structure options:
 * - individual_samples: One file per sample, all chr together
 * - individual_chr_split: One file per sample per chr
 * - merged_batch: Pre-merged batch (all samples, all chr)
 * - merged_chr_split: Pre-merged batch, split by chr
 */

Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header: true, sep: ',')
    .map { row -> 
        def platform_id = row.platform_id
        def batch_id = row.batch_id
        def input_path = row.input_path
        def file_type = row.file_type
        def build = row.build
        def file_structure = row.file_structure
        
        // Validate file_type
        if (!(file_type in ['plink', 'vcf'])) {
            error "Invalid file_type '${file_type}' for ${platform_id}/${batch_id}. Must be 'plink' or 'vcf'"
        }
        
        // Validate build
        if (!(build in ['hg19', 'hg38'])) {
            error "Invalid build '${build}' for ${platform_id}/${batch_id}. Must be 'hg19' or 'hg38'"
        }
        
        // Validate file_structure
        def valid_structures = ['individual_samples', 'individual_chr_split', 'merged_batch', 'merged_chr_split']
        if (!(file_structure in valid_structures)) {
            error "Invalid file_structure '${file_structure}' for ${platform_id}/${batch_id}. Must be one of: ${valid_structures.join(', ')}"
        }
        
        // Validate input_path exists
        if (file_structure in ['individual_samples', 'individual_chr_split']) {
            // Should be a directory
            if (!file(input_path).exists() || !file(input_path).isDirectory()) {
                error "Input path '${input_path}' for ${platform_id}/${batch_id} should be a directory but was not found or is not a directory"
            }
        } else {
            // Should be a file prefix (files should exist)
            def test_file = file_type == 'plink' ? "${input_path}.bed" : "${input_path}.vcf.gz"
            if (!file(test_file).exists()) {
                error "Input file '${test_file}' for ${platform_id}/${batch_id} not found"
            }
        }
        
        log.info "Loaded: ${platform_id}/${batch_id} | ${file_type} | ${build} | ${file_structure} | ${input_path}"
        
        return tuple(platform_id, batch_id, input_path, file_type, build, file_structure)
    }
    .set { sample_sheet_ch }

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
        log.info "Module 1 will:"
        log.info "  1. Discover and load samples based on file_structure"
        log.info "  2. Convert formats (VCF→PLINK, concatenate chr splits)"
        log.info "  3. Merge samples → batches (UNION merge, NO QC)"
        log.info "  4. Align to reference (hg19 or hg38)"
        log.info "  5. Liftover hg19 → hg38 (if needed)"
        log.info "  6. UNION merge batches → platform"
        log.info "  7. Service-specific validation (TOPMed + AnVIL)"
        log.info "  8. Light QC (first time QC happens)"
        log.info "  9. Create service-specific VCFs by chromosome"
        log.info ""
        
        PRE_IMPUTATION_QC(
            sample_sheet_ch
        )
        
        // Outputs from Module 1:
        // - imputation_vcfs: tuple(platform_id, service, chr, vcf, vcf_index)
        // - qc_logs, union_logs, vcf_logs, topmed_validation_logs, anvil_validation_logs
        module1_imputation_vcfs = PRE_IMPUTATION_QC.out.imputation_vcfs
        
    } else {
        error "Cannot skip Module 1 - it's required for all downstream modules"
    }
    
    // ------------------------------------------------------------------------
    // MODULE 2: Imputation via TOPMed and/or All of Us AnVIL
    // ------------------------------------------------------------------------
    if (!skip_modules.contains(2)) {
        log.info "=== Running Module 2: Imputation ==="
        
        // Group VCFs by platform and service for submission
        module1_imputation_vcfs
            .groupTuple(by: [0, 1])  // Group by platform_id and service
            .set { grouped_for_imputation }
        
        IMPUTATION(
            grouped_for_imputation,
            params.run_topmed,
            params.run_anvil,
            params.topmed_api_token,
            params.topmed_password,
            params.anvil_workspace,
            params.anvil_project
        )
        
        // Outputs: imputed VCFs per platform per service
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
        
        ✓ Module 1 outputs (per platform):
          - Platform-merged data: ${params.outdir}/<platform_id>/05_platform_merge/
          - Service-specific VCFs: ${params.outdir}/<platform_id>/08_*_vcfs/
        
        ✓ Analysis-ready datasets available in:
          - ${params.outdir}/module6/final_datasets/
        
        ✓ Ancestry results available in:
          - ${params.outdir}/module7/
        
        ✓ QC reports and logs in each module directory
        
        Next steps:
        1. Review QC reports in ${params.outdir}/module*/reports/
        2. Check Module 1 logs:
           - Union merge: ${params.outdir}/<platform_id>/05_platform_merge/*_union_merge.log
           - Validation: ${params.outdir}/<platform_id>/06_*_validation/*_validation.log
           - Light QC: ${params.outdir}/<platform_id>/07_*_qc/*_light_qc.log
        3. Check ancestry results in ${params.outdir}/module7/
        4. Use final datasets in ${params.outdir}/module6/final_datasets/
        5. See execution report: ${params.tracedir}/execution_report.html
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
        - Sample sheet format: Ensure CSV with correct columns
        - File paths: Verify all paths in sample sheet exist
        - Genome build: Check hg19/hg38 specified correctly
        - File structure: Verify file_structure matches actual data layout
        - Reference files: Ensure all reference files exist and are indexed
        - API token errors: Verify TOPMed/AnVIL credentials
        - Memory errors: Increase resources in nextflow.config
        - Container issues: Ensure Apptainer/Singularity is properly installed
        
        For help, see documentation in README.md and input_file_setup_v2.md
        ================================================================================
        """.stripIndent()
    }
}

workflow.onError {
    log.error "Pipeline execution stopped with error:"
    log.error "  ${workflow.errorMessage}"
    log.error "  ${workflow.errorReport}"
}
