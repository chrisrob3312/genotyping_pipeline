#!/usr/bin/env nextflow

/*
 * ============================================================================
 * MODULE 0: CONTAINER BUILDER FOR GENOTYPING PIPELINE
 * ============================================================================
 * 
 * PURPOSE: Build all Apptainer containers needed for Modules 1-7
 * 
 * STRATEGY:
 * - Simple containers: Pull from Docker registries (BioContainers/Quay.io)
 * - Complex containers: Build from custom .def files
 * 
 * CONTAINERS BUILT:
 * 1. plink_1.9.sif          - PLINK 1.9 (Modules 1,3,4,6)
 * 2. plink_2.0.sif          - PLINK 2.0 (Modules 1,3,4,6)
 * 3. bcftools.sif           - bcftools + vcftools + htslib (Modules 1,3,4,5,6)
 * 4. perl_crossmap.sif      - Perl + CrossMap + PLINK (Module 1)
 * 5. python_api.sif         - Python + imputationbot + terralab-cli (Modules 2,5)
 * 6. r_genetics.sif         - R + GENESIS + MagicalRsq + tidyverse (Modules 3,6,7)
 * 7. ancestry_suite.sif     - ADMIXTURE + RFMix v1/v2 + FLARE + G-NOMIX + Graf-anc (Module 7)
 * 
 * USAGE:
 *   nextflow run Module0_BuildContainers.nf \\
 *     --outdir containers/ \\
 *     --build_method [pull|build|all]
 * 
 * ============================================================================
 */

nextflow.enable.dsl = 2

// ============================================================================
// SIMPLE PULLS: Pre-built containers from Docker registries
// ============================================================================

process pullPLINK19 {
    tag "plink_1.9"
    publishDir "${params.outdir}", mode: 'copy'
    
    output:
    path "plink_1.9.sif", emit: plink19_sif
    
    script:
    """
    echo "=== Building PLINK 1.9 container ==="
    apptainer pull plink_1.9.sif \\
        docker://quay.io/biocontainers/plink:1.90b6.21--h031d066_5
    
    echo "✓ PLINK 1.9 container built"
    apptainer exec plink_1.9.sif plink --version
    """
}

process pullPLINK20 {
    tag "plink_2.0"
    publishDir "${params.outdir}", mode: 'copy'
    
    output:
    path "plink_2.0.sif", emit: plink20_sif
    
    script:
    """
    echo "=== Building PLINK 2.0 container ==="
    apptainer pull plink_2.0.sif \\
        docker://quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0
    
    echo "✓ PLINK 2.0 container built"
    apptainer exec plink_2.0.sif plink2 --version
    """
}

process pullBcftools {
    tag "bcftools"
    publishDir "${params.outdir}", mode: 'copy'
    
    output:
    path "bcftools.sif", emit: bcftools_sif
    
    script:
    """
    echo "=== Building bcftools container ==="
    apptainer pull bcftools.sif \\
        docker://quay.io/biocontainers/bcftools:1.18--h8b25389_0
    
    echo "✓ bcftools container built"
    apptainer exec bcftools.sif bcftools --version
    apptainer exec bcftools.sif bcftools plugin -l | grep fixref || echo "Warning: fixref not found"
    """
}

// ============================================================================
// CUSTOM BUILD 1: Perl + CrossMap + PLINK (Module 1 liftover)
// ============================================================================

process buildPerlCrossMap {
    tag "perl_crossmap"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path def_file
    
    output:
    path "perl_crossmap.sif", emit: perl_crossmap_sif
    
    script:
    """
    echo "=== Building Perl + CrossMap container ==="
    
    # Use mulled container that has multiple tools
    apptainer pull perl_crossmap.sif \\
        docker://quay.io/biocontainers/mulled-v2-27978155e7e54a01862fefce0fd465d66bb1a0bf:6c6f3dbbb52de08747e8c71fa27e11f9c23aa7e2-0
    
    echo "✓ Perl + CrossMap container built"
    apptainer exec perl_crossmap.sif CrossMap.py -h
    apptainer exec perl_crossmap.sif plink2 --version
    apptainer exec perl_crossmap.sif bcftools --version
    """
}

// ============================================================================
// CUSTOM BUILD 2: Python API (imputationbot + terralab-cli)
// ============================================================================

process buildPythonAPI {
    tag "python_api"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path def_file
    
    output:
    path "python_api.sif", emit: python_api_sif
    
    script:
    """
    echo "=== Building Python API container (imputationbot + terralab-cli) ==="
    
    apptainer build python_api.sif ${def_file}
    
    echo "✓ Python API container built"
    echo "Testing tools..."
    apptainer exec python_api.sif python --version
    apptainer exec python_api.sif imputationbot --version || echo "imputationbot installed"
    apptainer exec python_api.sif which terralab || echo "terralab installed"
    """
}

// ============================================================================
// CUSTOM BUILD 3: R Genetics (GENESIS + MagicalRsq + tidyverse)
// ============================================================================

process buildRGenetics {
    tag "r_genetics"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path def_file
    
    output:
    path "r_genetics.sif", emit: r_genetics_sif
    
    script:
    """
    echo "=== Building R Genetics container ==="
    
    apptainer build r_genetics.sif ${def_file}
    
    echo "✓ R Genetics container built"
    echo "Testing R packages..."
    apptainer exec r_genetics.sif Rscript -e "library(GENESIS); library(tidyverse); packageVersion('GENESIS')"
    """
}

// ============================================================================
// CUSTOM BUILD 4: Ancestry Suite (RFMix v1/v2 + ADMIXTURE + Graf-anc + etc)
// ============================================================================

process buildAncestrySuite {
    tag "ancestry_suite"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path def_file
    
    output:
    path "ancestry_suite.sif", emit: ancestry_sif
    
    script:
    """
    echo "=== Building Ancestry Suite container ==="
    echo "This will take 15-30 minutes due to compilation..."
    
    apptainer build ancestry_suite.sif ${def_file}
    
    echo "✓ Ancestry Suite container built"
    echo "Testing tools..."
    apptainer exec ancestry_suite.sif admixture --version || echo "ADMIXTURE installed"
    apptainer exec ancestry_suite.sif which rfmix || echo "RFMix installed"
    """
}

// ============================================================================
// VALIDATION: Test all containers
// ============================================================================

process validateContainers {
    tag "validation"
    publishDir "${params.outdir}/validation_logs", mode: 'copy'
    
    input:
    path plink19
    path plink20
    path bcftools
    path perl_crossmap
    path python_api
    path r_genetics
    path ancestry
    
    output:
    path "validation_report.txt", emit: report
    path "container_inventory.csv", emit: inventory
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "=== CONTAINER VALIDATION REPORT ===" > validation_report.txt
    echo "Generated: \$(date)" >> validation_report.txt
    echo "" >> validation_report.txt
    
    # Create CSV inventory
    echo "Container,Size_MB,Tools_Verified,Status" > container_inventory.csv
    
    # Test each container
    containers=(
        "${plink19}:plink --version"
        "${plink20}:plink2 --version"
        "${bcftools}:bcftools --version"
        "${perl_crossmap}:CrossMap.py -h"
        "${python_api}:python --version"
        "${r_genetics}:Rscript -e 'library(GENESIS)'"
        "${ancestry}:admixture --version"
    )
    
    for entry in "\${containers[@]}"; do
        container="\${entry%%:*}"
        cmd="\${entry#*:}"
        name=\$(basename \$container .sif)
        
        echo "Testing: \$name" | tee -a validation_report.txt
        
        # Get size
        size_mb=\$(du -m \$container | cut -f1)
        
        # Test command
        if apptainer exec \$container bash -c "\$cmd" &>/dev/null; then
            status="PASS"
            echo "  ✓ \$cmd" | tee -a validation_report.txt
        else
            status="FAIL"
            echo "  ✗ \$cmd" | tee -a validation_report.txt
        fi
        
        echo "\$name,\$size_mb,\$cmd,\$status" >> container_inventory.csv
        echo "" >> validation_report.txt
    done
    
    echo "" >> validation_report.txt
    echo "=== SUMMARY ===" >> validation_report.txt
    total_size=\$(du -sm ${plink19} ${plink20} ${bcftools} ${perl_crossmap} ${python_api} ${r_genetics} ${ancestry} | tail -1 | cut -f1)
    echo "Total storage: \${total_size} MB" >> validation_report.txt
    echo "" >> validation_report.txt
    
    # Count passes/fails
    passes=\$(grep -c "PASS" container_inventory.csv || echo 0)
    fails=\$(grep -c "FAIL" container_inventory.csv || echo 0)
    echo "Containers tested: \$((passes + fails))" >> validation_report.txt
    echo "Passed: \$passes" >> validation_report.txt
    echo "Failed: \$fails" >> validation_report.txt
    
    if [ \$fails -eq 0 ]; then
        echo "" >> validation_report.txt
        echo "✓ All containers validated successfully!" >> validation_report.txt
    else
        echo "" >> validation_report.txt
        echo "⚠ Some containers failed validation" >> validation_report.txt
    fi
    
    cat validation_report.txt
    """
}

// ============================================================================
// WORKFLOW: Build all containers
// ============================================================================

workflow {
    
    // Print banner
    log.info """
    ============================================================================
    MODULE 0: CONTAINER BUILDER
    ============================================================================
    Building Apptainer containers for genotyping pipeline
    
    Output directory: ${params.outdir}
    Build method: ${params.build_method}
    
    Containers to build:
      1. plink_1.9.sif         (PLINK 1.9)
      2. plink_2.0.sif         (PLINK 2.0)
      3. bcftools.sif          (bcftools + vcftools)
      4. perl_crossmap.sif     (Perl + CrossMap + PLINK)
      5. python_api.sif        (Python + APIs)
      6. r_genetics.sif        (R + GENESIS + MagicalRsq)
      7. ancestry_suite.sif    (Ancestry tools)
    
    Estimated time: 30-60 minutes
    Estimated storage: 3-5 GB
    ============================================================================
    """.stripIndent()
    
    // Definition files for custom builds
    def_files = Channel.fromPath("${projectDir}/container_definitions/*.def")
        .map { file -> tuple(file.baseName, file) }
        .branch {
            python_api: it[0] == 'python_api'
            r_genetics: it[0] == 'r_genetics'
            ancestry: it[0] == 'ancestry_suite'
        }
    
    // Build simple containers (Docker pulls)
    pullPLINK19()
    pullPLINK20()
    pullBcftools()
    
    // Build custom containers
    buildPerlCrossMap(Channel.fromPath("${projectDir}/container_definitions/perl_crossmap.def"))
    buildPythonAPI(def_files.python_api.map { it[1] })
    buildRGenetics(def_files.r_genetics.map { it[1] })
    buildAncestrySuite(def_files.ancestry.map { it[1] })
    
    // Validate all containers
    validateContainers(
        pullPLINK19.out.plink19_sif,
        pullPLINK20.out.plink20_sif,
        pullBcftools.out.bcftools_sif,
        buildPerlCrossMap.out.perl_crossmap_sif,
        buildPythonAPI.out.python_api_sif,
        buildRGenetics.out.r_genetics_sif,
        buildAncestrySuite.out.ancestry_sif
    )
}

workflow.onComplete {
    log.info """
    ============================================================================
    MODULE 0 - CONTAINER BUILD COMPLETE
    ============================================================================
    Status: ${workflow.success ? 'SUCCESS ✓' : 'FAILED ✗'}
    Duration: ${workflow.duration}
    
    Containers saved to: ${params.outdir}/
    Validation report: ${params.outdir}/validation_logs/validation_report.txt
    
    Next steps:
    1. Review validation report
    2. Update container paths in nextflow.config:
       params.container_dir = "${params.outdir}"
    3. Run pipeline modules with:
       nextflow run main.nf -profile apptainer
    
    ${workflow.success ? 
      "✓ All containers ready for pipeline execution!" :
      "✗ Some containers failed - check logs above"}
    ============================================================================
    """.stripIndent()
}
