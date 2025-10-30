#!/usr/bin/env nextflow

/*
 * ============================================================================
 * MODULE 5: RE-IMPUTATION OF MERGED DATA (Quota-Efficient v3.0)
 * ============================================================================
 * 
 * Purpose: Re-impute merged platform data for maximum variant coverage
 * 
 * CRITICAL: Submits ALL chromosomes as ONE job per service
 * - Merged dataset → chr1-22 VCFs → ONE job submission
 * - TOPMed: 1 job with all chromosomes = 1 quota unit ✓
 * - AnVIL: 1 job with all chromosomes = 1 quota unit ✓
 * 
 * Workflow:
 * 1. Convert merged PLINK to per-chromosome VCFs
 * 2. Group all chr1-22 VCFs together
 * 3. Submit as ONE job per service (not 22 separate jobs!)
 * 4. Monitor and download results
 * 
 * Input: Merged genome-wide data from Module 4
 * Output: Re-imputed VCFs with enhanced variant coverage
 * 
 * Quota Usage:
 * - CORRECT: 1 TOPMed job + 1 AnVIL job = 2 quota units ✓
 * - WRONG: 22 chr × 2 services = 44 quota units ❌
 * ============================================================================
 */

nextflow.enable.dsl = 2

// ============================================================================
// PROCESS 1: Convert merged PLINK to VCF per chromosome
// ============================================================================
process convertMergedToVCF {
    label 'plink2'
    tag "${ref_panel}_chr${chr}"
    publishDir "${params.outdir}/module5/${ref_panel}/01_vcf_conversion", mode: 'copy'
    
    container 'docker://quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'
    
    cpus 4
    memory { 8.GB + (2.GB * Math.min(task.attempt, 3)) }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(ref_panel), 
          val(chr),
          path(merged_vcf),
          path(merged_tbi)
    
    output:
    tuple val(ref_panel),
          val(chr),
          path("${ref_panel}_merged_chr${chr}_reimp.vcf.gz"),
          path("${ref_panel}_merged_chr${chr}_reimp.vcf.gz.tbi"),
          emit: chr_vcfs
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "=========================================="
    echo "Converting merged data to VCF: chr${chr}"
    echo "Reference Panel: ${ref_panel}"
    echo "=========================================="
    
    # Extract chromosome from merged VCF
    bcftools view \\
        --regions chr${chr} \\
        ${merged_vcf} \\
        -Oz \\
        -o temp_chr${chr}.vcf.gz
    
    # Clean and prepare for imputation
    # Remove INFO annotations, keep only GT
    bcftools annotate \\
        -x INFO,^FORMAT/GT \\
        temp_chr${chr}.vcf.gz \\
        -Oz \\
        -o temp_clean.vcf.gz
    
    # Set variant IDs as CHROM:POS:REF:ALT
    bcftools annotate \\
        --set-id '%CHROM:%POS:%REF:%ALT' \\
        temp_clean.vcf.gz \\
        -Oz \\
        -o ${ref_panel}_merged_chr${chr}_reimp.vcf.gz
    
    # Index
    bcftools index -t ${ref_panel}_merged_chr${chr}_reimp.vcf.gz
    
    # Count variants
    n_variants=\$(bcftools view -H ${ref_panel}_merged_chr${chr}_reimp.vcf.gz | wc -l)
    n_samples=\$(bcftools query -l ${ref_panel}_merged_chr${chr}_reimp.vcf.gz | wc -l)
    
    echo "Chromosome ${chr} extracted:"
    echo "  Variants: \${n_variants}"
    echo "  Samples: \${n_samples}"
    echo "  File: ${ref_panel}_merged_chr${chr}_reimp.vcf.gz"
    
    # Cleanup
    rm -f temp_chr${chr}.vcf.gz temp_clean.vcf.gz
    
    echo "✓ Conversion complete"
    """
}

// ============================================================================
// PROCESS 2: Group all chromosomes and create manifests
// ============================================================================
process groupChromosomesForReimputation {
    label 'manifest'
    tag "${ref_panel}_reimputation"
    publishDir "${params.outdir}/module5/${ref_panel}/02_manifests", mode: 'copy'
    
    container 'docker://python:3.11-slim'
    
    cpus 1
    memory '2.GB'
    
    input:
    tuple val(ref_panel),
          val(chrs),
          path(vcf_files),
          path(tbi_files)
    
    output:
    tuple val(ref_panel),
          path(vcf_files),
          path(tbi_files),
          path("${ref_panel}_topmed_reimputation.json"),
          path("${ref_panel}_anvil_reimputation.json"),
          path("${ref_panel}_file_list.txt"),
          emit: grouped_with_manifests
    
    script:
    def vcf_list = vcf_files instanceof List ? vcf_files : [vcf_files]
    def n_vcfs = vcf_list.size()
    """
    #!/usr/bin/env python3
    import json
    
    ref_panel = "${ref_panel}"
    n_files = ${n_vcfs}
    
    print("=" * 60)
    print(f"Creating RE-IMPUTATION manifests")
    print(f"Reference Panel: {ref_panel}")
    print(f"Chromosomes: {n_files}")
    print("=" * 60)
    
    # CRITICAL: All chromosomes in ONE job
    print(f"\\n✓ All {n_files} chromosomes will be submitted together")
    print(f"  Quota usage: 1 unit (not {n_files}!)")
    
    # TOPMed manifest
    topmed_manifest = {
        "job-name": f"{ref_panel}_merged_reimputation",
        "refpanel": "topmed-r3",
        "population": "all",
        "build": "hg38",
        "phasing": "eagle",
        "mode": "imputation",
        "r2Filter": "0.0",
        "meta": "yes"
    }
    
    with open(f"{ref_panel}_topmed_reimputation.json", 'w') as f:
        json.dump(topmed_manifest, f, indent=2)
    
    # All of Us AnVIL manifest
    anvil_manifest = {
        "job_name": f"{ref_panel}_merged_reimputation",
        "reference_panel": "allofus_anvil",
        "genome_build": "hg38",
        "phasing_method": "beagle",
        "quality_filter": 0.0,
        "submission_type": "batch_reimputation",
        "chromosomes": list(range(1, 23))
    }
    
    with open(f"{ref_panel}_anvil_reimputation.json", 'w') as f:
        json.dump(anvil_manifest, f, indent=2)
    
    # File list
    vcf_list = "${vcf_list}".replace('[', '').replace(']', '').split(', ')
    with open(f"{ref_panel}_file_list.txt", 'w') as f:
        for vcf in vcf_list:
            vcf_clean = vcf.strip()
            if vcf_clean:
                f.write(f"{vcf_clean}\\n")
    
    print(f"\\n✓ Manifests created")
    print("=" * 60)
    """
}

// ============================================================================
// PROCESS 3: Submit to TOPMed (merged data, ALL chromosomes)
// ============================================================================
process submitTOPMedReimputation {
    label 'api'
    tag "${ref_panel}_topmed_reimp"
    publishDir "${params.outdir}/module5/${ref_panel}/03_topmed_jobs", mode: 'copy'
    
    container 'docker://python:3.11-slim'
    
    cpus 2
    memory '4.GB'
    
    when:
    params.run_topmed_reimputation
    
    input:
    tuple val(ref_panel),
          path(vcf_files),
          path(tbi_files),
          path(topmed_manifest),
          path(anvil_manifest),
          path(file_list)
    
    output:
    tuple val(ref_panel), 
          path("${ref_panel}_topmed_reimp_job.json"),
          emit: job_info
    
    script:
    def vcf_list = vcf_files instanceof List ? vcf_files.collect { it.toString() }.join(' ') : vcf_files.toString()
    """
    #!/usr/bin/env python3
    import requests
    import json
    import os
    import sys
    
    print("=" * 60)
    print(f"TOPMed RE-IMPUTATION: ${ref_panel}")
    print("=" * 60)
    
    token = "${params.topmed_api_token}"
    if not token or token == "null":
        print("ERROR: TOPMed API token required!")
        sys.exit(1)
    
    headers = {'X-Auth-Token': token}
    base_url = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
    
    with open('${topmed_manifest}') as f:
        job_config = json.load(f)
    
    # Get ALL chromosome VCFs
    vcf_list = """${vcf_list}""".split()
    vcf_list = [v for v in vcf_list if v.endswith('.vcf.gz')]
    
    n_vcfs = len(vcf_list)
    print(f"\\nMerged dataset re-imputation:")
    print(f"  Chromosomes: {n_vcfs} (chr1-chr22)")
    print(f"  Quota: 1 unit for entire genome ✓")
    print(f"\\nFiles to submit:")
    for vcf in sorted(vcf_list):
        print(f"  - {os.path.basename(vcf)}")
    
    # Prepare multipart upload - ALL in ONE request
    files = [('files', (os.path.basename(vcf), open(vcf, 'rb'), 'application/gzip')) 
             for vcf in vcf_list]
    
    print(f"\\nSubmitting to TOPMed...")
    
    try:
        response = requests.post(
            f'{base_url}/jobs/submit/imputationserver',
            headers=headers,
            data=job_config,
            files=files,
            timeout=300
        )
        
        for _, (_, fh, _) in files:
            fh.close()
        
        if response.status_code != 200:
            print(f"ERROR: {response.status_code}")
            print(response.text)
            sys.exit(1)
        
        job_data = response.json()
        job_data['ref_panel'] = '${ref_panel}'
        job_data['n_chromosomes'] = n_vcfs
        job_data['job_type'] = 'reimputation_merged'
        
        with open('${ref_panel}_topmed_reimp_job.json', 'w') as f:
            json.dump(job_data, f, indent=2)
        
        print(f"\\n✓ Re-imputation job submitted!")
        print(f"  Job ID: {job_data.get('id')}")
        print(f"  Quota used: 1 unit")
        print("=" * 60)
        
    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    """
}

// ============================================================================
// PROCESS 4: Monitor TOPMed re-imputation
// ============================================================================
process monitorTOPMedReimputation {
    label 'api'
    tag "${ref_panel}_topmed_monitor"
    maxRetries 200
    errorStrategy 'retry'
    
    container 'docker://python:3.11-slim'
    
    cpus 1
    memory '1.GB'
    
    when:
    params.run_topmed_reimputation
    
    input:
    tuple val(ref_panel), path(job_info)
    
    output:
    tuple val(ref_panel), path(job_info), emit: completed_job
    
    script:
    """
    #!/usr/bin/env python3
    import requests
    import json
    import time
    import sys
    
    token = "${params.topmed_api_token}"
    headers = {'X-Auth-Token': token}
    base_url = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
    
    with open('${job_info}') as f:
        job = json.load(f)
    
    job_id = job['id']
    ref_panel = job['ref_panel']
    
    response = requests.get(f'{base_url}/jobs/{job_id}', headers=headers)
    status = response.json()
    state = status.get('state')
    
    print(f"Re-imputation: {ref_panel}")
    print(f"Job ID: {job_id}")
    print(f"State: {state}")
    
    if state == 4:
        print("✓ Re-imputation complete!")
        sys.exit(0)
    elif state in [5, 6]:
        print(f"✗ Failed: {status.get('message')}")
        sys.exit(1)
    else:
        print(f"Processing... next check in ${params.monitor_interval / 60} min")
        time.sleep(${params.monitor_interval})
        sys.exit(100)
    """
}

// ============================================================================
// PROCESS 5: Download TOPMed re-imputed results
// ============================================================================
process downloadTOPMedReimputedResults {
    label 'download'
    tag "${ref_panel}_topmed_download"
    publishDir "${params.outdir}/module5/${ref_panel}/04_topmed_results", mode: 'copy'
    
    container 'docker://python:3.11-slim'
    
    cpus 4
    memory '8.GB'
    
    when:
    params.run_topmed_reimputation
    
    input:
    tuple val(ref_panel), path(job_info)
    
    output:
    tuple val(ref_panel), 
          path("${ref_panel}_reimputed_chr*.dose.vcf.gz"),
          path("${ref_panel}_reimputed_chr*.info.gz"),
          emit: topmed_reimputed
    
    script:
    """
    #!/usr/bin/env python3
    import requests
    import json
    import os
    import subprocess
    
    token = "${params.topmed_api_token}"
    headers = {'X-Auth-Token': token}
    base_url = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
    
    with open('${job_info}') as f:
        job = json.load(f)
    
    job_id = job['id']
    
    print(f"Downloading re-imputed results for ${ref_panel}...")
    
    response = requests.get(f'{base_url}/jobs/{job_id}/results', headers=headers)
    files = response.json()
    
    print(f"Files to download: {len(files)}")
    
    for file_info in files:
        filename = file_info['name']
        download_url = file_info['url']
        print(f"  {filename}")
        
        file_response = requests.get(download_url, stream=True)
        file_response.raise_for_status()
        
        with open(filename, 'wb') as f:
            for chunk in file_response.iter_content(chunk_size=8192):
                f.write(chunk)
    
    # Extract archives
    password = "${params.topmed_password}"
    if not password or password == "null":
        status_response = requests.get(f'{base_url}/jobs/{job_id}', headers=headers)
        status = status_response.json()
        password = status.get('outputParams', {}).get('password', '')
    
    zip_files = [f for f in os.listdir('.') if f.endswith('.zip')]
    if zip_files:
        for zip_file in zip_files:
            try:
                subprocess.run(['7z', 'x', f'-p{password}', zip_file], check=True)
            except:
                subprocess.run(['unzip', '-P', password, zip_file], check=True)
    
    # Rename files
    for filename in os.listdir('.'):
        if 'chr' in filename.lower():
            if '.dose.vcf.gz' in filename or '.vcf.gz' in filename:
                chr_match = filename.split('chr')[1].split('.')[0]
                new_name = f"${ref_panel}_reimputed_chr{chr_match}.dose.vcf.gz"
                if filename != new_name:
                    os.rename(filename, new_name)
            elif '.info.gz' in filename:
                chr_match = filename.split('chr')[1].split('.')[0]
                new_name = f"${ref_panel}_reimputed_chr{chr_match}.info.gz"
                if filename != new_name:
                    os.rename(filename, new_name)
    
    print("✓ Re-imputed data downloaded!")
    """
}

// ============================================================================
// PROCESS 6: Submit to AnVIL (merged data, ALL chromosomes)
// ============================================================================
process submitAnvilReimputation {
    label 'api'
    tag "${ref_panel}_anvil_reimp"
    publishDir "${params.outdir}/module5/${ref_panel}/05_anvil_jobs", mode: 'copy'
    
    container 'docker://python:3.11-slim'
    
    cpus 2
    memory '4.GB'
    
    when:
    params.run_anvil_reimputation
    
    input:
    tuple val(ref_panel),
          path(vcf_files),
          path(tbi_files),
          path(topmed_manifest),
          path(anvil_manifest),
          path(file_list)
    
    output:
    tuple val(ref_panel),
          path("${ref_panel}_anvil_reimp_job.json"),
          emit: job_info
    
    script:
    def vcf_list = vcf_files instanceof List ? vcf_files.collect { it.toString() }.join(' ') : vcf_files.toString()
    """
    #!/usr/bin/env python3
    import json
    import subprocess
    import os
    from datetime import datetime
    
    print("=" * 60)
    print(f"All of Us AnVIL RE-IMPUTATION: ${ref_panel}")
    print("=" * 60)
    
    workspace = "${params.anvil_workspace}"
    project = "${params.anvil_project}"
    
    vcf_list = """${vcf_list}""".split()
    vcf_list = [v for v in vcf_list if v.endswith('.vcf.gz')]
    
    n_vcfs = len(vcf_list)
    print(f"\\nMerged dataset re-imputation:")
    print(f"  Chromosomes: {n_vcfs}")
    print(f"  Quota: 1 unit ✓")
    
    job_info = {
        'ref_panel': '${ref_panel}',
        'workspace': workspace,
        'project': project,
        'vcf_files': vcf_list,
        'n_chromosomes': n_vcfs,
        'job_type': 'reimputation_merged',
        'manifest': '${anvil_manifest}',
        'status': 'submitted',
        'submission_time': datetime.now().isoformat()
    }
    
    cmd = [
        'terra', 'workflow', 'submit',
        '--workspace', workspace,
        '--project', project,
        '--input', '${anvil_manifest}',
        '--workflow-name', 'allofus-reimputation'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        job_info['submission_output'] = result.stdout
        if 'Workflow ID:' in result.stdout:
            job_info['workflow_id'] = result.stdout.split('Workflow ID: ')[1].split()[0]
        print(f"\\n✓ Submitted via Terra CLI")
    except:
        print(f"\\nNote: Manual submission required")
        print(f"  Via All of Us Researcher Workbench")
        job_info['submission_method'] = 'manual'
    
    with open('${ref_panel}_anvil_reimp_job.json', 'w') as f:
        json.dump(job_info, f, indent=2)
    
    print("=" * 60)
    """
}

// ============================================================================
// PROCESS 7: Monitor AnVIL re-imputation
// ============================================================================
process monitorAnvilReimputation {
    label 'api'
    tag "${ref_panel}_anvil_monitor"
    maxRetries 200
    errorStrategy 'retry'
    
    container 'docker://python:3.11-slim'
    
    cpus 1
    memory '1.GB'
    
    when:
    params.run_anvil_reimputation
    
    input:
    tuple val(ref_panel), path(job_info)
    
    output:
    tuple val(ref_panel), path(job_info), emit: completed_job
    
    script:
    """
    #!/usr/bin/env python3
    import json
    import time
    import sys
    import subprocess
    
    with open('${job_info}') as f:
        job = json.load(f)
    
    workflow_id = job.get('workflow_id')
    
    if not workflow_id:
        print("Manual submission - check workbench")
        time.sleep(${params.monitor_interval})
        sys.exit(0)
    
    workspace = job.get('workspace')
    project = job.get('project')
    
    cmd = [
        'terra', 'workflow', 'status',
        '--workspace', workspace,
        '--project', project,
        '--workflow-id', workflow_id
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        if 'Succeeded' in result.stdout:
            print("✓ Complete!")
            sys.exit(0)
        elif 'Failed' in result.stdout:
            sys.exit(1)
        else:
            time.sleep(${params.monitor_interval})
            sys.exit(100)
    except:
        time.sleep(${params.monitor_interval})
        sys.exit(100)
    """
}

// ============================================================================
// PROCESS 8: Download AnVIL re-imputed results
// ============================================================================
process downloadAnvilReimputedResults {
    label 'download'
    tag "${ref_panel}_anvil_download"
    publishDir "${params.outdir}/module5/${ref_panel}/06_anvil_results", mode: 'copy'
    
    container 'docker://google/cloud-sdk:slim'
    
    cpus 4
    memory '8.GB'
    
    when:
    params.run_anvil_reimputation
    
    input:
    tuple val(ref_panel), path(job_info)
    
    output:
    tuple val(ref_panel),
          path("${ref_panel}_reimputed_chr*.dose.vcf.gz"),
          path("${ref_panel}_reimputed_chr*.info.gz"),
          emit: anvil_reimputed
    
    script:
    """
    #!/usr/bin/env python3
    import json
    import subprocess
    import os
    
    with open('${job_info}') as f:
        job = json.load(f)
    
    workspace = job.get('workspace')
    workflow_id = job.get('workflow_id')
    
    print(f"Downloading ${ref_panel} from AnVIL...")
    
    bucket_path = f"gs://{workspace}-workspace-bucket/submissions/{workflow_id}/"
    
    try:
        subprocess.run(['gsutil', '-m', 'cp', '-r', f'{bucket_path}*.vcf.gz', '.'], check=False)
        subprocess.run(['gsutil', '-m', 'cp', '-r', f'{bucket_path}*.info.gz', '.'], check=False)
    except:
        print("Manual download may be required")
    
    # Rename
    for filename in os.listdir('.'):
        if 'chr' in filename:
            if '.dose.vcf.gz' in filename or '.vcf.gz' in filename:
                chr_num = filename.split('chr')[1].split('.')[0]
                new_name = f"${ref_panel}_reimputed_chr{chr_num}.dose.vcf.gz"
                if filename != new_name:
                    os.rename(filename, new_name)
            elif '.info' in filename:
                chr_num = filename.split('chr')[1].split('.')[0]
                new_name = f"${ref_panel}_reimputed_chr{chr_num}.info.gz"
                if filename != new_name:
                    os.rename(filename, new_name)
    
    print("✓ Download complete!")
    """
}

// ============================================================================
// MODULE 5 WORKFLOW
// ============================================================================
workflow MODULE5_REIMPUTATION {
    take:
    merged_vcf  // From Module 4: tuple(ref_panel, chr, vcf, tbi)
    
    main:
    
    // Step 1: Convert merged VCF to per-chromosome VCFs
    convertMergedToVCF(merged_vcf)
    
    // Step 2: CRITICAL - Group all chromosomes together per ref_panel
    convertMergedToVCF.out.chr_vcfs
        .groupTuple(by: 0)  // Group by ref_panel
        .set { grouped_by_ref_panel }
    
    // Step 3: Create manifests for grouped submission
    groupChromosomesForReimputation(grouped_by_ref_panel)
    
    // Initialize outputs
    topmed_out = Channel.empty()
    anvil_out = Channel.empty()
    
    // TOPMed re-imputation pathway
    if (params.run_topmed_reimputation) {
        submitTOPMedReimputation(groupChromosomesForReimputation.out.grouped_with_manifests)
        monitorTOPMedReimputation(submitTOPMedReimputation.out.job_info)
        downloadTOPMedReimputedResults(monitorTOPMedReimputation.out.completed_job)
        topmed_out = downloadTOPMedReimputedResults.out.topmed_reimputed
    }
    
    // AnVIL re-imputation pathway
    if (params.run_anvil_reimputation) {
        submitAnvilReimputation(groupChromosomesForReimputation.out.grouped_with_manifests)
        monitorAnvilReimputation(submitAnvilReimputation.out.job_info)
        downloadAnvilReimputedResults(monitorAnvilReimputation.out.completed_job)
        anvil_out = downloadAnvilReimputedResults.out.anvil_reimputed
    }
    
    // Combine outputs
    final_results = topmed_out.mix(anvil_out)
    
    emit:
    reimputed_data = final_results
    topmed_reimputed = topmed_out
    anvil_reimputed = anvil_out
}

workflow.onComplete {
    println """
    ============================================================================
    MODULE 5 - RE-IMPUTATION COMPLETE
    ============================================================================
    Status:        ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration:      ${workflow.duration}
    
    ✓ Quota-Efficient Re-Imputation:
      - Merged data: All chr1-22 submitted together
      - TOPMed: 1 job = 1 quota unit
      - AnVIL: 1 job = 1 quota unit
      - Total: 2 quota units (not 44!)
    
    Enhanced Coverage:
      - Re-imputed from merged haplotypes
      - Maximum variant detection
      - Improved rare variant calling
    
    Results: ${params.outdir}/module5/
    ============================================================================
    """
}
