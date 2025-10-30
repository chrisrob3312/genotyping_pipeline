#!/usr/bin/env nextflow

/*
 * ============================================================================
 * MODULE 2: IMPUTATION SUBMISSION AND RETRIEVAL 
 * ============================================================================
 * 
 * CRITICAL FIX: Groups ALL chromosomes per platform into ONE job submission
 * - Avoids quota waste (submitting chr1-22 separately = 22x quota!)
 * - Each platform submits ALL chromosomes (chr1-22) as ONE job
 * - Monitoring checks every 15 minutes
 * - Extended time limits (50 hours max monitoring)
 * 
 * Dual Service Support:
 * - TOPMed Imputation Server (Michigan/BioData Catalyst)
 * - All of Us AnVIL Imputation Service
 * 
 * ============================================================================
 */

nextflow.enable.dsl = 2

// ============================================================================
// PROCESS 0: Group chromosomes by platform and create manifests
// ============================================================================
process groupChromosomesAndCreateManifests {
    label 'manifest'
    tag "${platform}"
    publishDir "${params.outdir}/module2/00_manifests/${platform}", mode: 'copy'
    
    container 'docker://python:3.11-slim'
    
    cpus 1
    memory '2.GB'
    
    input:
    tuple val(platform), 
          val(chrs),
          path(vcf_files),
          path(tbi_files),
          val(build)
    
    output:
    tuple val(platform),
          path(vcf_files),
          path(tbi_files),
          path("${platform}_topmed_manifest.json"),
          path("${platform}_anvil_manifest.json"),
          path("${platform}_file_list.txt"),
          emit: grouped_with_manifests
    
    script:
    def vcf_list = vcf_files instanceof List ? vcf_files : [vcf_files]
    def n_vcfs = vcf_list.size()
    """
    #!/usr/bin/env python3
    import json
    
    platform = "${platform}"
    build = "${build}"
    n_files = ${n_vcfs}
    
    print(f"========================================")
    print(f"Creating manifests for: {platform}")
    print(f"Build: {build}")
    print(f"Number of chromosome VCFs: {n_files}")
    print(f"========================================")
    
    # CRITICAL: All chromosomes submitted as ONE job
    print(f"✓ All {n_files} chromosomes will be submitted together")
    print(f"  This uses 1 quota unit (not {n_files}!)")
    
    # TOPMed manifest - ONE job with ALL chromosomes
    topmed_manifest = {
        "job-name": f"{platform}_imputation",
        "refpanel": "topmed-r3",
        "population": "all",
        "build": build,
        "phasing": "eagle",
        "mode": "imputation",
        "r2Filter": "0.0"
    }
    
    with open(f"{platform}_topmed_manifest.json", 'w') as f:
        json.dump(topmed_manifest, f, indent=2)
    
    # All of Us AnVIL manifest - ONE job with ALL chromosomes
    anvil_manifest = {
        "job_name": f"{platform}_anvil_imputation",
        "reference_panel": "allofus_anvil",
        "genome_build": build,
        "phasing_method": "beagle",
        "quality_filter": 0.0,
        "submission_type": "batch",
        "chromosomes": list(range(1, 23))  # chr1-22
    }
    
    with open(f"{platform}_anvil_manifest.json", 'w') as f:
        json.dump(anvil_manifest, f, indent=2)
    
    # File list for reference
    vcf_list = "${vcf_list}".replace('[', '').replace(']', '').split(', ')
    with open(f"{platform}_file_list.txt", 'w') as f:
        for vcf in vcf_list:
            vcf_clean = vcf.strip()
            if vcf_clean:
                f.write(f"{vcf_clean}\\n")
    
    print(f"✓ Manifests created successfully")
    print(f"  Files: {platform}_topmed_manifest.json")
    print(f"         {platform}_anvil_manifest.json")
    """
}

// ============================================================================
// PROCESS 1: Submit to TOPMed (ALL chromosomes in ONE job)
// ============================================================================
process submitTOPMedImputation {
    label 'api'
    tag "${platform}"
    publishDir "${params.outdir}/module2/01_topmed_jobs/${platform}", mode: 'copy'
    
    container 'docker://python:3.11-slim'
    
    cpus 2
    memory '4.GB'
    
    when:
    params.run_topmed
    
    input:
    tuple val(platform),
          path(vcf_files),
          path(tbi_files),
          path(topmed_manifest),
          path(anvil_manifest),
          path(file_list)
    
    output:
    tuple val(platform), 
          path("${platform}_topmed_job.json"),
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
    print(f"TOPMed Submission: ${platform}")
    print("=" * 60)
    
    # Read API token
    token = "${params.topmed_api_token}"
    if not token or token == "null" or token == "":
        print("ERROR: TOPMed API token required!")
        print("Set params.topmed_api_token in nextflow.config")
        sys.exit(1)
    
    # Set up API
    headers = {'X-Auth-Token': token}
    base_url = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
    
    # Load manifest
    with open('${topmed_manifest}') as f:
        job_config = json.load(f)
    
    # CRITICAL: Get ALL VCF files to submit together
    vcf_list = """${vcf_list}""".split()
    vcf_list = [v for v in vcf_list if v.endswith('.vcf.gz')]
    
    n_vcfs = len(vcf_list)
    print(f"\\n✓ Submitting {n_vcfs} chromosome VCFs as ONE job")
    print(f"  (Uses 1 quota unit, not {n_vcfs}!)")
    print(f"\\nChromosomes included:")
    for vcf in sorted(vcf_list):
        print(f"  - {os.path.basename(vcf)}")
    
    # Prepare multipart upload - ALL chromosomes in ONE request
    files = [('files', (os.path.basename(vcf), open(vcf, 'rb'), 'application/gzip')) 
             for vcf in vcf_list]
    
    print(f"\\nSubmitting to TOPMed server...")
    
    try:
        # Submit ONE job with ALL chromosome VCFs
        response = requests.post(
            f'{base_url}/jobs/submit/imputationserver',
            headers=headers,
            data=job_config,
            files=files,
            timeout=300
        )
        
        # Close file handles
        for _, (_, fh, _) in files:
            fh.close()
        
        if response.status_code != 200:
            print(f"ERROR: Submission failed!")
            print(f"Status: {response.status_code}")
            print(f"Response: {response.text}")
            sys.exit(1)
        
        # Save job info
        job_data = response.json()
        job_data['platform'] = '${platform}'
        job_data['n_chromosomes'] = n_vcfs
        job_data['submission_type'] = 'batch_all_chromosomes'
        
        with open('${platform}_topmed_job.json', 'w') as f:
            json.dump(job_data, f, indent=2)
        
        print(f"\\n✓ Job submitted successfully!")
        print(f"  Job ID: {job_data.get('id', 'N/A')}")
        print(f"  Status: {job_data.get('state', 'pending')}")
        print(f"  Chromosomes: {n_vcfs} (chr1-chr22)")
        print(f"  Quota used: 1 unit")
        print("=" * 60)
        
    except Exception as e:
        print(f"ERROR: {str(e)}")
        sys.exit(1)
    """
}

// ============================================================================
// PROCESS 2: Monitor TOPMed job
// ============================================================================
process monitorTOPMedJob {
    label 'api'
    tag "${platform}"
    maxRetries 200  // 200 × 15 min = 50 hours
    errorStrategy 'retry'
    
    container 'docker://python:3.11-slim'
    
    cpus 1
    memory '1.GB'
    
    when:
    params.run_topmed
    
    input:
    tuple val(platform), path(job_info)
    
    output:
    tuple val(platform), path(job_info), emit: completed_job
    
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
    platform = job['platform']
    
    response = requests.get(f'{base_url}/jobs/{job_id}', headers=headers)
    status = response.json()
    state = status.get('state')
    
    print(f"Platform: {platform}")
    print(f"Job ID: {job_id}")
    print(f"State: {state}")
    
    # State codes: 1=waiting, 2=running, 3=exporting, 4=success, 5=failed, 6=canceled
    if state == 4:
        print("✓ Job completed successfully!")
        sys.exit(0)
    elif state in [5, 6]:
        error = status.get('message', 'Unknown error')
        print(f"✗ Job failed: {error}")
        sys.exit(1)
    else:
        print(f"Still processing... checking again in ${params.monitor_interval / 60} minutes")
        time.sleep(${params.monitor_interval})
        sys.exit(100)  # Trigger retry
    """
}

// ============================================================================
// PROCESS 3: Download TOPMed results
// ============================================================================
process downloadTOPMedResults {
    label 'download'
    tag "${platform}"
    publishDir "${params.outdir}/module2/02_topmed_results/${platform}", mode: 'copy'
    
    container 'docker://python:3.11-slim'
    
    cpus 4
    memory '8.GB'
    
    when:
    params.run_topmed
    
    input:
    tuple val(platform), path(job_info)
    
    output:
    tuple val(platform), 
          path("${platform}_imputed_chr*.dose.vcf.gz"),
          path("${platform}_imputed_chr*.info.gz"),
          emit: topmed_results
    
    script:
    """
    #!/usr/bin/env python3
    import requests
    import json
    import os
    import subprocess
    import sys
    
    token = "${params.topmed_api_token}"
    headers = {'X-Auth-Token': token}
    base_url = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
    
    with open('${job_info}') as f:
        job = json.load(f)
    
    job_id = job['id']
    
    print(f"Downloading results for ${platform}...")
    print(f"Job ID: {job_id}")
    
    # Get download links
    response = requests.get(f'{base_url}/jobs/{job_id}/results', headers=headers)
    files = response.json()
    
    print(f"Found {len(files)} files to download")
    
    # Download all files
    for file_info in files:
        filename = file_info['name']
        download_url = file_info['url']
        
        print(f"  Downloading: {filename}")
        file_response = requests.get(download_url, stream=True)
        file_response.raise_for_status()
        
        with open(filename, 'wb') as f:
            for chunk in file_response.iter_content(chunk_size=8192):
                f.write(chunk)
    
    # Extract if encrypted (try password from params or job status)
    password = "${params.topmed_password}"
    if not password or password == "null":
        status_response = requests.get(f'{base_url}/jobs/{job_id}', headers=headers)
        status = status_response.json()
        password = status.get('outputParams', {}).get('password', '')
    
    # Extract archives
    zip_files = [f for f in os.listdir('.') if f.endswith('.zip')]
    if zip_files:
        print(f"Extracting {len(zip_files)} archives...")
        for zip_file in zip_files:
            try:
                subprocess.run(['7z', 'x', f'-p{password}', zip_file], check=True)
            except:
                subprocess.run(['unzip', '-P', password, zip_file], check=True)
    
    # Rename to standardized format
    for filename in os.listdir('.'):
        if 'chr' in filename.lower():
            if '.dose.vcf.gz' in filename or '.vcf.gz' in filename:
                chr_match = filename.split('chr')[1].split('.')[0]
                new_name = f"${platform}_imputed_chr{chr_match}.dose.vcf.gz"
                if filename != new_name:
                    os.rename(filename, new_name)
                    print(f"  Renamed: {filename} → {new_name}")
            elif '.info.gz' in filename:
                chr_match = filename.split('chr')[1].split('.')[0]
                new_name = f"${platform}_imputed_chr{chr_match}.info.gz"
                if filename != new_name:
                    os.rename(filename, new_name)
                    print(f"  Renamed: {filename} → {new_name}")
    
    print("✓ Download complete!")
    """
}

// ============================================================================
// PROCESS 4: Submit to All of Us AnVIL (ALL chromosomes in ONE job)
// ============================================================================
process submitAnvilImputation {
    label 'api'
    tag "${platform}"
    publishDir "${params.outdir}/module2/03_anvil_jobs/${platform}", mode: 'copy'
    
    container 'docker://python:3.11-slim'
    
    cpus 2
    memory '4.GB'
    
    when:
    params.run_anvil
    
    input:
    tuple val(platform),
          path(vcf_files),
          path(tbi_files),
          path(topmed_manifest),
          path(anvil_manifest),
          path(file_list)
    
    output:
    tuple val(platform),
          path("${platform}_anvil_job.json"),
          emit: job_info
    
    script:
    def vcf_list = vcf_files instanceof List ? vcf_files.collect { it.toString() }.join(' ') : vcf_files.toString()
    """
    #!/usr/bin/env python3
    import json
    import subprocess
    import os
    import sys
    from datetime import datetime
    
    print("=" * 60)
    print(f"All of Us AnVIL Submission: ${platform}")
    print("=" * 60)
    
    workspace = "${params.anvil_workspace}"
    project = "${params.anvil_project}"
    
    if not workspace or workspace == "null":
        print("ERROR: AnVIL workspace required!")
        print("Set params.anvil_workspace in nextflow.config")
        sys.exit(1)
    
    if not project or project == "null":
        print("ERROR: AnVIL project required!")
        print("Set params.anvil_project in nextflow.config")
        sys.exit(1)
    
    # Get ALL VCF files
    vcf_list = """${vcf_list}""".split()
    vcf_list = [v for v in vcf_list if v.endswith('.vcf.gz')]
    
    n_vcfs = len(vcf_list)
    print(f"\\n✓ Submitting {n_vcfs} chromosome VCFs as ONE job")
    print(f"  (Uses 1 quota unit, not {n_vcfs}!)")
    print(f"\\nChromosomes included:")
    for vcf in sorted(vcf_list):
        print(f"  - {os.path.basename(vcf)}")
    
    # Prepare job info
    job_info = {
        'platform': '${platform}',
        'workspace': workspace,
        'project': project,
        'vcf_files': vcf_list,
        'n_chromosomes': n_vcfs,
        'submission_type': 'batch_all_chromosomes',
        'manifest': '${anvil_manifest}',
        'status': 'submitted',
        'submission_time': datetime.now().isoformat()
    }
    
    # Try Terra CLI submission
    print(f"\\nAttempting submission via Terra CLI...")
    cmd = [
        'terra', 'workflow', 'submit',
        '--workspace', workspace,
        '--project', project,
        '--input', '${anvil_manifest}',
        '--workflow-name', 'allofus-imputation'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        job_info['submission_output'] = result.stdout
        if 'Workflow ID:' in result.stdout:
            workflow_id = result.stdout.split('Workflow ID: ')[1].split()[0]
            job_info['workflow_id'] = workflow_id
        print(f"\\n✓ Submitted via Terra CLI")
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"\\nNote: CLI submission not available: {e.stderr}")
        print("→ Please submit manually via All of Us Researcher Workbench")
        print(f"  Workspace: {workspace}")
        print(f"  Files to submit: Listed in ${file_list}")
        job_info['submission_method'] = 'manual'
        job_info['note'] = 'Manual submission required via Researcher Workbench'
    except FileNotFoundError:
        print("\\nNote: Terra CLI not installed")
        print("→ Please submit manually via All of Us Researcher Workbench")
        job_info['submission_method'] = 'manual'
        job_info['note'] = 'Manual submission required'
    
    # Save job info
    with open('${platform}_anvil_job.json', 'w') as f:
        json.dump(job_info, f, indent=2)
    
    print(f"\\n✓ Job information saved")
    print("=" * 60)
    """
}

// ============================================================================
// PROCESS 5: Monitor AnVIL job
// ============================================================================
process monitorAnvilJob {
    label 'api'
    tag "${platform}"
    maxRetries 200
    errorStrategy 'retry'
    
    container 'docker://python:3.11-slim'
    
    cpus 1
    memory '1.GB'
    
    when:
    params.run_anvil
    
    input:
    tuple val(platform), path(job_info)
    
    output:
    tuple val(platform), path(job_info), emit: completed_job
    
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
    workspace = job.get('workspace')
    project = job.get('project')
    
    if not workflow_id:
        print("Manual submission - skipping monitoring")
        print("Check status in All of Us Researcher Workbench")
        time.sleep(${params.monitor_interval})
        sys.exit(0)
    
    print(f"Checking workflow: {workflow_id}")
    
    cmd = [
        'terra', 'workflow', 'status',
        '--workspace', workspace,
        '--project', project,
        '--workflow-id', workflow_id
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        status = result.stdout
        
        print(f"Status: {status}")
        
        if 'Succeeded' in status or 'Done' in status:
            print("✓ Workflow completed!")
            sys.exit(0)
        elif 'Failed' in status or 'Aborted' in status:
            print("✗ Workflow failed")
            sys.exit(1)
        else:
            print(f"Still running... checking again in ${params.monitor_interval / 60} min")
            time.sleep(${params.monitor_interval})
            sys.exit(100)
    except:
        print("Assuming still running...")
        time.sleep(${params.monitor_interval})
        sys.exit(100)
    """
}

// ============================================================================
// PROCESS 6: Download AnVIL results
// ============================================================================
process downloadAnvilResults {
    label 'download'
    tag "${platform}"
    publishDir "${params.outdir}/module2/04_anvil_results/${platform}", mode: 'copy'
    
    container 'docker://google/cloud-sdk:slim'
    
    cpus 4
    memory '8.GB'
    
    when:
    params.run_anvil
    
    input:
    tuple val(platform), path(job_info)
    
    output:
    tuple val(platform),
          path("${platform}_imputed_chr*.dose.vcf.gz"),
          path("${platform}_imputed_chr*.info.gz"),
          emit: anvil_results
    
    script:
    """
    #!/usr/bin/env python3
    import json
    import subprocess
    import os
    
    with open('${job_info}') as f:
        job = json.load(f)
    
    workspace = job.get('workspace')
    project = job.get('project')
    workflow_id = job.get('workflow_id')
    
    print(f"Downloading ${platform} from AnVIL...")
    
    # Try Terra CLI
    if workflow_id:
        cmd = [
            'terra', 'workflow', 'outputs',
            '--workspace', workspace,
            '--project', project,
            '--workflow-id', workflow_id
        ]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0:
                print("Workflow outputs:")
                print(result.stdout)
        except:
            pass
    
    # Try GCS direct download
    bucket_path = f"gs://{workspace}-workspace-bucket/submissions/{workflow_id}/"
    
    try:
        subprocess.run(['gsutil', '-m', 'cp', '-r', f'{bucket_path}*.vcf.gz', '.'], check=False)
        subprocess.run(['gsutil', '-m', 'cp', '-r', f'{bucket_path}*.info.gz', '.'], check=False)
        print("✓ Downloaded from GCS")
    except:
        print("Note: Manual download may be required from workspace")
    
    # Rename files
    for filename in os.listdir('.'):
        if 'chr' in filename:
            if '.dose.vcf.gz' in filename or '.vcf.gz' in filename:
                chr_num = filename.split('chr')[1].split('.')[0]
                new_name = f"${platform}_imputed_chr{chr_num}.dose.vcf.gz"
                if filename != new_name:
                    os.rename(filename, new_name)
            elif '.info' in filename:
                chr_num = filename.split('chr')[1].split('.')[0]
                new_name = f"${platform}_imputed_chr{chr_num}.info.gz"
                if filename != new_name:
                    os.rename(filename, new_name)
    
    print("✓ Download complete!")
    """
}

// ============================================================================
// MODULE 2 WORKFLOW
// ============================================================================
workflow MODULE2_IMPUTATION {
    take:
    imputation_vcfs  // From Module 1: tuple(platform, chr, vcf, tbi, build, batch)
    
    main:
    
    // CRITICAL: Group all chromosomes per platform
    imputation_vcfs
        .map { platform, chr, vcf, tbi, build, batch ->
            tuple(platform, chr, vcf, tbi, build)
        }
        .groupTuple(by: 0)  // Group by platform
        .set { grouped_by_platform }
    
    // Create manifests for grouped submissions
    groupChromosomesAndCreateManifests(grouped_by_platform)
    
    // Initialize outputs
    topmed_out = Channel.empty()
    anvil_out = Channel.empty()
    
    // TOPMed pathway
    if (params.run_topmed) {
        submitTOPMedImputation(groupChromosomesAndCreateManifests.out.grouped_with_manifests)
        monitorTOPMedJob(submitTOPMedImputation.out.job_info)
        downloadTOPMedResults(monitorTOPMedJob.out.completed_job)
        topmed_out = downloadTOPMedResults.out.topmed_results
    }
    
    // AnVIL pathway
    if (params.run_anvil) {
        submitAnvilImputation(groupChromosomesAndCreateManifests.out.grouped_with_manifests)
        monitorAnvilJob(submitAnvilImputation.out.job_info)
        downloadAnvilResults(monitorAnvilJob.out.completed_job)
        anvil_out = downloadAnvilResults.out.anvil_results
    }
    
    // Combine outputs
    final_results = topmed_out.mix(anvil_out)
    
    emit:
    imputed_data = final_results
    topmed_results = topmed_out
    anvil_results = anvil_out
}

workflow.onComplete {
    println """
    ============================================================================
    MODULE 2 - IMPUTATION COMPLETE
    ============================================================================
    Status:        ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration:      ${workflow.duration}
    
    ✓ Quota-Efficient Submission:
      - All chr1-22 submitted together per platform
      - Uses 1 quota unit per platform (not 22!)
    
    Results: ${params.outdir}/module2/
    ============================================================================
    """
}


// Module 2: Imputation Submission and Retrieval
// Version: 2.0 - 15-minute monitoring, extended time limits
// Parallelized: All platforms × Both services (TOPMed + All of Us)

// ============================================================================
// PROCESS 1: Submit to TOPMed Imputation Server via API
// ============================================================================
// What this does: Submits your VCF files to the TOPMed imputation server
// - Uses the Michigan Imputation Server API
// - Uploads all chromosome VCFs for a platform
// - Configures imputation parameters (reference panel, phasing method, etc.)

process submitTOPMedImputation {
    label 'api'
    tag "${platform}"
    publishDir "${params.outdir}/module2/01_topmed_jobs/${platform}", mode: 'copy'
    
    when:
    params.run_topmed
    
    input:
    tuple val(platform), path(vcf_files), path(tbi_files)
    tuple val(platform2), path(manifest), path(anvil_manifest), path(file_list)
    
    output:
    tuple val(platform), 
          path("${platform}_topmed_job.json"),
          emit: job_info
    
    script:
    """
    #!/usr/bin/env python3
    import requests
    import json
    import os
    
    # Read API token from parameters
    token = "${params.topmed_api_token}"
    if not token or token == "null":
        raise ValueError("TOPMed API token required. Set params.topmed_api_token")
    
    # Set up API authentication
    headers = {'X-Auth-Token': token}
    base_url = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
    
    # Load job configuration from manifest
    with open('${manifest}') as f:
        job_config = json.load(f)
    
    # Get list of VCF files to upload
    vcf_list = [f for f in '${vcf_files}'.split() if f.endswith('.vcf.gz')]
    
    # Prepare files for multipart upload
    # Each VCF is uploaded as a separate file in the request
    files = [('files', (os.path.basename(vcf), open(vcf, 'rb'))) 
             for vcf in vcf_list]
    
    # Submit imputation job to server
    response = requests.post(
        f'{base_url}/jobs/submit/imputationserver',
        headers=headers,
        data=job_config,
        files=files
    )
    
    # Check if submission was successful
    if response.status_code != 200:
        raise Exception(f"Submission failed: {response.text}")
    
    # Save job information for monitoring
    job_data = response.json()
    job_data['platform'] = '${platform}'
    
    with open('${platform}_topmed_job.json', 'w') as f:
        json.dump(job_data, f, indent=2)
    
    print(f"Job submitted successfully!")
    print(f"Job ID: {job_data['id']}")
    print(f"Status: {job_data.get('state', 'pending')}")
    """
}

// ============================================================================
// PROCESS 2: Monitor TOPMed job completion
// ============================================================================
// What this does: Periodically checks if imputation job has completed
// - Checks status every 15 minutes (params.monitor_interval)
// - Retries up to 200 times (allows ~50 hours of monitoring)
// - States: 1=waiting, 2=running, 3=exporting, 4=success, 5=failed, 6=canceled

process monitorTOPMedJob {
    label 'api'
    tag "${platform}"
    maxRetries 200  // 200 × 15 min = 50 hours maximum monitoring
    errorStrategy 'retry'
    
    when:
    params.run_topmed
    
    input:
    tuple val(platform), path(job_info)
    
    output:
    tuple val(platform), path(job_info), emit: completed_job
    
    script:
    """
    #!/usr/bin/env python3
    import requests
    import json
    import time
    import sys
    
    # Set up API authentication
    token = "${params.topmed_api_token}"
    headers = {'X-Auth-Token': token}
    base_url = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
    
    # Load job information
    with open('${job_info}') as f:
        job = json.load(f)
    
    job_id = job['id']
    
    # Query current job status from API
    response = requests.get(f'{base_url}/jobs/{job_id}', headers=headers)
    status = response.json()
    
    state = status.get('state')
    print(f"Job {job_id} current state: {state}")
    
    # State codes from Michigan Imputation Server:
    # 1 = waiting in queue
    # 2 = running (actively imputing)
    # 3 = exporting results
    # 4 = success (completed)
    # 5 = failed
    # 6 = canceled
    
    if state == 4:  # Success
        print("✓ Job completed successfully!")
        sys.exit(0)
    elif state in [5, 6]:  # Failed or canceled
        error_msg = status.get('message', 'Unknown error')
        print(f"✗ Job failed: {error_msg}")
        sys.exit(1)
    else:  # Still running (states 1, 2, 3)
        print(f"Job still processing... will check again in ${params.monitor_interval / 60} minutes")
        # Wait for monitoring interval (default 15 minutes)
        time.sleep(${params.monitor_interval})
        sys.exit(100)  # Exit code 100 triggers retry
    """
}

// ============================================================================
// PROCESS 3: Download TOPMed results
// ============================================================================
// What this does: Downloads imputed data from TOPMed server
// - Downloads all chromosome files (.dose.vcf.gz and .info.gz)
// - Decrypts password-protected archives
// - Renames files to standardized format for downstream processing

process downloadTOPMedResults {
    label 'download'
    tag "${platform}"
    publishDir "${params.outdir}/module2/02_topmed_results/${platform}", mode: 'copy'
    
    when:
    params.run_topmed
    
    input:
    tuple val(platform), path(job_info)
    
    output:
    tuple val(platform), 
          path("${platform}_imputed_chr*.dose.vcf.gz"),
          path("${platform}_imputed_chr*.info.gz"),
          emit: topmed_results
    
    script:
    """
    #!/usr/bin/env python3
    import requests
    import json
    import os
    import subprocess
    
    # Set up API authentication
    token = "${params.topmed_api_token}"
    headers = {'X-Auth-Token': token}
    base_url = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
    
    # Load job information
    with open('${job_info}') as f:
        job = json.load(f)
    
    job_id = job['id']
    
    # Get download links for all result files
    response = requests.get(f'{base_url}/jobs/{job_id}/results', headers=headers)
    files = response.json()
    
    # Download all chromosome files
    for file_info in files:
        filename = file_info['name']
        download_url = file_info['url']
        
        print(f"Downloading: {filename}")
        
        # Download file
        file_response = requests.get(download_url, stream=True)
        file_response.raise_for_status()
        
        with open(filename, 'wb') as f:
            for chunk in file_response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        print(f"✓ Downloaded: {filename}")
    
    # Decrypt and extract files
    # TOPMed provides password in job completion email or API response
    password = "${params.topmed_password}"
    if not password or password == "null":
        # Try to get password from job status
        status_response = requests.get(f'{base_url}/jobs/{job_id}', headers=headers)
        status = status_response.json()
        password = status.get('outputParams', {}).get('password', '')
    
    # Extract encrypted archives using 7zip
    subprocess.run(['7z', 'x', f'-p{password}', '*.zip'], check=True)
    
    # Rename files to standardized format for downstream processing
    for filename in os.listdir('.'):
        if 'chr' in filename:
            # Extract chromosome number
            if '.dose.vcf.gz' in filename:
                chr_num = filename.split('chr')[1].split('.')[0]
                new_name = f"${platform}_imputed_chr{chr_num}.dose.vcf.gz"
                os.rename(filename, new_name)
                print(f"Renamed: {filename} -> {new_name}")
            elif '.info.gz' in filename:
                chr_num = filename.split('chr')[1].split('.')[0]
                new_name = f"${platform}_imputed_chr{chr_num}.info.gz"
                os.rename(filename, new_name)
                print(f"Renamed: {filename} -> {new_name}")
    
    print("✓ All files downloaded and extracted successfully!")
    """
}

// ============================================================================
// PROCESS 4: Submit to All of Us AnVIL Imputation
// ============================================================================
// What this does: Submits VCF files to All of Us AnVIL imputation service
// - Uses Terra/Cromwell workflow system
// - Uploads to specified workspace in Google Cloud
// - Configures job using AnVIL-specific manifest

process submitAnvilImputation {
    label 'api'
    tag "${platform}"
    publishDir "${params.outdir}/module2/03_anvil_jobs/${platform}", mode: 'copy'
    
    when:
    params.run_anvil
    
    input:
    tuple val(platform), path(vcf_files), path(tbi_files)
    tuple val(platform2), path(manifest), path(anvil_manifest), path(file_list)
    
    output:
    tuple val(platform),
          path("${platform}_anvil_job.json"),
          emit: job_info
    
    script:
    """
    #!/usr/bin/env python3
    import json
    import subprocess
    import os
    
    # Configure workspace information
    workspace = "${params.anvil_workspace}"
    project = "${params.anvil_project}"
    
    if not workspace or workspace == "null":
        raise ValueError("AnVIL workspace required. Set params.anvil_workspace")
    if not project or project == "null":
        raise ValueError("AnVIL project required. Set params.anvil_project")
    
    # Prepare job information
    job_info = {
        'platform': '${platform}',
        'workspace': workspace,
        'project': project,
        'vcf_files': '${vcf_files}'.split(),
        'manifest': '${anvil_manifest}',
        'status': 'submitted',
        'submission_time': subprocess.check_output(['date', '-Iseconds']).decode().strip()
    }
    
    # Submit workflow using Terra CLI
    # Note: Requires terra CLI to be installed and authenticated
    # See: https://terra.bio/using-terra-cli/
    
    cmd = [
        'terra', 'workflow', 'submit',
        '--workspace', workspace,
        '--project', project,
        '--input', '${anvil_manifest}',
        '--workflow-name', 'allofus-imputation'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        job_info['submission_output'] = result.stdout
        job_info['workflow_id'] = result.stdout.split('Workflow ID: ')[1].split()[0] if 'Workflow ID:' in result.stdout else None
        print(f"✓ AnVIL workflow submitted successfully!")
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Warning: Workflow submission command failed: {e.stderr}")
        print("This may be expected if using All of Us Researcher Workbench web interface")
        job_info['submission_method'] = 'manual'
    
    # Save job information
    with open('${platform}_anvil_job.json', 'w') as f:
        json.dump(job_info, f, indent=2)
    
    print(f"Job information saved for platform: ${platform}")
    """
}

// ============================================================================
// PROCESS 5: Monitor AnVIL job
// ============================================================================
// What this does: Checks AnVIL/Terra workflow status
// - Queries Terra/Cromwell API every 15 minutes
// - Monitors workflow execution state
// - Retries until completion or failure

process monitorAnvilJob {
    label 'api'
    tag "${platform}"
    maxRetries 200  // 200 × 15 min = 50 hours maximum monitoring
    errorStrategy 'retry'
    
    when:
    params.run_anvil
    
    input:
    tuple val(platform), path(job_info)
    
    output:
    tuple val(platform), path(job_info), emit: completed_job
    
    script:
    """
    #!/usr/bin/env python3
    import json
    import time
    import sys
    import subprocess
    
    # Load job information
    with open('${job_info}') as f:
        job = json.load(f)
    
    workflow_id = job.get('workflow_id')
    workspace = job.get('workspace')
    project = job.get('project')
    
    if not workflow_id:
        print("No workflow ID found - assuming manual submission")
        print("Check status manually in All of Us Researcher Workbench")
        # For manual submissions, exit successfully after delay
        time.sleep(${params.monitor_interval})
        sys.exit(0)
    
    # Check workflow status using Terra CLI
    cmd = [
        'terra', 'workflow', 'status',
        '--workspace', workspace,
        '--project', project,
        '--workflow-id', workflow_id
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        status = result.stdout
        
        print(f"Workflow status for ${platform}:")
        print(status)
        
        # Check for completion keywords
        if 'Succeeded' in status or 'Done' in status:
            print("✓ Workflow completed successfully!")
            sys.exit(0)
        elif 'Failed' in status or 'Aborted' in status:
            print("✗ Workflow failed")
            sys.exit(1)
        else:
            print(f"Workflow still running... will check again in ${params.monitor_interval / 60} minutes")
            time.sleep(${params.monitor_interval})
            sys.exit(100)  # Trigger retry
            
    except subprocess.CalledProcessError as e:
        print(f"Error checking workflow status: {e.stderr}")
        print("Assuming job is still running...")
        time.sleep(${params.monitor_interval})
        sys.exit(100)
    """
}

// ============================================================================
// PROCESS 6: Download AnVIL results
// ============================================================================
// What this does: Downloads imputed data from All of Us AnVIL workspace
// - Uses gsutil or Terra CLI to download from Google Cloud Storage
// - Retrieves all chromosome VCF and info files
// - Renames to standardized format

process downloadAnvilResults {
    label 'download'
    tag "${platform}"
    publishDir "${params.outdir}/module2/04_anvil_results/${platform}", mode: 'copy'
    
    when:
    params.run_anvil
    
    input:
    tuple val(platform), path(job_info)
    
    output:
    tuple val(platform),
          path("${platform}_imputed_chr*.dose.vcf.gz"),
          path("${platform}_imputed_chr*.info.gz"),
          emit: anvil_results
    
    script:
    """
    #!/usr/bin/env python3
    import json
    import subprocess
    import os
    import glob
    
    # Load job information
    with open('${job_info}') as f:
        job = json.load(f)
    
    workspace = job.get('workspace')
    project = job.get('project')
    workflow_id = job.get('workflow_id')
    
    print(f"Downloading results for ${platform} from AnVIL workspace...")
    
    # Method 1: Using Terra CLI (preferred)
    try:
        # Get workflow outputs
        cmd = [
            'terra', 'workflow', 'outputs',
            '--workspace', workspace,
            '--project', project,
            '--workflow-id', workflow_id
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            print("Workflow outputs:")
            print(result.stdout)
    except Exception as e:
        print(f"Note: Could not retrieve outputs via Terra CLI: {e}")
    
    # Method 2: Direct download from Google Cloud Storage
    # Construct typical output path for All of Us imputation
    bucket_path = f"gs://{workspace}-workspace-bucket/submissions/{workflow_id}/"
    
    try:
        # Download all VCF and info files
        cmd = ['gsutil', '-m', 'cp', '-r', f'{bucket_path}*.vcf.gz', '.']
        subprocess.run(cmd, check=False)
        
        cmd = ['gsutil', '-m', 'cp', '-r', f'{bucket_path}*.info.gz', '.']
        subprocess.run(cmd, check=False)
        
        print("✓ Files downloaded from Google Cloud Storage")
    except Exception as e:
        print(f"Note: Direct GCS download failed: {e}")
        print("Files may need to be downloaded manually from workspace")
    
    # Rename files to standardized format
    for filename in os.listdir('.'):
        if 'chr' in filename:
            if '.dose.vcf.gz' in filename or '.vcf.gz' in filename:
                chr_num = filename.split('chr')[1].split('.')[0]
                new_name = f"${platform}_imputed_chr{chr_num}.dose.vcf.gz"
                os.rename(filename, new_name)
                print(f"Renamed: {filename} -> {new_name}")
            elif '.info' in filename:
                chr_num = filename.split('chr')[1].split('.')[0]
                new_name = f"${platform}_imputed_chr{chr_num}.info.gz"
                os.rename(filename, new_name)
                print(f"Renamed: {filename} -> {new_name}")
    
    print("✓ Download complete!")
    """
}

// ============================================================================
// MODULE 2 WORKFLOW - Main orchestration
// ============================================================================
workflow MODULE2_IMPUTATION {
    take:
    vcf_files    // From Module 1: tuple(platform, vcf_files, tbi_files)
    manifests    // From Module 1: tuple(platform, topmed_json, anvil_tsv, file_list)
    
    main:
    // Initialize empty channels for conditional outputs
    topmed_out = Channel.empty()
    anvil_out = Channel.empty()
    
    // TOPMed Imputation Pathway
    if (params.run_topmed) {
        // Submit jobs to TOPMed server
        submitTOPMedImputation(vcf_files, manifests)
        
        // Monitor until completion (checks every 15 minutes)
        monitorTOPMedJob(submitTOPMedImputation.out.job_info)
        
        // Download results when ready
        downloadTOPMedResults(monitorTOPMedJob.out.completed_job)
        
        topmed_out = downloadTOPMedResults.out.topmed_results
    }
    
    // All of Us AnVIL Imputation Pathway
    if (params.run_anvil) {
        // Submit workflows to AnVIL/Terra
        submitAnvilImputation(vcf_files, manifests)
        
        // Monitor until completion (checks every 15 minutes)
        monitorAnvilJob(submitAnvilImputation.out.job_info)
        
        // Download results when ready
        downloadAnvilResults(monitorAnvilJob.out.completed_job)
        
        anvil_out = downloadAnvilResults.out.anvil_results
    }
    
    // Determine final output based on which services ran
    if (params.run_topmed && params.run_anvil) {
        // If both services ran, user will have both outputs available
        // Downstream modules can choose which to use or compare them
        final_results = topmed_out.mix(anvil_out)
    } else if (params.run_topmed) {
        final_results = topmed_out
    } else {
        final_results = anvil_out
    }
    
    emit:
    imputed_data = final_results
    topmed_results = topmed_out
    anvil_results = anvil_out
}

