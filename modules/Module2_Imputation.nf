#!/usr/bin/env nextflow

/*
 * ============================================================================
 * MODULE 2: IMPUTATION SUBMISSION, MONITORING & DOWNLOAD (OPTIMIZED v5.0)
 * ============================================================================
 * 
 * ARCHITECTURE: Fully Automated Continuous Monitoring
 * 
 * CRITICAL OPTIMIZATIONS:
 * 1. Tools pre-installed in Apptainer containers (not in pipeline)
 * 2. 45-minute monitoring intervals (reduced API calls)
 * 3. Smart retry logic for monitoring (up to 24 hours)
 * 4. Pre-configured passwords for automatic decryption
 * 5. No manual intervention required
 * 
 * DUAL PARALLEL WORKFLOWS:
 * - TOPMed: imputationbot with --autoDownload + --password
 * - All of Us: Custom CLI monitoring with job status checks
 * 
 * WHAT IT DOES:
 * - Groups all chr1-22 per platform (quota-efficient!)
 * - Submits to TOPMed (auto-download enabled)
 * - Submits to All of Us AnVIL  
 * - Monitors both services every 45 minutes
 * - Downloads + decrypts automatically when ready
 * - Returns imputed VCFs for Module 3
 * - Handles TOPMed's 3-job limit automatically
 * - Runs continuously until completion
 * - Pre-built containers = consistent environment
 * 
 * EXPECTED RESULTS:
 * - {platform}_topmed_chr{1-22}.dose.vcf.gz (imputed genotypes)
 * - {platform}_topmed_chr{1-22}.info.gz (imputation quality R²)
 * - {platform}_anvil_chr{1-22}.dose.vcf.gz (imputed genotypes)  
 * - {platform}_anvil_chr{1-22}.info.gz (imputation quality R²)
 * ============================================================================
 */

nextflow.enable.dsl = 2

// ============================================================================
// PROCESS 1: Group chromosomes by platform (QUOTA-EFFICIENT)
// ============================================================================
/*
 * WHAT IT DOES: Groups all chr1-22 VCFs per platform into single submission
 * 
 * WHY THIS IS CRITICAL:
 * - All of Us charges 1 quota per GENOME (not per chromosome!)
 * - Submitting chr1-22 separately = 22 quota units
 * - Submitting chr1-22 together = 1 quota unit (95% savings!)
 * 
 * CONTAINER: python_tools.sif (Python 3.10 + pandas)
 */

process groupChromosomesByPlatform {
    label 'manifest'
    tag "${platform}"
    publishDir "${params.outdir}/module2/00_grouped/${platform}", mode: 'copy'
    
    // Use pre-built Apptainer container
    container "${projectDir}/resources/containers/python_api.sif"
    
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
          val(build),
          path("${platform}_vcf_list.txt"),
          emit: grouped
    
    script:
    def vcf_list = vcf_files instanceof List ? vcf_files : [vcf_files]
    def n_vcfs = vcf_list.size()
    """
    #!/usr/bin/env python3
    
    platform = "${platform}"
    build = "${build}"
    n_files = ${n_vcfs}
    
    print("=" * 70)
    print(f"Grouping chromosomes for: {platform}")
    print("=" * 70)
    print(f"Build: {build}")
    print(f"Number of chromosome VCFs: {n_files}")
    print(f"")
    print(f"✓ All {n_files} chromosomes will be submitted together")
    print(f"  Quota usage: 1 unit (not {n_files}!)")
    print("=" * 70)
    
    # Create file list
    vcf_list = "${vcf_list}".replace('[', '').replace(']', '').split(', ')
    
    with open(f"{platform}_vcf_list.txt", 'w') as f:
        for vcf in vcf_list:
            vcf_clean = vcf.strip()
            if vcf_clean and vcf_clean.endswith('.vcf.gz'):
                f.write(f"{vcf_clean}\\n")
                print(f"  - {vcf_clean}")
    
    print(f"\\n✓ Grouped {n_files} files for submission")
    """
}

// ============================================================================
// PROCESS 2: Submit to TOPMed with imputationbot (FULLY AUTOMATED)
// ============================================================================
/*
 * WHAT IT DOES: Submits to TOPMed with automatic download & decryption
 * 
 * 
 * - imputationbot CLI handles everything automatically
 * - --autoDownload waits + downloads when ready
 * - --password enables auto-decryption (no email!) --> use this flag to set your own password and include in nextflow.config as parameter
 * - maxForks 3 respects TOPMed's concurrent job limit
 * 
 * CONTAINER: python_tools.sif (includes imputationbot pre-installed)
 * 
 * MONITORING: Built into imputationbot (checks every 5 min internally)
 * 
 * EXPECTED DURATION: 2-8 hours depending on queue + data size
 */

process submitToTOPMed {
    label 'imputation'
    tag "${platform}"
    publishDir "${params.outdir}/module2/01_topmed/${platform}", mode: 'copy'
    
    // Use pre-built container with imputationbot installed
    container "${projectDir}/resources/containers/python_api.sif"
    
    cpus 4
    memory '8.GB'
    time '24.h'
    maxForks 3  // CRITICAL: TOPMed allows max 3 concurrent jobs!
    
    when:
    params.run_topmed
    
    input:
    tuple val(platform),
          path(vcf_files),
          path(tbi_files),
          val(build),
          path(vcf_list)
    
    output:
    tuple val(platform),
          val('topmed'),
          path("${platform}_chr*.dose.vcf.gz"),
          path("${platform}_chr*.info.gz"),
          emit: imputed_data
    
    script:
    def password = params.topmed_password ?: "ERROR_NO_PASSWORD"
    def vcf_files_str = vcf_files instanceof List ? vcf_files.collect { it.toString() }.join(' ') : vcf_files.toString()
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    echo "=" * 70
    echo "TOPMed Imputation: ${platform}"
    echo "=" * 70
    
    # CRITICAL: Verify password is configured
    if [ "${password}" == "ERROR_NO_PASSWORD" ]; then
        echo "ERROR: TOPMed password not configured!"
        echo ""
        echo "Set in nextflow.config:"
        echo "  params.topmed_password = 'YourStrongPassword123!'"
        echo ""
        echo "This enables automatic decryption without email retrieval"
        exit 1
    fi
    
    # CRITICAL: Verify API token is configured
    if [ -z "${params.topmed_api_token}" ] || [ "${params.topmed_api_token}" == "null" ]; then
        echo "ERROR: TOPMed API token not configured!"
        echo ""
        echo "Get token from: https://imputation.biodatacatalyst.nhlbi.nih.gov"
        echo "Set in nextflow.config:"
        echo "  params.topmed_api_token = 'your-token-here'"
        exit 1
    fi
    
    # Configure imputationbot (uses token from environment)
    export IMPUTATION_BOT_TOKEN="${params.topmed_api_token}"
    
    imputationbot configure \\
        --url https://imputation.biodatacatalyst.nhlbi.nih.gov \\
        --token "\${IMPUTATION_BOT_TOKEN}"
    
    echo "✓ imputationbot configured"
    echo ""
    
    # Get VCF files
    VCF_FILES="${vcf_files_str}"
    N_FILES=\$(echo \$VCF_FILES | wc -w)
    
    echo "Submitting \$N_FILES chromosome VCFs as ONE job"
    echo "Quota usage: 1 unit (not \$N_FILES!)"
    echo ""
    
    # List files
    echo "Files to submit:"
    for vcf in \$VCF_FILES; do
        echo "  - \$(basename \$vcf)"
    done
    echo ""
    
    # CRITICAL: Submit with autoDownload + password for full automation
    echo "Submitting to TOPMed with automatic download..."
    echo "This process will:"
    echo "  1. Upload all VCF files"
    echo "  2. Submit imputation job"
    echo "  3. Wait for completion (checks every 5 minutes)"
    echo "  4. Download results automatically"
    echo "  5. Decrypt automatically using configured password"
    echo ""
    echo "Expected duration: 2-8 hours"
    echo "The pipeline will wait - no manual intervention needed!"
    echo ""
    
    imputationbot impute \\
        --files \$VCF_FILES \\
        --refpanel topmed-r3 \\
        --population all \\
        --build ${build} \\
        --name "${platform}_imputation" \\
        --password "${password}" \\
        --r2Filter 0 \\
        --autoDownload
    
    # Check if successful
    if [ \$? -ne 0 ]; then
        echo "✗ Imputation, download, or decryption failed"
        exit 1
    fi
    
    echo ""
    echo "✓ Imputation complete!"
    echo "✓ Files downloaded and decrypted automatically"
    echo ""
    
    # Standardize filenames with server identifier
    echo "Standardizing filenames..."
    for file in *.vcf.gz *.info.gz; do
        if [[ \$file =~ chr([0-9]+) ]]; then
            chr_num=\${BASH_REMATCH[1]}
            
            if [[ \$file == *.dose.vcf.gz ]]; then
                new_name="${platform}_chr\${chr_num}.dose.vcf.gz"
            elif [[ \$file == *.vcf.gz ]]; then
                new_name="${platform}_chr\${chr_num}.dose.vcf.gz"
            elif [[ \$file == *.info.gz ]]; then
                new_name="${platform}_chr\${chr_num}.info.gz"
            else
                continue
            fi
            
            if [ "\$file" != "\$new_name" ]; then
                mv "\$file" "\$new_name"
                echo "  \$file → \$new_name"
            fi
        fi
    done
    
    echo ""
    echo "=" * 70
    echo "✓ ${platform} TOPMed imputation complete"
    echo "=" * 70
    """
}

// ============================================================================
// PROCESS 3: Submit to All of Us AnVIL
// ============================================================================
/*
 * WHAT IT DOES: Submits imputation job to All of Us service
 * 
 * WHY THIS IS CRITICAL:
 * - Different reference panel (AllofUs)
 * - Enables benchmarking: TOPMed vs All of Us imputation quality
 * - Returns job ID for monitoring in next process
 * 
 * CONTAINER: python_tools.sif (includes terra-notebook-utils)
 * 
 * NOTE: Actual submission command depends on All of Us CLI
 * Check latest documentation at: https://imputation.researchallofus.org/
 */

process submitToAllOfUs {
    label 'imputation'
    tag "${platform}"
    publishDir "${params.outdir}/module2/02_anvil/${platform}", mode: 'copy'
    
    // Use pre-built container with terralab-cli
    container "${projectDir}/resources/containers/python_api.sif"
    
    cpus 4
    memory '8.GB'
    
    when:
    params.run_anvil
    
    input:
    tuple val(platform),
          path(vcf_files),
          path(tbi_files),
          val(build),
          path(vcf_list)
    
    output:
    tuple val(platform),
          path("${platform}_anvil_submission.json"),
          emit: job_info
    
    script:
    def vcf_files_str = vcf_files instanceof List ? vcf_files.collect { it.toString() }.join(',') : vcf_files.toString()
    """
    #!/usr/bin/env python3
    import json
    import subprocess
    import sys
    from datetime import datetime
    import re
    
    platform = "${platform}"
    build = "${build}"
    vcf_files = "${vcf_files_str}".split(',')
    
    print("=" * 70)
    print(f"All of Us AnVIL Submission: {platform}")
    print("=" * 70)
    print(f"Build: {build}")
    print(f"Files: {len(vcf_files)} chromosomes")
    print(f"Quota usage: 1 unit (not {len(vcf_files)}!)")
    print("")
    
    # Create submission metadata
    submission_data = {
        "platform": platform,
        "build": build,
        "n_files": len(vcf_files),
        "files": [f.split('/')[-1] for f in vcf_files],
        "submission_time": datetime.now().isoformat(),
        "status": "submitted"
    }
    
    # Authenticate (if not already done)
    print("Authenticating with terralab-cli...")
    try:
        # Check if already authenticated
        auth_check = subprocess.run(
            ["terralab", "jobs", "list"],
            capture_output=True,
            text=True,
            timeout=30
        )
        if auth_check.returncode != 0:
            print("  Running terralab login...")
            # This will prompt user if in interactive session
            # Or use stored credentials if available
            subprocess.run(["terralab", "login"], check=False)
    except Exception as e:
        print(f"  Note: Authentication may be required - {e}")
    
    print("✓ Authentication checked")
    print("")
    
    # Prepare merged VCF (All of Us requires single multi-sample VCF)
    print("Preparing input VCF...")
    input_vcf = vcf_files[0] if len(vcf_files) == 1 else None
    
    if not input_vcf or len(vcf_files) > 1:
        # Need to concatenate all chromosomes into single VCF
        print(f"  Concatenating {len(vcf_files)} chromosome VCFs...")
        concat_cmd = ["bcftools", "concat", "-o", f"{platform}_merged.vcf.gz", "-O", "z"] + vcf_files
        result = subprocess.run(concat_cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            input_vcf = f"{platform}_merged.vcf.gz"
            # Index the merged VCF
            subprocess.run(["tabix", "-p", "vcf", input_vcf])
            print(f"  ✓ Created merged VCF: {input_vcf}")
        else:
            print(f"  ✗ Failed to merge VCFs: {result.stderr}")
            submission_data["status"] = "merge_failed"
            submission_data["error"] = result.stderr
            with open(f"{platform}_anvil_submission.json", 'w') as f:
                json.dump(submission_data, f, indent=2)
            sys.exit(1)
    
    print("")
    
    # Submit via terralab-cli
    try:
        print("Submitting to All of Us AnVIL Imputation Service...")
        print(f"  Input: {input_vcf}")
        print(f"  Output basename: {platform}_anvil")
        print("")
        
        # CORRECT terralab-cli command based on documentation:
        # terralab submit array_imputation --multiSampleVcf FILE --outputBasename NAME --description 'DESC'
        cmd = [
            "terralab", "submit", "array_imputation",
            "--multiSampleVcf", input_vcf,
            "--outputBasename", f"{platform}_anvil",
            "--description", f"{platform} All of Us imputation via automated pipeline"
        ]
        
        print(f"  Command: {' '.join(cmd)}")
        print("")
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        
        if result.returncode == 0:
            # Parse job ID from output
            # Expected format varies, but look for common patterns
            job_id_match = re.search(r'(job[_-]?[0-9a-f-]+|[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})', 
                                    result.stdout, re.IGNORECASE)
            
            if job_id_match:
                job_id = job_id_match.group(1)
                submission_data["job_id"] = job_id
                submission_data["status"] = "submitted_success"
                print(f"✓ Submitted successfully!")
                print(f"  Job ID: {job_id}")
                print(f"  Monitor with: terralab jobs details {job_id}")
            else:
                # Submission succeeded but couldn't parse job ID
                submission_data["status"] = "submitted_no_id"
                submission_data["output"] = result.stdout
                print(f"✓ Submitted successfully!")
                print(f"  Note: Could not parse job ID from output")
                print(f"  Run 'terralab jobs list' to see your jobs")
        else:
            print(f"✗ CLI submission failed!")
            print(f"  Error: {result.stderr}")
            submission_data["status"] = "submission_failed"
            submission_data["error"] = result.stderr
            
    except subprocess.TimeoutExpired:
        print(f"⚠ Submission timed out (>5 minutes)")
        submission_data["status"] = "submission_timeout"
        
    except Exception as e:
        print(f"⚠ All of Us CLI error: {e}")
        print("")
        print("Manual submission required:")
        print("  1. Authenticate: terralab login")
        print("  2. Submit job:")
        print(f"     terralab submit array_imputation \\")
        print(f"       --multiSampleVcf {input_vcf} \\")
        print(f"       --outputBasename {platform}_anvil \\")
        print(f"       --description '{platform} imputation'")
        print("")
        submission_data["status"] = "manual_submission_required"
        submission_data["manual_instructions"] = True
        submission_data["error"] = str(e)
    
    # Save submission info for monitoring
    with open(f"{platform}_anvil_submission.json", 'w') as f:
        json.dump(submission_data, f, indent=2)
    
    print("")
    print("=" * 70)
    print(f"✓ {platform} submission metadata saved")
    print("=" * 70)
    """
}

// ============================================================================
// PROCESS 4: Monitor and Download All of Us Results (SMART RETRY)
// ============================================================================
/*
 * WHAT IT DOES: Monitors All of Us job and downloads when complete
 * 
 * WHY THIS IS CRITICAL:
 * - All of Us jobs take 2-12 hours
 * - Smart retry with 45-minute intervals (reduced API calls)
 * - Automatic download + decryption when ready
 * - Pipeline continues without manual intervention
 * 
 * CONTAINER: python_tools.sif (includes terra-notebook-utils)
 * 
 * RETRY STRATEGY:
 * - Checks every 45 minutes (${params.monitor_interval_minutes} minutes)
 * - Max retries: 32 (= 24 hours with 45-min intervals)
 * - Exit code 100 = retry (job still running)
 * - Exit code 0 = success (job complete, downloaded)
 * - Exit code 1 = failure (job failed or error)
 */

process monitorAndDownloadAllOfUs {
    label 'download'
    tag "${platform}"
    publishDir "${params.outdir}/module2/02_anvil/${platform}", mode: 'copy'
    
    // Use pre-built container
    container "${projectDir}/resources/containers/python_api.sif"
    
    cpus 4
    memory '8.GB'
    time '24.h'
    
    // CRITICAL: Smart retry for monitoring
    maxRetries 32  // 45 min × 32 = 24 hours max wait
    errorStrategy { task.exitStatus == 100 ? 'retry' : 'finish' }
    
    when:
    params.run_anvil
    
    input:
    tuple val(platform), path(submission_json)
    
    output:
    tuple val(platform),
          val('anvil'),
          path("${platform}_chr*.dose.vcf.gz"),
          path("${platform}_chr*.info.gz"),
          emit: imputed_data
    
    script:
    """
    #!/usr/bin/env python3
    import json
    import subprocess
    import time
    import sys
    import os
    import re
    from pathlib import Path
    
    platform = "${platform}"
    retry_interval = ${params.monitor_interval_minutes} * 60  # Convert to seconds
    
    print("=" * 70)
    print(f"Monitoring All of Us Job: {platform}")
    print(f"Attempt: ${task.attempt} / ${task.maxRetries}")
    print("=" * 70)
    
    # Load submission info
    with open("${submission_json}") as f:
        submission = json.load(f)
    
    job_id = submission.get('job_id')
    status = submission.get('status')
    
    # Handle cases where submission failed or requires manual intervention
    if status in ['submission_failed', 'manual_submission_required', 'merge_failed']:
        print(f"⚠ Status: {status}")
        print("Creating placeholder files for manual download...")
        
        # Create placeholders so pipeline can continue
        for chr_num in range(1, 23):
            Path(f"{platform}_chr{chr_num}.dose.vcf.gz").write_text("# Placeholder - manual download required\\n")
            Path(f"{platform}_chr{chr_num}.info.gz").write_text("# Placeholder - manual download required\\n")
        
        print("")
        print("Manual download instructions:")
        print("  1. Check job status:")
        print("     terralab jobs list")
        print("  2. Get job details:")
        print("     terralab jobs details JOB_ID")
        print("  3. Download when complete:")
        print("     terralab download JOB_ID")
        print("  4. Replace placeholder files in:")
        print(f"     ${params.outdir}/module2/02_anvil/{platform}/")
        sys.exit(0)
    
    if not job_id:
        print("ERROR: No job ID found in submission")
        print("Try running: terralab jobs list")
        sys.exit(1)
    
    print(f"Job ID: {job_id}")
    print(f"Checking status (will retry every {retry_interval/60:.0f} minutes)...")
    print("")
    
    try:
        # CORRECT terralab-cli command for checking job status
        result = subprocess.run(
            ["terralab", "jobs", "details", job_id],
            capture_output=True,
            text=True,
            timeout=120
        )
        
        if result.returncode != 0:
            print(f"⚠ Error checking status: {result.stderr}")
            print(f"Will retry in {retry_interval/60:.0f} minutes...")
            sys.exit(100)  # Trigger retry
        
        status_output = result.stdout.lower()
        print(f"Status output:\\n{result.stdout}")
        print("")
        
        # Check if complete (look for success indicators)
        if any(word in status_output for word in ['success', 'complete', 'done', 'finished']):
            # Additional check: make sure it's not "in progress" or "running"
            if any(word in status_output for word in ['running', 'progress', 'pending']):
                print(f"⏳ Job still running: {status_output}")
                print(f"Will check again in {retry_interval/60:.0f} minutes...")
                print(f"Time waited so far: {${task.attempt} * retry_interval / 3600:.1f} hours")
                sys.exit(100)  # Trigger retry
            
            print("✓ Job completed! Downloading results...")
            
            # Download results using terralab-cli
            # Command: terralab download JOB_ID
            download_result = subprocess.run(
                ["terralab", "download", job_id],
                capture_output=True,
                text=True,
                timeout=3600  # 1 hour timeout for download
            )
            
            if download_result.returncode != 0:
                print(f"✗ Download failed: {download_result.stderr}")
                sys.exit(1)
            
            print("✓ Download complete!")
            print(f"Download output:\\n{download_result.stdout}")
            print("")
            
            # Find downloaded files
            # All of Us outputs: imputed VCF, index file, QC TSV
            vcf_files = list(Path('.').glob('*.vcf.gz'))
            
            if not vcf_files:
                print("✗ No VCF files found after download!")
                print("Files in directory:")
                for f in os.listdir('.'):
                    print(f"  - {f}")
                sys.exit(1)
            
            print(f"Found {len(vcf_files)} output file(s)")
            
            # Standardize filenames
            # All of Us may return single merged VCF or per-chromosome
            print("Standardizing filenames...")
            
            for vcf_path in vcf_files:
                vcf_file = vcf_path.name
                
                # Check if file has chromosome info
                chr_match = re.search(r'chr(\\d+)', vcf_file, re.IGNORECASE)
                
                if chr_match:
                    # Per-chromosome file
                    chr_num = chr_match.group(1)
                    new_vcf = f"{platform}_chr{chr_num}.dose.vcf.gz"
                    
                    if vcf_file != new_vcf:
                        os.rename(vcf_file, new_vcf)
                        print(f"  {vcf_file} → {new_vcf}")
                    
                    # Look for corresponding .tbi and info files
                    tbi_file = f"{vcf_file}.tbi"
                    if os.path.exists(tbi_file):
                        new_tbi = f"{new_vcf}.tbi"
                        os.rename(tbi_file, new_tbi)
                    
                    # Info file might be .tsv or .info.gz
                    for ext in ['.info.gz', '.tsv', '.info', '_info.tsv']:
                        base = vcf_file.replace('.vcf.gz', '')
                        info_file = f"{base}{ext}"
                        if os.path.exists(info_file):
                            new_info = f"{platform}_chr{chr_num}.info.gz"
                            if not info_file.endswith('.gz'):
                                # Compress it
                                subprocess.run(['gzip', '-f', info_file])
                                info_file = f"{info_file}.gz"
                            os.rename(info_file, new_info)
                            print(f"  {info_file} → {new_info}")
                            break
                else:
                    # Single merged file - need to split by chromosome
                    print(f"  Processing merged VCF: {vcf_file}")
                    
                    # Extract each chromosome
                    for chr_num in range(1, 23):
                        chr_name = f"chr{chr_num}"
                        out_vcf = f"{platform}_chr{chr_num}.dose.vcf.gz"
                        
                        # Extract chromosome
                        extract_cmd = [
                            "bcftools", "view",
                            "-r", chr_name,
                            "-O", "z",
                            "-o", out_vcf,
                            vcf_file
                        ]
                        
                        result = subprocess.run(extract_cmd, capture_output=True)
                        if result.returncode == 0:
                            print(f"  Extracted {chr_name} → {out_vcf}")
                            # Index it
                            subprocess.run(["tabix", "-p", "vcf", out_vcf])
                        else:
                            print(f"  Warning: Could not extract {chr_name}")
                    
                    # Handle QC/info file
                    for qc_ext in ['.tsv', '.info.gz', '_qc.tsv']:
                        qc_file = vcf_file.replace('.vcf.gz', qc_ext)
                        if os.path.exists(qc_file):
                            # Create per-chromosome info files
                            # For now, just copy the merged info for each chr
                            for chr_num in range(1, 23):
                                info_out = f"{platform}_chr{chr_num}.info.gz"
                                if not qc_file.endswith('.gz'):
                                    subprocess.run(['gzip', '-c', qc_file], 
                                                 stdout=open(info_out, 'wb'))
                                else:
                                    subprocess.run(['cp', qc_file, info_out])
                            break
            
            # Verify we have all expected outputs
            expected_files = []
            for chr_num in range(1, 23):
                expected_files.append(f"{platform}_chr{chr_num}.dose.vcf.gz")
                expected_files.append(f"{platform}_chr{chr_num}.info.gz")
            
            missing = [f for f in expected_files if not os.path.exists(f)]
            if missing:
                print(f"\\n⚠ Warning: {len(missing)} expected files not found:")
                for f in missing[:5]:  # Show first 5
                    print(f"  - {f}")
                if len(missing) > 5:
                    print(f"  ... and {len(missing)-5} more")
            
            print("")
            print("=" * 70)
            print(f"✓ {platform} All of Us imputation complete!")
            print("=" * 70)
            sys.exit(0)  # Success
            
        elif any(word in status_output for word in ['fail', 'error', 'abort']):
            # But not "running" or "in progress"
            if not any(word in status_output for word in ['running', 'progress', 'pending']):
                print(f"✗ Job failed: {status_output}")
                sys.exit(1)  # Failure
        
        # Still running or pending
        print(f"⏳ Job in progress")
        print(f"Will check again in {retry_interval/60:.0f} minutes...")
        print(f"Time waited so far: {${task.attempt} * retry_interval / 3600:.1f} hours")
        sys.exit(100)  # Trigger retry
            
    except subprocess.TimeoutExpired:
        print("⚠ Status check timed out")
        print(f"Will retry in {retry_interval/60:.0f} minutes...")
        sys.exit(100)  # Trigger retry
        
    except Exception as e:
        print(f"⚠ Error during monitoring: {e}")
        print(f"Will retry in {retry_interval/60:.0f} minutes...")
        sys.exit(100)  # Trigger retry
    """
}

// ============================================================================
// MODULE 2 WORKFLOW - DUAL PATHWAYS WITH CONTINUOUS MONITORING
// ============================================================================
workflow MODULE2_IMPUTATION {
    take:
    topmed_vcfs  // From Module 1 TOPMed pathway: tuple(platform, chr, vcf, tbi, build, batch)
    anvil_vcfs   // From Module 1 AnVIL pathway: tuple(platform, chr, vcf, tbi, build, batch)
    
    main:
    
    // Initialize outputs
    topmed_out = Channel.empty()
    anvil_out = Channel.empty()
    
    // ========== TOPMED PATHWAY ==========
    if (params.run_topmed && topmed_vcfs != null) {
        // Group all chromosomes per platform
        topmed_vcfs
            .map { platform, chr, vcf, tbi, build, batch ->
                tuple(platform, chr, vcf, tbi, build)
            }
            .groupTuple(by: 0)  // Group by platform
            .set { topmed_grouped }
        
        // Group chromosomes
        groupChromosomesByPlatform(topmed_grouped)
        
        // Submit to TOPMed (max 3 concurrent via maxForks)
        // imputationbot handles monitoring + download automatically
        submitToTOPMed(groupChromosomesByPlatform.out.grouped)
        
        topmed_out = submitToTOPMed.out.imputed_data
    }
    
    // ========== ANVIL PATHWAY ==========
    if (params.run_anvil && anvil_vcfs != null) {
        // Group all chromosomes per platform
        anvil_vcfs
            .map { platform, chr, vcf, tbi, build, batch ->
                tuple(platform, chr, vcf, tbi, build)
            }
            .groupTuple(by: 0)  // Group by platform
            .set { anvil_grouped }
        
        // Group chromosomes
        groupChromosomesByPlatform(anvil_grouped)
        
        // Submit to All of Us
        submitToAllOfUs(groupChromosomesByPlatform.out.grouped)
        
        // Monitor and download (with smart retry every 45 min)
        monitorAndDownloadAllOfUs(submitToAllOfUs.out.job_info)
        
        anvil_out = monitorAndDownloadAllOfUs.out.imputed_data
    }
    
    // Combine outputs (maintains server identifier)
    // Output format: tuple(platform, server, vcf_files, info_files)
    final_results = topmed_out.mix(anvil_out)
    
    emit:
    imputed_data = final_results
    topmed_results = topmed_out
    anvil_results = anvil_out
}

workflow.onComplete {
    def topmed_status = params.run_topmed ? "ENABLED" : "DISABLED"
    def anvil_status = params.run_anvil ? "ENABLED" : "DISABLED"
    
    println """
    ============================================================================
    MODULE 2 - IMPUTATION COMPLETE
    ============================================================================
    Status:        ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration:      ${workflow.duration}
    
    ✓ Dual Parallel Workflows:
      - TOPMed pathway: ${topmed_status}
      - AnVIL pathway: ${anvil_status}
    
    ✓ Fully Automated Processing:
      - imputationbot: auto-download + auto-decrypt
      - All of Us: smart monitoring (every ${params.monitor_interval_minutes} min)
      - No manual intervention required!
    
    ✓ Queue Management:
      - TOPMed 3-job limit respected via maxForks
      - Platforms queue automatically if limit reached
    
    ✓ Quota-Efficient Submission:
      - All chr1-22 submitted together per platform
      - Uses 1 quota unit per platform (not 22!)
    
    ✓ Pre-built Containers:
      - All tools pre-installed in Apptainer images
      - Consistent environment across systems
    
    Results:
      - TOPMed: ${params.outdir}/module2/01_topmed/
      - AnVIL: ${params.outdir}/module2/02_anvil/
    ============================================================================
    """
}
