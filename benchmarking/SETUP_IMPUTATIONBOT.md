# Setting Up Imputationbot for Automated Benchmarking

This guide covers setting up automated submission to imputation servers.

## Overview

The benchmark matrix runs 24 combinations (6 approaches × 4 servers):

| Server | Panel | Description |
|--------|-------|-------------|
| **TOPMed** | TOPMed r2 | Most diverse, best for admixed |
| **All of Us** | TOPMed-based | 50%+ non-EUR, diverse |
| **Michigan HRC** | HRC r1.1 | Largest, EUR-focused |
| **Michigan 1KG** | 1000G Phase 3 | Global diversity |

## 1. Install Imputationbot

```bash
# Option 1: Quick install
curl -sL imputationbot.now.sh | bash

# Option 2: Via pip
pip install imputationbot

# Verify installation
imputationbot --version
```

## 2. Get API Tokens

### TOPMed Imputation Server

1. Go to https://imputation.biodatacatalyst.nhlbi.nih.gov
2. Create an account or log in
3. Go to Profile → API Access
4. Generate new token
5. Save token securely

### Michigan Imputation Server

1. Go to https://imputationserver.sph.umich.edu
2. Create an account or log in
3. Go to Profile → API Access
4. Generate new token
5. Save token securely

### All of Us (Terra/AnVIL)

```bash
# Install terralab CLI
pip install terralab-cli

# Authenticate (opens browser)
terralab login

# Verify authentication
terralab auth status
```

## 3. Configure Imputationbot

```bash
# Add TOPMed instance
imputationbot add-instance topmed \
    --url https://imputation.biodatacatalyst.nhlbi.nih.gov \
    --token YOUR_TOPMED_TOKEN

# Add Michigan instance
imputationbot add-instance michigan \
    --url https://imputationserver.sph.umich.edu \
    --token YOUR_MICHIGAN_TOKEN

# List configured instances
imputationbot instances
```

## 4. Test Submission

```bash
# Test TOPMed submission (small file)
imputationbot impute \
    --instance topmed \
    --files test_chr22.vcf.gz \
    --refpanel topmed-r2 \
    --population mixed \
    --name test_job

# Test Michigan HRC
imputationbot impute \
    --instance michigan \
    --files test_chr22.vcf.gz \
    --refpanel hrc-r1.1 \
    --population mixed \
    --name test_job_hrc

# Test Michigan 1KG
imputationbot impute \
    --instance michigan \
    --files test_chr22.vcf.gz \
    --refpanel 1000g-phase3-v5 \
    --population mixed \
    --name test_job_1kg
```

## 5. Run Full Benchmark

### Option A: Environment Variables

```bash
export TOPMED_TOKEN="your_topmed_token"
export TOPMED_PASSWORD="your_password"
export MICHIGAN_TOKEN="your_michigan_token"

./run_full_benchmark_matrix.sh \
    --input test_data \
    --output benchmark_results
```

### Option B: Config File

Create `benchmark_config.yaml`:

```yaml
credentials:
  topmed:
    token: "your_topmed_token"
    password: "your_password"
  michigan:
    token: "your_michigan_token"
  allofus:
    # Uses terralab auth
    authenticated: true

settings:
  threads: 8
  poll_interval: 300  # 5 minutes
```

```bash
./run_full_benchmark_matrix.sh \
    --input test_data \
    --config benchmark_config.yaml
```

## 6. Monitor Jobs

```bash
# List all jobs
imputationbot jobs

# Check specific job
imputationbot jobs --id job-abc123

# Download results when complete
imputationbot download --id job-abc123 --output results/
```

## 7. Reference Panels Available

### TOPMed Server
- `topmed-r2` - TOPMed Freeze 8 (97,256 samples, diverse)

### Michigan Server
- `hrc-r1.1` - Haplotype Reference Consortium (64,976 samples, EUR-focused)
- `1000g-phase3-v5` - 1000 Genomes Phase 3 (2,504 samples, global)
- `caapa` - CAAPA (African ancestry)
- `genome-asia-pilot` - GenomeAsia (Asian ancestry)

### All of Us
- TOPMed-based panel with enhanced diversity

## 8. Troubleshooting

### Authentication Issues

```bash
# Refresh TOPMed token
imputationbot update-instance topmed --token NEW_TOKEN

# Re-authenticate terralab
terralab logout
terralab login
```

### Job Failures

```bash
# Check job logs
imputationbot logs --id job-abc123

# Common issues:
# - VCF format errors: Ensure proper chromosome naming
# - Sample mismatch: Check sample IDs are consistent
# - Reference mismatch: Verify genome build (hg19/hg38)
```

### Network Issues

```bash
# Retry failed download
imputationbot download --id job-abc123 --output results/ --retry 3
```

## 9. Best Practices

1. **Start small**: Test with chr22 only first
2. **Monitor costs**: Some servers have usage limits
3. **Save tokens securely**: Don't commit to git
4. **Check job status**: Don't submit duplicates
5. **Validate inputs**: Run bcftools check before submitting

## 10. Expected Timings

Approximate wall-clock times for 2,504 samples (1KG):

| Server | Per-chromosome | Full genome |
|--------|---------------|-------------|
| TOPMed | 1-2 hours | 8-16 hours |
| Michigan HRC | 30-60 min | 4-8 hours |
| Michigan 1KG | 20-40 min | 3-6 hours |
| All of Us | 1-2 hours | 8-16 hours |

Note: Queue times vary by server load.
