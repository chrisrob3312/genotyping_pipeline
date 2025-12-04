# GWAS Simulation Guide for Benchmarking

This guide covers software and methods for simulating phenotypes with controlled genetic architecture to test imputation approaches fairly.

## Why Simulate Instead of Using Published GWAS?

Published GWAS effect sizes are **EUR-biased** because:
1. ~80% of GWAS participants are European ancestry
2. Causal variants may be in different LD with tag SNPs across ancestries
3. Effect sizes may differ by ancestry (G×E, G×G interactions)
4. Ancestry-specific causal variants are missed entirely

**Solution:** Simulate ground truth where you control the genetic architecture.

---

## Recommended Tools

### 1. GCTA-GREML (Most Widely Used)

**Best for:** Simulating phenotypes from real genotypes with specified heritability.

```bash
# Installation
wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip
unzip gcta-*.zip

# Simulate phenotype with h2 = 0.5, 1000 causal variants
gcta64 --bfile 1kg_data \
    --simu-qt \
    --simu-causal-loci causal_snps.txt \
    --simu-hsq 0.5 \
    --simu-k 0.1 \
    --out simulated_pheno

# causal_snps.txt format (SNP_ID, effect_size):
# rs123456  0.05
# rs789012  -0.03
```

**Documentation:** https://yanglab.westlake.edu.cn/software/gcta/#GREMLpowerCalculator

---

### 2. PhenotypeSimulator (R Package)

**Best for:** Complex trait simulation with multiple components (genetics, covariates, noise).

```r
# Installation
install.packages("PhenotypeSimulator")

library(PhenotypeSimulator)

# Simulate with ancestry-specific effects
genotypes <- read_genotypes("1kg_data")

# Define genetic effects with population structure
pheno <- simulatePhenotypes(
    N = nrow(genotypes),
    P = 1,  # 1 phenotype
    genVar = 0.5,  # Genetic variance (h2)
    noiseVar = 0.5,
    genotypes = genotypes,
    # Add population stratification
    pcaEffects = TRUE,
    pcs = ancestry_pcs
)
```

**Documentation:** https://cran.r-project.org/package=PhenotypeSimulator

---

### 3. msprime + SLiM (Population Genetics Simulation)

**Best for:** Simulating genotypes with realistic demography and selection, including admixed populations.

```python
# Installation
pip install msprime

import msprime

# Simulate African-European admixture (like AMR populations)
demography = msprime.Demography()
demography.add_population(name="AFR", initial_size=14000)
demography.add_population(name="EUR", initial_size=10000)
demography.add_population(name="AMR", initial_size=1000)
demography.add_population(name="ANC", initial_size=14000)

# Admixture event
demography.add_admixture(
    time=10,  # 10 generations ago
    derived="AMR",
    ancestral=["AFR", "EUR"],
    proportions=[0.3, 0.7]  # 30% AFR, 70% EUR ancestry
)

# Out of Africa
demography.add_population_split(
    time=2000, derived=["AFR", "EUR"], ancestral="ANC"
)

# Simulate
ts = msprime.sim_ancestry(
    samples={"AFR": 500, "EUR": 500, "AMR": 500},
    demography=demography,
    sequence_length=1e6,
    recombination_rate=1e-8
)

# Add mutations
ts = msprime.sim_mutations(ts, rate=1e-8)
```

**Documentation:** https://tskit.dev/msprime/docs/stable/

---

### 4. simplePHENOTYPES (Bioconductor)

**Best for:** Multi-trait simulation with epistasis and pleiotropy.

```r
# Installation
BiocManager::install("simplePHENOTYPES")

library(simplePHENOTYPES)

# Create phenotype with ancestry-specific effects
create_phenotypes(
    geno_file = "genotypes.hmp.txt",
    add_QTN_num = 50,           # 50 additive QTN
    add_effect = 0.1,           # Effect size
    h2 = 0.5,                   # Heritability
    rep = 10,                   # 10 replicates
    model = "A",                # Additive model
    output_dir = "sim_output"
)
```

---

### 5. Tractor's Admixture Simulation

**Best for:** Simulating phenotypes where effects depend on local ancestry.

The Tractor team provides simulation scripts:

```bash
# Clone Tractor
git clone https://github.com/Atkinson-Lab/Tractor.git

# Their simulation approach (from paper supplements):
# 1. Simulate admixed genomes with known local ancestry
# 2. Assign causal effects that depend on LAI
# 3. Compare Tractor vs standard GWAS power
```

**Paper:** Atkinson et al. 2021, Nature Genetics
**GitHub:** https://github.com/Atkinson-Lab/Tractor

---

## Recommended Simulation Scenarios for Your Benchmarking

### Scenario 1: AFR-Enriched Variants

Tests whether your pipeline recovers variants common in AFR but rare in EUR.

```bash
# Using GCTA with custom causal SNP list
# Select SNPs where AFR_AF > 0.1 and EUR_AF < 0.05

# Get AFR-enriched SNPs from gnomAD
bcftools query -f '%ID\t%INFO/AF_afr\t%INFO/AF_nfe\n' gnomad.vcf.gz | \
    awk '$2 > 0.1 && $3 < 0.05 {print $1, 0.05}' > afr_enriched_causal.txt

# Simulate
gcta64 --bfile 1kg_data \
    --simu-qt \
    --simu-causal-loci afr_enriched_causal.txt \
    --simu-hsq 0.5 \
    --out pheno_afr_enriched
```

### Scenario 2: Ancestry-Specific Effect Sizes

Same variant, different effect by ancestry (tests Tractor value).

```python
# Using custom Python simulation
import numpy as np
import pandas as pd

def simulate_ancestry_specific(genotypes, ancestry, n_causal=50, h2=0.5):
    """
    Simulate phenotype where effect size depends on ancestry.
    """
    n_samples, n_snps = genotypes.shape

    # Select causal variants
    causal_idx = np.random.choice(n_snps, n_causal, replace=False)

    # Base effects
    beta_base = np.random.normal(0, 1, n_causal)

    # Ancestry-specific multipliers
    ancestry_mult = {
        'AFR': np.random.uniform(0.5, 1.5, n_causal),
        'EUR': np.ones(n_causal),  # Reference
        'AMR': np.random.uniform(0.7, 1.3, n_causal)
    }

    # Calculate genetic values
    G = np.zeros(n_samples)
    for i, idx in enumerate(causal_idx):
        for j in range(n_samples):
            anc = ancestry[j]
            mult = ancestry_mult.get(anc, 1.0)[i]
            G[j] += genotypes[j, idx] * beta_base[i] * mult

    # Add noise for target h2
    var_G = np.var(G)
    var_E = var_G * (1 - h2) / h2
    E = np.random.normal(0, np.sqrt(var_E), n_samples)

    return G + E
```

### Scenario 3: Local Ancestry-Dependent Effects

Effects only manifest on specific local ancestry background.

```python
# Requires phased genotypes + LAI calls
def simulate_lai_dependent(haplotypes, local_ancestry, n_causal=30, h2=0.4):
    """
    Effect only present when variant is on AFR local ancestry tract.
    Standard GWAS will see diluted effect.
    Tractor GWAS should detect full ancestry-specific effect.
    """
    n_samples = haplotypes.shape[0] // 2  # Diploid
    n_snps = haplotypes.shape[1]

    causal_idx = np.random.choice(n_snps, n_causal, replace=False)
    beta = np.random.normal(0, 1, n_causal)

    G = np.zeros(n_samples)
    for i, idx in enumerate(causal_idx):
        for j in range(n_samples):
            # Haplotype 1
            if haplotypes[2*j, idx] == 1:  # Has alt allele
                if local_ancestry[2*j, idx] == 0:  # On AFR tract
                    G[j] += beta[i]
                # No effect if on EUR/NAT tract

            # Haplotype 2
            if haplotypes[2*j+1, idx] == 1:
                if local_ancestry[2*j+1, idx] == 0:
                    G[j] += beta[i]

    # Add noise
    var_G = np.var(G) if np.var(G) > 0 else 1
    var_E = var_G * (1 - h2) / h2
    E = np.random.normal(0, np.sqrt(var_E), n_samples)

    return G + E
```

### Scenario 4: Rare Variant Effects

Tests imputation quality for rare variants.

```bash
# Select rare variants (MAF 0.1-1%)
bcftools view -q 0.001:minor -Q 0.01:minor 1kg.vcf.gz -Oz -o rare_variants.vcf.gz

# Extract rare variant IDs
bcftools query -f '%ID\n' rare_variants.vcf.gz > rare_snp_list.txt

# Assign larger effects to rare variants (observed in nature)
awk '{print $1, (rand()-0.5)*0.2}' rare_snp_list.txt | head -100 > rare_causal.txt

# Simulate
gcta64 --bfile 1kg_data \
    --simu-qt \
    --simu-causal-loci rare_causal.txt \
    --simu-hsq 0.3 \
    --out pheno_rare
```

---

## Benchmarking Workflow

```
┌─────────────────────────────────────────────────────────────────┐
│  1. SIMULATE GROUND TRUTH                                       │
│     - Choose scenario (AFR-enriched, LAI-specific, etc.)        │
│     - Generate phenotypes with known causal variants            │
│     - Record: causal SNPs, effect sizes, ancestry effects       │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  2. IMPUTE WITH DIFFERENT APPROACHES                            │
│     - Approach A: Traditional (Intersect → QC → Impute)         │
│     - Approach F: Ours (Union → Impute → MagicalRsq-X → QC)     │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  3. RUN GWAS                                                    │
│     - Standard GWAS on both imputed datasets                    │
│     - Tractor GWAS (if LAI available)                           │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  4. COMPARE POWER                                               │
│     - % causal variants reaching significance                   │
│     - Effect size correlation with truth                        │
│     - Power by ancestry group                                   │
│     - Standard vs Tractor comparison                            │
└─────────────────────────────────────────────────────────────────┘
```

---

## Key References

1. **GCTA:** Yang et al. 2011, AJHG - https://doi.org/10.1016/j.ajhg.2010.11.011
2. **Tractor:** Atkinson et al. 2021, Nature Genetics - https://doi.org/10.1038/s41588-020-00766-y
3. **msprime:** Kelleher et al. 2016, PLoS Comput Biol - https://doi.org/10.1371/journal.pcbi.1004842
4. **PRS portability:** Martin et al. 2019, Nature Genetics - https://doi.org/10.1038/s41588-019-0379-x
5. **PAGE Study:** Wojcik et al. 2019, Nature - https://doi.org/10.1038/s41586-019-1310-4

---

## Quick Start Commands

```bash
# 1. Install GCTA
wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip
unzip gcta-*.zip && chmod +x gcta_1.94.1*

# 2. Download 1KG data (if not already)
./benchmarking/download_benchmark_data.sh -o benchmarking/test_data

# 3. Create causal SNP list (AFR-enriched example)
# (Manually or using gnomAD frequencies)
echo "rs123456 0.05" > causal_afr.txt
echo "rs789012 -0.03" >> causal_afr.txt

# 4. Simulate
gcta64 --bfile benchmarking/test_data/genotypes/1kg_omni \
    --simu-qt \
    --simu-causal-loci causal_afr.txt \
    --simu-hsq 0.5 \
    --out benchmarking/test_data/phenotypes/sim_afr_enriched

# 5. Run benchmarking with simulated phenotype
# (Use your pipeline with the simulated .phen file)
```
