# Ancestry Reference Panel Directory

This directory contains ancestry reference panels for Module 7 (Ancestry Estimation).

## Required Files for Global Ancestry

| File | Description |
|------|-------------|
| `reference_panel.rds` | Main reference panel for ADMIXTURE and global ancestry |

## Required Files for Local Ancestry Inference (LAI)

### RFMix v2 (Default LAI Tool)

| File | Description |
|------|-------------|
| `reference_haplotypes.vcf.gz` | Phased reference haplotypes from multiple populations |
| `reference_sample_map.txt` | Tab-separated file mapping samples to populations |

**Sample Map Format:**
```
sample1    EUR
sample2    AFR
sample3    EAS
```

### FLARE

| File | Description |
|------|-------------|
| `reference_haplotypes.vcf.gz` | Same as RFMix v2 (can share) |
| `flare_panels.txt` | Two-column file: sample_id → panel_name |

**Panels File Format:**
```
sample1    EUR
sample2    AFR
```

### G-NOMIX

| File | Description |
|------|-------------|
| `gnomix_models/` | Directory containing pre-trained models per chromosome |

Pre-trained models can be downloaded from:
- https://github.com/AI-sandbox/gnomix (see pretrained_gnomix_models/)

Alternatively, provide reference data to train new models:
- Reference VCF and sample map (same format as RFMix v2)

### RFMix v1 (Legacy)

| File | Description |
|------|-------------|
| `rfmix1_reference_alleles.txt` | Binary alleles file (0/1 format) |
| `rfmix1_reference_classes.txt` | Population labels for each haplotype |

**Note:** RFMix v1 uses a different input format than v2. Converting from VCF to binary alleles is handled automatically.

## Genetic Maps

Genetic maps should be placed in `resources/genetic_maps/`:

| File | Description |
|------|-------------|
| `chr1.genetic_map.txt` | Genetic map for chromosome 1 |
| `chr2.genetic_map.txt` | Genetic map for chromosome 2 |
| ... | |
| `chr22.genetic_map.txt` | Genetic map for chromosome 22 |

**Genetic Map Format (tab-separated):**
```
chromosome    position    rate    cM
22            16050075    0       0.0
22            16050115    0.0001  0.00001
```

## Reference Sources

### 1000 Genomes Project
- Website: https://www.internationalgenome.org/data
- Phased haplotypes for 26 populations
- Use GRCh38/hg38 version

### Human Genome Diversity Project (HGDP)
- Website: https://www.internationalgenome.org/data-portal/data-collection/hgdp
- 54 worldwide populations

### gnomAD
- Website: https://gnomad.broadinstitute.org/
- Large-scale reference for diverse populations

### TOPMed Freeze 8
- High-depth WGS reference panel
- Contact dbGaP for access

## Creating Reference Files

### From 1000 Genomes VCF:

```bash
# Download 1KG phase 3 data
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# Create sample map from population info
awk -F'\t' 'NR>1 {print $1"\t"$7}' integrated_call_samples_v3.20130502.ALL.panel > reference_sample_map.txt

# Subset to reference samples only
bcftools view -S reference_samples.txt ALL.chr22.vcf.gz -Oz -o reference_haplotypes_chr22.vcf.gz
```

### For G-NOMIX Training:

```bash
# If training new models, ensure you have:
# 1. Reference VCF with phased haplotypes
# 2. Sample map file
# 3. Genetic map

python3 gnomix.py \
    query.vcf.gz \
    genetic_map.txt \
    output_dir/ \
    22 \
    False \
    reference.vcf.gz \
    sample_map.txt
```

## Validation

Before running LAI, verify:
1. Reference and query VCFs have overlapping variants
2. Sample map contains all reference samples
3. Genetic map covers the chromosome range
4. Reference data is phased (using Beagle, SHAPEIT, or Eagle)

## Memory Requirements

| Tool | Memory (per chromosome) |
|------|-------------------------|
| RFMix v2 | 8-16 GB |
| FLARE | 4-8 GB |
| G-NOMIX | 4-8 GB |
| RFMix v1 | 8-16 GB |
