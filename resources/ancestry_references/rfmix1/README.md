# RFMix v1 Reference Directory

This directory contains reference data in RFMix v1 (legacy) format.

## Directory Structure

```
rfmix1/
├── rfmix1_reference_chr1.alleles        # Binary alleles for chr1
├── rfmix1_reference_chr1.classes        # Population classes for chr1
├── rfmix1_reference_chr1.snp_locations  # Genetic positions (cM) for chr1
├── rfmix1_reference_chr2.alleles
├── rfmix1_reference_chr2.classes
├── rfmix1_reference_chr2.snp_locations
└── ... (per chromosome 1-22)
```

## File Formats

### .alleles
Binary alleles file, one haplotype per line:
- Each line: space-separated 0/1 values for each SNP
- Two lines per diploid sample (one per haplotype)

### .classes
Population labels file:
- Space-separated population codes (one per haplotype)
- Example: `AFR AFR EUR EUR EAS EAS ...`

### .snp_locations
Genetic map positions in centiMorgans:
- One value per line (or space-separated on single line)
- Must match SNP count in .alleles file

## Creating Reference Files

Use the convert_vcf_to_rfmix1.py helper script:

```bash
for chr in {1..22}; do
    python3 helper_scripts/convert_vcf_to_rfmix1.py \
        --vcf /path/to/reference_chr${chr}.vcf.gz \
        --sample-map /path/to/sample_populations.txt \
        --genetic-map /path/to/genetic_maps/chr${chr}.map \
        --output-prefix resources/ancestry_references/rfmix1/rfmix1_reference \
        --chromosome ${chr} \
        --verbose
done
```

## References

- RFMix v1: https://sites.google.com/site/rfmixlocalancestryinference/
- armartin/ancestry_pipeline: https://github.com/armartin/ancestry_pipeline

## Note

RFMix v1 is a legacy format. For new analyses, consider using:
- **RFMix v2** (default) - Modern implementation with better accuracy
- **FLARE** - Fast, memory-efficient
- **G-NOMIX** - Neural network-based, fastest option

These tools use the shared reference in `../lai_reference/`.
