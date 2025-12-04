# LAI Reference Directory

This directory contains shared reference data for Local Ancestry Inference (LAI) tools:
- **RFMix v2** (default)
- **FLARE**
- **G-NOMIX**

## Directory Structure

```
lai_reference/
├── reference_haplotypes.vcf.gz      # Phased VCF (shared by all tools)
├── reference_haplotypes.vcf.gz.tbi  # VCF index
├── rfmix2_sample_map.txt            # RFMix v2 sample map (tab-separated)
├── flare_panels.txt                 # FLARE panel file (space-separated)
└── gnomix_sample_map.tsv            # G-NOMIX sample map (tab-separated with header)
```

## File Formats

### reference_haplotypes.vcf.gz
Phased VCF containing reference panel haplotypes. Used by all three tools.

### rfmix2_sample_map.txt
Tab-separated file (no header):
```
SampleID	Population
NA19700	YRI
NA12878	CEU
HG00403	CHS
```

### flare_panels.txt
Space-separated file (no header):
```
NA19700 YRI
NA12878 CEU
HG00403 CHS
```

### gnomix_sample_map.tsv
Tab-separated file **with header**:
```
Sample	Population
NA19700	YRI
NA12878	CEU
HG00403	CHS
```

## Creating Reference Files

Use the format_reference_panel.sh helper script:

```bash
./helper_scripts/format_reference_panel.sh \
    --input /path/to/reference.vcf.gz \
    --sample-map /path/to/sample_populations.txt \
    --genetic-map-dir /path/to/genetic_maps/ \
    --output-dir /path/to/formatted_references

# Then copy to pipeline:
cp -r /path/to/formatted_references/lai_reference/* resources/ancestry_references/lai_reference/
```

## Recommended Reference Panels

1. **HGDP-1KG** (Human Genome Diversity Project + 1000 Genomes)
   - ~4000 samples from diverse global populations
   - Publicly available phased haplotypes

2. **Custom panels** (e.g., HGDP-1KG + MX Biobank)
   - Add population-specific samples for improved accuracy
   - Ensure phasing quality matches reference

## See Also

- `../rfmix1/` - RFMix v1 format files (legacy, per-chromosome)
- `../gnomix_models/` - Pre-trained G-NOMIX models
- `../../genetic_maps/` - Genetic recombination maps
