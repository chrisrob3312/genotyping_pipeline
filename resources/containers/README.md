# Container Directory

This directory holds Apptainer/Singularity container images (.sif files) for the pipeline.

## Required Containers

The following containers will be built by Module 0 or need to be placed here:

- `plink.sif` - PLINK 1.9 and 2.0
- `perl_vcftools.sif` - Perl, VCFtools, CrossMap
- `r_genetics.sif` - R with genetics packages (GENESIS, MagicalRsq)
- `python_tools.sif` - Python with imputationbot, terralab-cli
- `ancestry_tools.sif` - ADMIXTURE, RFMix, FLARE, G-NOMIX, GRAF-anc

## Building Containers

Run Module 0 to build all containers:

```bash
nextflow run modules/Module0_Apptainer_Build.nf -c config/nextflow_module0.config
```
