# Reference Files Directory

This directory should contain reference genome files and variant references for the pipeline.

## Required Files

### Reference Genomes
- `hg19.fa` - Human reference genome GRCh37/hg19
- `hg19.fa.fai` - FAI index for hg19 reference
- `hg38.fa` - Human reference genome GRCh38/hg38
- `hg38.fa.fai` - FAI index for hg38 reference

### Liftover Chain File
- `hg19ToHg38.over.chain.gz` - Chain file for lifting hg19 coordinates to hg38

### TOPMed Reference
- `PASS.Variants.TOPMed_freeze10_hg38.tab.gz` - TOPMed variant reference for strand checking

## Download Sources

### Reference Genomes
Download from UCSC or NCBI:
- hg19: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
- hg38: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

After downloading, decompress and index:
```bash
gunzip hg19.fa.gz hg38.fa.gz
samtools faidx hg19.fa
samtools faidx hg38.fa
```

### Liftover Chain
Download from UCSC:
```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
```

### TOPMed Reference
Download from TOPMed:
```bash
# From TOPMed Imputation Server reference files
wget https://bravo.sph.umich.edu/freeze5/hg38/PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz
# Or use freeze 10 if available
```

## Usage

These files are automatically used by the pipeline through nextflow.config parameters:
- `params.hg19_fasta`
- `params.hg38_fasta`
- `params.liftover_chain`
- `params.topmed_reference`
