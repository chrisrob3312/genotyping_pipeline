# GRAF-anc Data Directory

This directory contains the reference data files required for GRAF-anc ancestry estimation.

## Required Files

| File | Description |
|------|-------------|
| `AncSnpPopAFs.txt` | Ancestry SNP population allele frequencies |

## Downloading GRAF-anc Data

Download the required data files from the official GRAF-anc repository:

```bash
# Clone GRAF-anc repository
git clone https://github.com/jimmy-penn/grafanc.git

# Copy data file to this directory
cp grafanc/data/AncSnpPopAFs.txt resources/grafanc_data/
```

## File Format

### AncSnpPopAFs.txt

Tab-separated file containing:
- SNP identifiers
- Population-specific allele frequencies for ancestry-informative markers
- Used to calculate genetic distance axes (GD1-GD6)

## Building GRAF-anc

GRAF-anc requires htslib. To build from source:

```bash
cd grafanc

# Install dependencies (Ubuntu/Debian)
sudo apt-get install libhts-dev

# Build
make

# Test
./test_grafanc.pl
```

## Container Usage

The pipeline's `ancestry_suite.sif` container includes GRAF-anc with required dependencies.
Data files must be mounted or included in the container.

## Reference

- GitHub: https://github.com/jimmy-penn/grafanc
- GRAF-anc uses 3-digit ancestry group codes:
  - 1XX = African
  - 2XX = Middle East/North Africa
  - 3XX = European
  - 4XX = South Asian
  - 5XX = East Asian
  - 6XX = American
  - 7XX = Oceania
  - 8XX = Multi-ancestry
