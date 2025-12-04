#!/usr/bin/env python3
"""
Convert VCF to RFMix v1 Input Format

This script converts a phased VCF file to the three input files required by RFMix v1:
  1. .alleles   - Binary alleles (0/1) for each haplotype, one line per haplotype
  2. .classes   - Population labels for each haplotype (space-separated)
  3. .snp_locations - Genetic positions in cM for each marker

Reference: https://sites.google.com/site/rfmixlocalancestryinference/
See also: https://github.com/armartin/ancestry_pipeline

Usage:
    python convert_vcf_to_rfmix1.py \
        --vcf reference_haplotypes.vcf.gz \
        --sample-map sample_populations.txt \
        --genetic-map genetic_map.txt \
        --output-prefix rfmix1_reference \
        --chromosome 1

Author: Genotyping Pipeline
"""

import argparse
import gzip
import sys
from collections import OrderedDict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert VCF to RFMix v1 format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Output files:
  {prefix}_chr{chr}.alleles       Binary alleles (one haplotype per line)
  {prefix}_chr{chr}.classes       Population labels (space-separated)
  {prefix}_chr{chr}.snp_locations Genetic positions in cM

Sample map format (tab-separated):
  SampleID    Population
  NA19700     AFR
  NA12878     EUR

Genetic map format (tab or space-separated, with header):
  position    rate_cM_per_Mb    genetic_position_cM
  16050       0.0147            0.0000

  OR (PLINK format):
  chr    position    genetic_position_cM    rsid
        """
    )

    parser.add_argument("--vcf", "-v", required=True,
                        help="Input phased VCF file (can be gzipped)")
    parser.add_argument("--sample-map", "-s", required=True,
                        help="Sample to population mapping file (tab-separated)")
    parser.add_argument("--genetic-map", "-g", required=True,
                        help="Genetic map file for interpolating cM positions")
    parser.add_argument("--output-prefix", "-o", required=True,
                        help="Output file prefix")
    parser.add_argument("--chromosome", "-c", required=True,
                        help="Chromosome to process")
    parser.add_argument("--min-maf", type=float, default=0.0,
                        help="Minimum minor allele frequency filter [default: 0.0]")
    parser.add_argument("--populations", "-p", nargs="+", default=None,
                        help="Populations to include (default: all)")
    parser.add_argument("--verbose", action="store_true",
                        help="Print progress messages")

    return parser.parse_args()


def load_sample_map(filename):
    """
    Load sample to population mapping.
    Returns dict: sample_id -> population
    """
    sample_pop = OrderedDict()

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) >= 2:
                sample_id = parts[0]
                population = parts[1]
                sample_pop[sample_id] = population

    return sample_pop


def load_genetic_map(filename):
    """
    Load genetic map for interpolating genetic positions.
    Returns list of tuples: [(physical_pos, genetic_pos_cM), ...]

    Supports multiple formats:
    - HapMap format: position rate_cM_per_Mb genetic_position_cM
    - PLINK format: chr position genetic_position_cM rsid
    - Simple format: position genetic_position_cM
    """
    positions = []

    with open(filename, 'r') as f:
        header = f.readline().strip().lower()

        # Detect format from header
        if 'rate' in header or 'cm_per_mb' in header.replace(' ', '_'):
            # HapMap format: position, rate, genetic_pos
            format_type = 'hapmap'
        elif header.startswith('chr') or 'chromosome' in header:
            # PLINK format: chr, position, genetic_pos, rsid
            format_type = 'plink'
        else:
            # Simple format or no header
            format_type = 'simple'
            # Re-read as data if first line looks like data
            if header.replace('.', '').replace('\t', '').replace(' ', '').isdigit():
                parts = header.split()
                if len(parts) >= 2:
                    positions.append((int(parts[0]), float(parts[-1])))

        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()

            try:
                if format_type == 'hapmap':
                    # position, rate, genetic_pos
                    phys_pos = int(parts[0])
                    gen_pos = float(parts[2])
                elif format_type == 'plink':
                    # chr, position, genetic_pos, rsid
                    phys_pos = int(parts[1])
                    gen_pos = float(parts[2])
                else:
                    # position, genetic_pos
                    phys_pos = int(parts[0])
                    gen_pos = float(parts[-1])

                positions.append((phys_pos, gen_pos))
            except (ValueError, IndexError):
                continue

    # Sort by physical position
    positions.sort(key=lambda x: x[0])

    return positions


def interpolate_genetic_position(phys_pos, genetic_map):
    """
    Interpolate genetic position (cM) for a given physical position.
    Uses linear interpolation between flanking markers.
    """
    if not genetic_map:
        return 0.0

    # Binary search for flanking positions
    left = 0
    right = len(genetic_map) - 1

    # Handle edge cases
    if phys_pos <= genetic_map[0][0]:
        return genetic_map[0][1]
    if phys_pos >= genetic_map[-1][0]:
        return genetic_map[-1][1]

    # Find flanking markers
    while right - left > 1:
        mid = (left + right) // 2
        if genetic_map[mid][0] <= phys_pos:
            left = mid
        else:
            right = mid

    # Linear interpolation
    pos1, gen1 = genetic_map[left]
    pos2, gen2 = genetic_map[right]

    if pos2 == pos1:
        return gen1

    fraction = (phys_pos - pos1) / (pos2 - pos1)
    return gen1 + fraction * (gen2 - gen1)


def process_vcf(vcf_file, sample_map, genetic_map, chromosome,
                output_prefix, min_maf=0.0, populations=None, verbose=False):
    """
    Process VCF and write RFMix v1 format files.
    """
    # Open VCF
    if vcf_file.endswith('.gz'):
        vcf_handle = gzip.open(vcf_file, 'rt')
    else:
        vcf_handle = open(vcf_file, 'r')

    # Parse header to get sample order
    vcf_samples = []
    for line in vcf_handle:
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            fields = line.strip().split('\t')
            vcf_samples = fields[9:]
            break

    if verbose:
        print(f"Found {len(vcf_samples)} samples in VCF")

    # Filter samples by population
    selected_samples = []
    selected_pops = []
    selected_indices = []

    for i, sample in enumerate(vcf_samples):
        if sample in sample_map:
            pop = sample_map[sample]
            if populations is None or pop in populations:
                selected_samples.append(sample)
                selected_pops.append(pop)
                selected_indices.append(i)

    if verbose:
        print(f"Selected {len(selected_samples)} samples from {len(set(selected_pops))} populations")
        pop_counts = {}
        for p in selected_pops:
            pop_counts[p] = pop_counts.get(p, 0) + 1
        for p, c in sorted(pop_counts.items()):
            print(f"  {p}: {c} samples ({c*2} haplotypes)")

    if not selected_samples:
        print("ERROR: No samples found matching sample map and population filter",
              file=sys.stderr)
        sys.exit(1)

    # Initialize haplotype storage
    # Each sample has 2 haplotypes, so we have n_samples * 2 haplotypes
    n_haplotypes = len(selected_samples) * 2
    haplotypes = [[] for _ in range(n_haplotypes)]
    snp_positions = []

    # Process variants
    n_variants = 0
    n_skipped_maf = 0
    n_skipped_multiallelic = 0

    for line in vcf_handle:
        if line.startswith('#'):
            continue

        fields = line.strip().split('\t')
        chrom = fields[0].replace('chr', '')

        # Skip if not target chromosome
        if chrom != str(chromosome).replace('chr', ''):
            continue

        pos = int(fields[1])
        ref = fields[3]
        alt = fields[4]

        # Skip multiallelic sites
        if ',' in alt:
            n_skipped_multiallelic += 1
            continue

        # Get genotypes for selected samples
        genotypes = fields[9:]

        # Calculate allele frequency for MAF filter
        alt_count = 0
        total_count = 0

        for idx in selected_indices:
            gt = genotypes[idx].split(':')[0]
            if '|' in gt:
                a1, a2 = gt.split('|')
            elif '/' in gt:
                a1, a2 = gt.split('/')
            else:
                continue

            if a1 in ('0', '1'):
                alt_count += int(a1)
                total_count += 1
            if a2 in ('0', '1'):
                alt_count += int(a2)
                total_count += 1

        if total_count > 0:
            af = alt_count / total_count
            maf = min(af, 1 - af)

            if maf < min_maf:
                n_skipped_maf += 1
                continue

        # Get genetic position
        gen_pos = interpolate_genetic_position(pos, genetic_map)
        snp_positions.append(gen_pos)

        # Extract haplotypes
        hap_idx = 0
        for idx in selected_indices:
            gt = genotypes[idx].split(':')[0]

            if '|' in gt:
                a1, a2 = gt.split('|')
            elif '/' in gt:
                a1, a2 = gt.split('/')
            else:
                a1, a2 = '.', '.'

            # Convert to binary (0=ref, 1=alt, 9=missing)
            def to_binary(a):
                if a == '0':
                    return '0'
                elif a == '1':
                    return '1'
                else:
                    return '9'

            haplotypes[hap_idx].append(to_binary(a1))
            haplotypes[hap_idx + 1].append(to_binary(a2))
            hap_idx += 2

        n_variants += 1

        if verbose and n_variants % 10000 == 0:
            print(f"  Processed {n_variants} variants...")

    vcf_handle.close()

    if verbose:
        print(f"\nVariant summary:")
        print(f"  Total variants processed: {n_variants}")
        print(f"  Skipped (multiallelic): {n_skipped_multiallelic}")
        print(f"  Skipped (MAF < {min_maf}): {n_skipped_maf}")

    # Write output files
    chr_suffix = f"_chr{chromosome}"

    # 1. Write alleles file (one haplotype per line)
    alleles_file = f"{output_prefix}{chr_suffix}.alleles"
    if verbose:
        print(f"\nWriting {alleles_file}...")

    with open(alleles_file, 'w') as f:
        for hap in haplotypes:
            f.write(' '.join(hap) + '\n')

    # 2. Write classes file (population labels, 2 per sample for each haplotype)
    classes_file = f"{output_prefix}{chr_suffix}.classes"
    if verbose:
        print(f"Writing {classes_file}...")

    with open(classes_file, 'w') as f:
        # Each sample has 2 haplotypes, so repeat each population label twice
        labels = []
        for pop in selected_pops:
            labels.extend([pop, pop])
        f.write(' '.join(labels) + '\n')

    # 3. Write SNP locations file (genetic positions in cM)
    snp_loc_file = f"{output_prefix}{chr_suffix}.snp_locations"
    if verbose:
        print(f"Writing {snp_loc_file}...")

    with open(snp_loc_file, 'w') as f:
        for pos in snp_positions:
            f.write(f"{pos:.6f}\n")

    if verbose:
        print(f"\nConversion complete!")
        print(f"  Haplotypes: {n_haplotypes}")
        print(f"  Variants: {n_variants}")
        print(f"  Output files:")
        print(f"    - {alleles_file}")
        print(f"    - {classes_file}")
        print(f"    - {snp_loc_file}")


def main():
    args = parse_args()

    if args.verbose:
        print("=" * 60)
        print("VCF to RFMix v1 Format Converter")
        print("=" * 60)
        print(f"Input VCF: {args.vcf}")
        print(f"Sample map: {args.sample_map}")
        print(f"Genetic map: {args.genetic_map}")
        print(f"Chromosome: {args.chromosome}")
        print(f"Output prefix: {args.output_prefix}")
        if args.populations:
            print(f"Populations: {', '.join(args.populations)}")
        print("=" * 60)

    # Load reference data
    if args.verbose:
        print("\nLoading sample map...")
    sample_map = load_sample_map(args.sample_map)
    if args.verbose:
        print(f"  Loaded {len(sample_map)} samples")

    if args.verbose:
        print("\nLoading genetic map...")
    genetic_map = load_genetic_map(args.genetic_map)
    if args.verbose:
        print(f"  Loaded {len(genetic_map)} positions")

    # Process VCF
    if args.verbose:
        print("\nProcessing VCF...")

    process_vcf(
        vcf_file=args.vcf,
        sample_map=sample_map,
        genetic_map=genetic_map,
        chromosome=args.chromosome,
        output_prefix=args.output_prefix,
        min_maf=args.min_maf,
        populations=args.populations,
        verbose=args.verbose
    )


if __name__ == "__main__":
    main()
