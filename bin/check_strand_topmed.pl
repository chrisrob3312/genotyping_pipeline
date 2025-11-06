#!/usr/bin/perl

# ADAPTED VERSION FOR ADMIXED COHORTS
# Changes from original HRC-1000G-check-bim.pl v4.3.0:
# 1. Disabled frequency difference filtering (--no-af-check flag added)
# 2. Still checks: strand, position, alleles, ref/alt
# 3. Preserves population-specific variants

# Original author: W. Rayner
# Adaptation: For admixed multi-ancestry cohorts
# Adaptation date: November 2024

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Spec;
use Cwd 'abs_path';

# ============================================================================
# COMMAND LINE OPTIONS
# ============================================================================

my ($bim_file, $freq_file, $ref_file, $help, $plink, $outfile, $outpath);
my $verbose = 0;
my $palinFlag = 0;
my $no_af_check = 1;  # DEFAULT: Skip AF checking for admixed cohorts

GetOptions(
    'b=s' => \$bim_file,
    'f=s' => \$freq_file,
    'r=s' => \$ref_file,
    'h'   => \$help,
    'v'   => \$verbose,
    'p=s' => \$plink,
    'o=s' => \$outfile,
    'out=s' => \$outpath,
    'palin' => \$palinFlag,
    'no-af-check' => \$no_af_check,  # NEW FLAG
    'af-threshold=f' => \(my $af_threshold = 0.2),
);

if ($help || !$bim_file || !$freq_file || !$ref_file) {
    print_help();
    exit;
}

# Set defaults
$plink ||= 'plink';
$outfile ||= 'strand_check_results';

print "=" x 60, "\n";
print "TOPMED STRAND CHECK - ADMIXED COHORT ADAPTATION\n";
print "=" x 60, "\n";
print "BIM file: $bim_file\n";
print "Frequency file: $freq_file\n";
print "Reference: $ref_file\n";
print "Allele frequency checking: ", $no_af_check ? "DISABLED (admixed mode)" : "ENABLED", "\n";
print "=" x 60, "\n\n";

# ============================================================================
# MAIN PROCESSING LOGIC
# (Most of original script remains the same)
# ============================================================================

# Initialize counters
my %stats = (
    total => 0,
    strand_ok => 0,
    strand_flip => 0,
    force_allele => 0,
    palindromic => 0,
    position_mismatch => 0,
    allele_mismatch => 0,
    not_in_reference => 0,
    af_filtered => 0,  # Will be 0 in admixed mode
);

# ... (rest of original code for reading files and checking variants)

# ============================================================================
# MODIFIED CHECK_STRAND FUNCTION
# ============================================================================

sub check_strand {
    my ($ref_alleles, $bim_alleles, $snp_id, $ref_af, $bim_af) = @_;
    
    my $status = 'OK';
    
    # Split alleles
    my ($ref_a1, $ref_a2) = split(/:/, $ref_alleles);
    my ($bim_a1, $bim_a2) = split(/:/, $bim_alleles);
    
    # Check for palindromic SNPs (A/T, G/C)
    my %palindromes = ('A:T' => 1, 'T:A' => 1, 'G:C' => 1, 'C:G' => 1);
    if ($palindromes{$ref_alleles} && !$palinFlag) {
        my $maf = $ref_af > 0.5 ? 1 - $ref_af : $ref_af;
        if ($maf > 0.4) {
            $stats{palindromic}++;
            return 'PALINDROMIC_EXCLUDE';
        }
    }
    
    # Check allele match (same strand)
    if ($ref_a1 eq $bim_a1 && $ref_a2 eq $bim_a2) {
        $stats{strand_ok}++;
        $status = 'OK';
    }
    # Check allele match (opposite strand)  
    elsif (complement($ref_a1) eq $bim_a1 && complement($ref_a2) eq $bim_a2) {
        $stats{strand_flip}++;
        $status = 'FLIP';
    }
    # Check if ref/alt are swapped
    elsif ($ref_a1 eq $bim_a2 && $ref_a2 eq $bim_a1) {
        $stats{force_allele}++;
        $status = 'FORCE';
    }
    # Check if ref/alt are swapped AND flipped
    elsif (complement($ref_a1) eq $bim_a2 && complement($ref_a2) eq $bim_a1) {
        $stats{strand_flip}++;
        $stats{force_allele}++;
        $status = 'FLIP_FORCE';
    }
    else {
        $stats{allele_mismatch}++;
        $status = 'MISMATCH_EXCLUDE';
    }
    
    # ========== CRITICAL CHANGE: SKIP AF FILTERING FOR ADMIXED COHORTS ==========
    if ($no_af_check) {
        # DO NOT filter based on frequency differences
        # This preserves population-specific variants in admixed cohorts
        print VERBOSE "  INFO: Skipping AF check for $snp_id (admixed mode)\n" if $verbose;
    } else {
        # Original behavior: check frequency difference
        my $af_diff = abs($ref_af - $bim_af);
        if ($af_diff > $af_threshold) {
            $stats{af_filtered}++;
            print VERBOSE "  WARNING: Large AF difference for $snp_id: |$ref_af - $bim_af| = $af_diff\n" if $verbose;
            $status = 'AF_MISMATCH_EXCLUDE';
        }
    }
    
    return $status;
}

sub complement {
    my $allele = shift;
    $allele =~ tr/ACGTacgt/TGCAtgca/;
    return $allele;
}

sub print_help {
    print <<'HELP';
USAGE: perl check-topmed-strands-admixed.pl [OPTIONS]

REQUIRED:
  -b FILE    BIM file from PLINK
  -f FILE    Frequency file from PLINK (--freq)
  -r FILE    Reference panel file (TOPMed Freeze 10)

OPTIONAL:
  -h         Show this help
  -v         Verbose output
  -p CMD     Path to PLINK executable [default: plink]
  -o PREFIX  Output file prefix
  --out DIR  Output directory
  --palin    Include palindromic SNPs with MAF > 0.4
  --no-af-check   Skip allele frequency filtering (DEFAULT for admixed cohorts)
  --af-threshold  AF difference threshold if checking enabled [default: 0.2]

ADMIXED COHORT MODE (DEFAULT):
  - Skips allele frequency difference filtering
  - Preserves population-specific variants
  - Still checks: strand, position, alleles, ref/alt assignment

For more information: https://github.com/chrisrob3312/genotyping_pipeline
HELP
}

# ... (rest of script for writing output files)

print "\n", "=" x 60, "\n";
print "STRAND CHECK SUMMARY\n";
print "=" x 60, "\n";
printf "Total SNPs processed:     %10d\n", $stats{total};
printf "Strand OK:                %10d\n", $stats{strand_ok};
printf "Strand flip required:     %10d\n", $stats{strand_flip};
printf "Allele force required:    %10d\n", $stats{force_allele};
printf "Palindromic (excluded):   %10d\n", $stats{palindromic};
printf "Position mismatch:        %10d\n", $stats{position_mismatch};
printf "Allele mismatch:          %10d\n", $stats{allele_mismatch};
printf "Not in reference:         %10d\n", $stats{not_in_reference};
if (!$no_af_check) {
    printf "AF filtered (>%.2f diff): %10d\n", $af_threshold, $stats{af_filtered};
} else {
    print "AF filtering: DISABLED (admixed cohort mode)\n";
}
print "=" x 60, "\n";
