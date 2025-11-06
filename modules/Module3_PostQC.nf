// Module 3: Post-Imputation QC
// UPDATED: Added VCF indexing and formatting at the beginning
// Uses R package MagicalRsq for variant quality filtering
// Parallelized across all platforms

// ============================================================================
// PROCESS 0: Index imputed VCF files from Module 2
// ============================================================================
process indexImputedVCFs {
    label 'bcftools'
    tag "${platform}_chr${chr}"
    publishDir "${params.outdir}/module3/00_indexed/${platform}", mode: 'copy'
    
    input:
    tuple val(platform), val(chr), path(vcf), path(info)
    
    output:
    tuple val(platform), val(chr), 
          path(vcf), 
          path("${vcf}.tbi"), 
          path(info),
          emit: indexed_vcfs
    
    script:
    """
    # Create tabix index for the VCF file
    # This is required because Module 2 imputation servers don't return indexed files
    echo "Indexing ${vcf}..."
    bcftools index -t ${vcf}
    
    # Verify index was created
    if [ ! -f "${vcf}.tbi" ]; then
        echo "ERROR: Failed to create index for ${vcf}"
        exit 1
    fi
    
    echo "Successfully indexed ${vcf}"
    """
}

// ============================================================================
// PROCESS 1: Filter variants with MagicalRsq (R² >= 0.3)
// ============================================================================
process magicalRsqFilter {
    label 'R'
    tag "${platform}_chr${chr}"
    publishDir "${params.outdir}/module3/01_rsq_filtered/${platform}", mode: 'copy'
    
    input:
    tuple val(platform), val(chr), path(vcf), path(tbi), path(info)
    
    output:
    tuple val(platform), val(chr),
          path("${platform}_chr${chr}_rsq0.3.vcf.gz"),
          path("${platform}_chr${chr}_rsq0.3.vcf.gz.tbi"),
          path("${platform}_chr${chr}_rsq_stats.txt"),
          emit: filtered_vcf
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(MagicalRsq)
    library(data.table)
    
    cat("Processing chromosome ${chr} for platform ${platform}\\n")
    cat("Input VCF: ${vcf}\\n")
    cat("Input INFO: ${info}\\n")
    
    # Read imputation info file
    info <- fread("${info}")
    cat("Read", nrow(info), "variants from info file\\n")
    
    # Calculate MagicalRsq (empirical R²)
    # This accounts for LD and provides better quality metric than standard Rsq
    # MagicalRsq-X features include:
    # - LD scores from TOP-LD (4 populations, long-range ±1Mb and short-range ±100kb)
    # - Recombination rates from 1000 Genomes
    # - Ancestry-specific MAF
    # - S/HIC features (selection/haplotype scores)
    # NO ANCESTRY INPUT REQUIRED - model uses all population LD scores automatically
    
    cat("Calculating MagicalRsq values...\\n")
    info\$magical_rsq <- calculate_magical_rsq(
        info\$Rsq,           # Server-reported R²
        info\$MAF,           # Minor allele frequency
        info\$Genotyped      # Whether variant was genotyped
    )
    
    # Filter: MagicalRsq >= 0.3
    pass_variants <- info[magical_rsq >= 0.3, ]
    cat("Variants passing filter:", nrow(pass_variants), "/", nrow(info), "\\n")
    
    # Save statistics
    stats <- data.frame(
        chromosome = "${chr}",
        platform = "${platform}",
        total_variants = nrow(info),
        pass_variants = nrow(pass_variants),
        filtered_variants = nrow(info) - nrow(pass_variants),
        pct_pass = 100 * nrow(pass_variants) / nrow(info),
        mean_rsq = mean(info\$Rsq),
        mean_magical_rsq = mean(info\$magical_rsq),
        mean_rsq_pass = mean(pass_variants\$magical_rsq)
    )
    write.table(stats, "${platform}_chr${chr}_rsq_stats.txt", 
                quote=FALSE, row.names=FALSE, sep="\\t")
    
    # Write variant IDs to keep
    write.table(pass_variants\$ID, "variants_to_keep.txt", 
                quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    # Filter VCF with bcftools
    cat("Filtering VCF with bcftools...\\n")
    system("bcftools view -i 'ID=@variants_to_keep.txt' ${vcf} -Oz -o ${platform}_chr${chr}_rsq0.3.vcf.gz")
    system("bcftools index -t ${platform}_chr${chr}_rsq0.3.vcf.gz")
    
    cat("\\n=== MagicalRsq Filtering Complete ===\\n")
    cat("Filtered from", nrow(info), "to", nrow(pass_variants), "variants\\n")
    cat("Percentage retained:", sprintf("%.2f%%", 100 * nrow(pass_variants) / nrow(info)), "\\n")
    """
}

// ============================================================================
// PROCESS 2: Calculate sample and genotype call rates
// ============================================================================
process calculateCallRates {
    label 'bcftools'
    tag "${platform}_chr${chr}"
    publishDir "${params.outdir}/module3/02_call_rates/${platform}", mode: 'copy'
    
    input:
    tuple val(platform), val(chr), path(vcf), path(tbi), path(stats)
    
    output:
    tuple val(platform), val(chr), path(vcf), path(tbi),
          path("${platform}_chr${chr}_sample_call_rate.txt"),
          path("${platform}_chr${chr}_variant_call_rate.txt"),
          emit: with_call_rates
    
    script:
    """
    echo "Calculating call rates for ${platform} chromosome ${chr}..."
    
    # Calculate per-sample call rate
    # Call rate = (# of called genotypes) / (# of total genotypes)
    bcftools query -f '[%SAMPLE\\t%GT\\n]' ${vcf} | \\
        awk '{
            sample=\$1; gt=\$2
            total[sample]++
            if (gt != "./." && gt != ".|.") called[sample]++
        } END {
            for (s in total) {
                call_rate = (called[s] > 0) ? called[s]/total[s] : 0
                print s, call_rate
            }
        }' > ${platform}_chr${chr}_sample_call_rate.txt
    
    # Calculate per-variant call rate
    # F_MISSING is fraction missing, so call rate = 1 - F_MISSING
    bcftools query -f '%ID\\t%F_MISSING\\n' ${vcf} | \\
        awk '{print \$1, 1-\$2}' > ${platform}_chr${chr}_variant_call_rate.txt
    
    echo "Call rate calculation complete"
    echo "Sample call rates written to: ${platform}_chr${chr}_sample_call_rate.txt"
    echo "Variant call rates written to: ${platform}_chr${chr}_variant_call_rate.txt"
    """
}

// ============================================================================
// PROCESS 3: Filter samples with <95% call rate
// ============================================================================
process filterLowCallRateSamples {
    label 'bcftools'
    tag "${platform}"
    publishDir "${params.outdir}/module3/03_sample_filtered/${platform}", mode: 'copy'
    
    input:
    tuple val(platform), path(vcfs), path(tbis), path(sample_call_rates)
    
    output:
    tuple val(platform),
          path("${platform}_chr*_sample_qc.vcf.gz"),
          path("${platform}_chr*_sample_qc.vcf.gz.tbi"),
          path("${platform}_removed_samples.txt"),
          path("${platform}_sample_qc_stats.txt"),
          emit: sample_filtered
    
    script:
    """
    echo "Filtering low call rate samples for ${platform}..."
    
    # Merge all chromosome call rates to get genome-wide call rate per sample
    # Average call rate across all chromosomes
    cat ${sample_call_rates} | \\
        awk '{sample[\$1]+=\$2; count[\$1]++} 
             END {for (s in sample) print s, sample[s]/count[s]}' | \\
        sort -k2,2n > ${platform}_genome_wide_call_rates.txt
    
    # Identify samples with <95% call rate
    awk '\$2 < 0.95 {print \$1}' ${platform}_genome_wide_call_rates.txt > ${platform}_removed_samples.txt
    
    # Count samples
    total_samples=\$(wc -l < ${platform}_genome_wide_call_rates.txt)
    removed_samples=\$(wc -l < ${platform}_removed_samples.txt)
    retained_samples=\$((total_samples - removed_samples))
    
    # Create stats file
    cat > ${platform}_sample_qc_stats.txt <<EOF
Metric\tValue
Total samples\t\${total_samples}
Samples removed (<95% call rate)\t\${removed_samples}
Samples retained\t\${retained_samples}
Percentage retained\t\$(awk "BEGIN {printf \\"%.2f\\", 100*\${retained_samples}/\${total_samples}}")%
EOF
    
    echo "Sample QC Statistics:"
    cat ${platform}_sample_qc_stats.txt
    
    # Filter each chromosome VCF to remove low call rate samples
    for vcf in ${vcfs}; do
        chr=\$(basename \$vcf | sed 's/.*chr\\([0-9XY]*\\).*/\\1/')
        out="${platform}_chr\${chr}_sample_qc.vcf.gz"
        
        echo "Processing chromosome \${chr}..."
        
        if [ -s ${platform}_removed_samples.txt ]; then
            # Remove samples listed in file
            bcftools view -S ^${platform}_removed_samples.txt \$vcf -Oz -o \$out
        else
            # No samples to remove, just copy
            echo "No samples to remove for chr\${chr}"
            cp \$vcf \$out
        fi
        
        # Index the filtered VCF
        bcftools index -t \$out
    done
    
    echo "\\n=== Sample Filtering Complete ==="
    echo "Removed \${removed_samples} samples with <95% call rate"
    echo "Retained \${retained_samples} samples"
    """
}

// ============================================================================
// PROCESS 4: Filter variants with <95% call rate
// ============================================================================
process filterLowCallRateVariants {
    label 'bcftools'
    tag "${platform}_chr${chr}"
    publishDir "${params.outdir}/module3/04_variant_filtered/${platform}", mode: 'copy'
    
    input:
    tuple val(platform), val(chr), path(vcf), path(tbi), path(variant_call_rate)
    
    output:
    tuple val(platform), val(chr),
          path("${platform}_chr${chr}_qc.vcf.gz"),
          path("${platform}_chr${chr}_qc.vcf.gz.tbi"),
          path("${platform}_chr${chr}_removed_variants.txt"),
          path("${platform}_chr${chr}_variant_qc_stats.txt"),
          emit: variant_filtered
    
    script:
    """
    echo "Filtering low call rate variants for ${platform} chromosome ${chr}..."
    
    # Count total variants
    total_variants=\$(wc -l < ${variant_call_rate})
    
    # Identify variants with <95% call rate
    awk '\$2 < 0.95 {print \$1}' ${variant_call_rate} > ${platform}_chr${chr}_removed_variants.txt
    
    removed_variants=\$(wc -l < ${platform}_chr${chr}_removed_variants.txt)
    retained_variants=\$((total_variants - removed_variants))
    
    # Create stats file
    cat > ${platform}_chr${chr}_variant_qc_stats.txt <<EOF
Metric\tValue
Chromosome\t${chr}
Total variants\t\${total_variants}
Variants removed (<95% call rate)\t\${removed_variants}
Variants retained\t\${retained_variants}
Percentage retained\t\$(awk "BEGIN {printf \\"%.2f\\", 100*\${retained_variants}/\${total_variants}}")%
EOF
    
    echo "Variant QC Statistics for chr${chr}:"
    cat ${platform}_chr${chr}_variant_qc_stats.txt
    
    # Filter VCF
    if [ -s ${platform}_chr${chr}_removed_variants.txt ]; then
        echo "Removing \${removed_variants} variants from chr${chr}..."
        bcftools view -e 'ID=@${platform}_chr${chr}_removed_variants.txt' \\
            ${vcf} -Oz -o ${platform}_chr${chr}_qc.vcf.gz
    else
        echo "No variants to remove for chr${chr}"
        cp ${vcf} ${platform}_chr${chr}_qc.vcf.gz
    fi
    
    # Index filtered VCF
    bcftools index -t ${platform}_chr${chr}_qc.vcf.gz
    
    echo "\\n=== Variant Filtering Complete for chr${chr} ==="
    echo "Removed \${removed_variants} variants with <95% call rate"
    echo "Retained \${retained_variants} variants"
    """
}

// ============================================================================
// PROCESS 5: Generate comprehensive QC summary report
// ============================================================================
process generateQCSummary {
    label 'R'
    tag "${platform}"
    publishDir "${params.outdir}/module3/05_qc_reports/${platform}", mode: 'copy'
    
    input:
    tuple val(platform), path(rsq_stats), path(sample_stats), path(variant_stats)
    
    output:
    tuple val(platform), 
          path("${platform}_module3_qc_report.html"),
          path("${platform}_module3_qc_summary.txt"),
          emit: qc_report
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(data.table)
    library(ggplot2)
    library(knitr)
    
    cat("\\n=== Generating QC Summary Report for ${platform} ===\\n")
    
    # ========== Read all statistics files ==========
    
    # 1. MagicalRsq filtering statistics
    rsq_files <- list.files(pattern="_rsq_stats.txt", full.names=TRUE)
    rsq_data <- rbindlist(lapply(rsq_files, fread))
    
    total_variants_input <- sum(rsq_data\$total_variants)
    total_variants_pass_rsq <- sum(rsq_data\$pass_variants)
    total_variants_filtered_rsq <- sum(rsq_data\$filtered_variants)
    pct_pass_rsq <- 100 * total_variants_pass_rsq / total_variants_input
    
    # 2. Sample filtering statistics
    sample_stats <- fread("${sample_stats}", header=TRUE)
    
    # 3. Variant filtering statistics
    variant_files <- list.files(pattern="_variant_qc_stats.txt", full.names=TRUE)
    variant_data <- rbindlist(lapply(variant_files, fread))
    
    total_variants_after_sample_filter <- sum(variant_data[, as.numeric(gsub(",", "", get("Total variants")))])
    total_variants_removed_call_rate <- sum(variant_data[, as.numeric(gsub(",", "", get("Variants removed (<95% call rate)")))])
    final_variant_count <- sum(variant_data[, as.numeric(gsub(",", "", get("Variants retained")))])
    
    # ========== Create comprehensive summary ==========
    
    summary <- data.frame(
        Step = c(
            "1. Initial imputed variants",
            "2. After MagicalRsq filtering (R² ≥ 0.3)",
            "3. After sample QC (call rate ≥ 95%)",
            "4. After variant QC (call rate ≥ 95%)",
            "5. Final QC'd variant count"
        ),
        Variants = c(
            format(total_variants_input, big.mark=","),
            format(total_variants_pass_rsq, big.mark=","),
            "—",
            "—",
            format(final_variant_count, big.mark=",")
        ),
        Removed = c(
            "—",
            format(total_variants_filtered_rsq, big.mark=","),
            "—",
            format(total_variants_removed_call_rate, big.mark=","),
            "—"
        ),
        Percentage = c(
            "100%",
            sprintf("%.2f%%", pct_pass_rsq),
            "—",
            sprintf("%.2f%%", 100 * final_variant_count / total_variants_pass_rsq),
            sprintf("%.2f%%", 100 * final_variant_count / total_variants_input)
        )
    )
    
    # Add sample statistics
    sample_summary <- data.frame(
        Metric = c(
            "Total samples",
            "Samples removed (<95% call rate)",
            "Samples retained",
            "Sample retention rate"
        ),
        Value = c(
            sample_stats[Metric == "Total samples", Value],
            sample_stats[Metric == "Samples removed (<95% call rate)", Value],
            sample_stats[Metric == "Samples retained", Value],
            sample_stats[Metric == "Percentage retained", Value]
        )
    )
    
    # ========== Write summary text file ==========
    
    sink("${platform}_module3_qc_summary.txt")
    cat("================================================================================\\n")
    cat("MODULE 3 POST-IMPUTATION QC SUMMARY\\n")
    cat("Platform: ${platform}\\n")
    cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")
    cat("================================================================================\\n\\n")
    
    cat("VARIANT FILTERING SUMMARY:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    print(kable(summary, align='lrrr'))
    cat("\\n")
    
    cat("SAMPLE FILTERING SUMMARY:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    print(kable(sample_summary, align='lr'))
    cat("\\n")
    
    cat("KEY METRICS:\\n")
    cat("  • Overall variant retention:", sprintf("%.2f%%", 100 * final_variant_count / total_variants_input), "\\n")
    cat("  • MagicalRsq effectiveness:", 
        sprintf("%.1fM variants filtered", total_variants_filtered_rsq/1e6), "\\n")
    cat("  • Sample retention:", sample_stats[Metric == "Percentage retained", Value], "\\n")
    cat("  • Final dataset:", format(final_variant_count, big.mark=","), 
        "variants across", sample_stats[Metric == "Samples retained", Value], "samples\\n")
    cat("\\n================================================================================\\n")
    sink()
    
    # ========== Generate HTML report ==========
    
    cat("\\nGenerating HTML report...\\n")
    
    # Create R Markdown content
    rmd_content <- sprintf('
---
title: "Module 3 Post-Imputation QC Report"
subtitle: "Platform: %s"
date: "`r format(Sys.time(), \\"%%Y-%%m-%%d %%H:%%M:%%S\\")`"
output: 
  html_document:
    theme: cosmo
    toc: true
    toc_float: true
    code_folding: hide
---

# Executive Summary

This report summarizes the quality control performed on imputed genotype data for **%s**.

## QC Pipeline Overview

The Module 3 QC pipeline consists of four sequential filtering steps:

1. **MagicalRsq Filtering**: Removes poorly imputed variants (empirical R² < 0.3)
2. **Sample QC**: Removes samples with low genotype call rates (<95%%)
3. **Variant QC**: Removes variants with low sample call rates (<95%%)
4. **Final QC Dataset**: High-quality, analysis-ready genotype data

---

# Variant Filtering Results

```{r variant_summary, echo=FALSE}
summary_table <- read.table("${platform}_module3_qc_summary.txt", 
                             skip=7, nrows=5, header=FALSE, sep="\\t",
                             col.names=c("Step", "Variants", "Removed", "Percentage"))
knitr::kable(summary_table, format="html", align="lrrr",
             caption="Variant Filtering Pipeline Summary")
```

## Key Findings

- **Initial variants**: %s
- **Final variants**: %s (%.2f%%%% retention)
- **MagicalRsq filtered**: %s variants
- **Call rate filtered**: %s variants

---

# Sample Filtering Results

```{r sample_summary, echo=FALSE}
sample_table <- read.table("${platform}_module3_qc_summary.txt",
                           skip=16, nrows=4, header=FALSE, sep="\\t",
                           col.names=c("Metric", "Value"))
knitr::kable(sample_table, format="html", align="lr",
             caption="Sample Quality Control Summary")
```

---

# Per-Chromosome Statistics

```{r per_chr_stats, echo=FALSE}
rsq_data <- read.table(textConnection("%s"), header=TRUE, sep="\\t")
knitr::kable(rsq_data[, c("chromosome", "total_variants", "pass_variants", 
                          "filtered_variants", "pct_pass")],
             format="html", digits=2,
             col.names=c("Chr", "Total", "Passed", "Filtered", "%%% Passed"),
             caption="MagicalRsq Filtering by Chromosome")
```

---

# Recommendations

## ✓ Quality Indicators
- Sample retention rate: **%s**
- Variant retention rate: **%.2f%%**
- Mean MagicalRsq for passing variants: High quality expected

## → Next Steps
1. **Module 4**: Merge platforms (if multiple)
2. **Module 5**: Re-imputation of merged data
3. **Module 6**: Final QC and relatedness filtering
4. **Module 7**: Ancestry estimation

---

# Technical Details

## MagicalRsq-X Features
- LD scores from TOP-LD (4 populations)
  - Long-range: ±1 Mb window
  - Short-range: ±100 kb window
- Recombination rates (1000 Genomes)
- Ancestry-specific MAF
- S/HIC features (selection/haplotype)

**Note**: MagicalRsq-X automatically incorporates ancestry-aware LD patterns without requiring ancestry input.

## Filtering Thresholds
- **MagicalRsq**: ≥ 0.3 (empirical R²)
- **Sample call rate**: ≥ 95%%
- **Variant call rate**: ≥ 95%%

---

*Report generated: `r Sys.time()`*
', 
    "${platform}", "${platform}",
    format(total_variants_input, big.mark=","),
    format(final_variant_count, big.mark=","),
    100 * final_variant_count / total_variants_input,
    format(total_variants_filtered_rsq, big.mark=","),
    format(total_variants_removed_call_rate, big.mark=","),
    paste(capture.output(write.table(rsq_data, sep="\\t", quote=FALSE, row.names=FALSE)), collapse="\\n"),
    sample_stats[Metric == "Percentage retained", Value],
    100 * final_variant_count / total_variants_input
    )
    
    writeLines(rmd_content, "${platform}_report.Rmd")
    
    # Render to HTML
    rmarkdown::render("${platform}_report.Rmd", 
                      output_file="${platform}_module3_qc_report.html",
                      quiet=TRUE)
    
    cat("\\n=== QC Report Generation Complete ===\\n")
    cat("Summary text file:", "${platform}_module3_qc_summary.txt", "\\n")
    cat("HTML report:", "${platform}_module3_qc_report.html", "\\n")
    """
}

// ============================================================================
// MODULE 3 WORKFLOW
// ============================================================================
workflow MODULE3_POSTQC {
    take:
    imputed_data  // From Module 2: tuple(platform, vcf_files, info_files)
    
    main:
    
    // ========== STEP 0: Index VCF files ==========
    // Flatten to per-chromosome channels for parallelization
    imputed_data
        .flatMap { platform, vcfs, infos ->
            def vcf_list = vcfs instanceof List ? vcfs : [vcfs].flatten()
            def info_list = infos instanceof List ? infos : [infos].flatten()
            
            vcf_list.collect { vcf ->
                def chr = (vcf.name =~ /chr(\d+|X|Y)/)[0][1]
                def info = info_list.find { it.name.contains("chr\${chr}") }
                tuple(platform, chr, vcf, info)
            }
        }
        .set { per_chr_data }
    
    // Create index files (.tbi) for all VCFs
    indexImputedVCFs(per_chr_data)
    
    // ========== STEP 1: MagicalRsq filtering ==========
    magicalRsqFilter(indexImputedVCFs.out.indexed_vcfs)
    
    // ========== STEP 2: Calculate call rates ==========
    calculateCallRates(magicalRsqFilter.out.filtered_vcf)
    
    // ========== STEP 3: Filter low call rate samples ==========
    // Group by platform for genome-wide sample filtering
    calculateCallRates.out.with_call_rates
        .map { platform, chr, vcf, tbi, sample_cr, variant_cr ->
            tuple(platform, vcf, tbi, sample_cr)
        }
        .groupTuple()
        .set { grouped_for_sample_filter }
    
    filterLowCallRateSamples(grouped_for_sample_filter)
    
    // ========== STEP 4: Filter low call rate variants ==========
    // Per-chromosome variant filtering
    calculateCallRates.out.with_call_rates
        .map { platform, chr, vcf, tbi, sample_cr, variant_cr ->
            tuple(platform, chr, vcf, tbi, variant_cr)
        }
        .set { for_variant_filter }
    
    filterLowCallRateVariants(for_variant_filter)
    
    // ========== STEP 5: Generate QC summary report ==========
    // Collect all QC statistics for report generation
    magicalRsqFilter.out.filtered_vcf
        .map { platform, chr, vcf, tbi, stats -> tuple(platform, stats) }
        .groupTuple()
        .join(
            filterLowCallRateSamples.out.sample_filtered
                .map { platform, vcfs, tbis, removed, stats -> tuple(platform, stats) }
        )
        .join(
            filterLowCallRateVariants.out.variant_filtered
                .map { platform, chr, vcf, tbi, removed, stats -> tuple(platform, stats) }
                .groupTuple()
        )
        .set { for_report }
    
    generateQCSummary(for_report)
    
    // ========== Prepare output channels ==========
    // Group final QC'd data by platform for Module 4
    filterLowCallRateVariants.out.variant_filtered
        .map { platform, chr, vcf, tbi, removed, stats ->
            tuple(platform, vcf, tbi)
        }
        .groupTuple()
        .set { qc_complete_data }
    
    emit:
    qc_data = qc_complete_data  // For Module 4: tuple(platform, vcfs[], tbis[])
    qc_reports = generateQCSummary.out.qc_report  // QC reports
}
