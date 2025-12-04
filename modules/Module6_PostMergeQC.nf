
/*
================================================================================
    MODULE 6: POST-MERGE QC (Second Pass)
================================================================================
    Comprehensive QC after re-imputation with:
    - Flexible input (from Module 5 OR Module 4)
    - MagicalRsq-X filtering (R² > 0.3 OR genotyped)
    - Call rate filtering (95% variant and sample)
    - HWE filtering (p > 1e-6)
    - Extreme heterozygosity removal (±3 SD)
    - Relatedness filtering (1st-degree: kinship > 0.177)
    - Separate TOPMed vs All of Us tracks
    - Two output branches: MAF > 0.01 and no MAF filter
    - Top 15 PCs calculation
    - Imputation performance benchmarking
================================================================================
*/

// ============================================================================
// PROCESS DEFINITIONS
// ============================================================================

/*
 * Process 1: Determine Input Source
 * Checks whether input is from Module 5 (re-imputed) or Module 4 (merged only)
 */
process CHECK_INPUT_SOURCE {
    tag "${sample_id}"
    label 'process_low'
    
    input:
    tuple val(sample_id), val(source), path(vcf), path(index), val(server)
    
    output:
    tuple val(sample_id), val(source), path(vcf), path(index), val(server), emit: validated
    
    script:
    """
    echo "Processing ${sample_id} from ${source} (Server: ${server})"
    
    # Validate VCF file
    bcftools query -l ${vcf} | head -5
    
    # Check if this is post-reimputation (Module 5) or post-merge (Module 4)
    if [[ "${source}" == "module5" ]]; then
        echo "Input from Module 5: Re-imputed merged data"
    elif [[ "${source}" == "module4" ]]; then
        echo "Input from Module 4: Merged platform data (skipped Module 5)"
    else
        echo "ERROR: Unknown source ${source}"
        exit 1
    fi
    """
}

/*
 * Process 2: Apply MagicalRsq-X Filtering
 * Uses pre-trained models to filter variants with R² > 0.3 OR genotyped variants
 * Separate processing for TOPMed vs All of Us
 */
process MAGICALRSQX_FILTER {
    tag "${sample_id}_chr${chr}_${server}"
    label 'process_medium'
    publishDir "${params.outdir}/module6/01_magicalrsq/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(source), path(vcf), path(index), val(server), val(chr)
    path(magicalrsqx_script)
    path(magicalrsqx_models)
    
    output:
    tuple val(sample_id), path("${sample_id}_chr${chr}_${server}_filtered.vcf.gz"), 
          path("${sample_id}_chr${chr}_${server}_filtered.vcf.gz.csi"), 
          val(server), val(chr), emit: filtered_vcf
    path("${sample_id}_chr${chr}_${server}_magicalrsq_stats.txt"), emit: stats
    path("${sample_id}_chr${chr}_${server}_rsq_comparison.pdf"), emit: plot
    
    script:
    def ancestry_flag = params.primary_ancestry == 'EUR' ? '--ancestry EUR' : '--ancestry AFR'
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    echo "=== MagicalRsq-X Filtering for ${sample_id} chr${chr} (${server}) ==="
    
    # Step 1: Calculate MagicalRsq-X scores using pre-trained models
    Rscript ${magicalrsqx_script} \\
        --vcf ${vcf} \\
        --models ${magicalrsqx_models} \\
        ${ancestry_flag} \\
        --output ${sample_id}_chr${chr}_${server}_magicalrsqx.txt \\
        --threads ${task.cpus}
    
    # Step 2: Apply filtering criteria
    # Keep variants if:
    # - MagicalRsq-X >= 0.3 (well-imputed)
    # - OR genotyped (no imputation quality score)
    
    bcftools view -h ${vcf} > header.txt
    
    bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/RSQ\\n' ${vcf} > variants.txt
    
    # Merge with MagicalRsq-X scores and filter
    Rscript -e "
    library(data.table)
    
    # Read variant info
    vars <- fread('variants.txt', 
                  col.names = c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'RSQ'))
    
    # Read MagicalRsq-X scores
    magic <- fread('${sample_id}_chr${chr}_${server}_magicalrsqx.txt')
    
    # Merge
    merged <- merge(vars, magic, by = c('CHROM', 'POS', 'ID'), all.x = TRUE)
    
    # Filter logic:
    # Keep if MagicalRsq-X >= 0.3 OR is genotyped (RSQ is NA or RSQ == 1)
    merged[, keep := (MagicalRsqX >= 0.3 | is.na(RSQ) | RSQ >= 0.99)]
    
    # Save filtering stats
    cat(sprintf('Total variants: %d\\n', nrow(merged)), 
        file = '${sample_id}_chr${chr}_${server}_magicalrsq_stats.txt')
    cat(sprintf('Genotyped variants: %d\\n', sum(is.na(merged\$RSQ) | merged\$RSQ >= 0.99)), 
        file = '${sample_id}_chr${chr}_${server}_magicalrsq_stats.txt', append = TRUE)
    cat(sprintf('Well-imputed (MagicalRsq-X >= 0.3): %d\\n', 
                sum(!is.na(merged\$MagicalRsqX) & merged\$MagicalRsqX >= 0.3, na.rm=TRUE)),
        file = '${sample_id}_chr${chr}_${server}_magicalrsq_stats.txt', append = TRUE)
    cat(sprintf('Poorly-imputed (filtered out): %d\\n', sum(!merged\$keep)), 
        file = '${sample_id}_chr${chr}_${server}_magicalrsq_stats.txt', append = TRUE)
    cat(sprintf('Variants passing filter: %d\\n', sum(merged\$keep)), 
        file = '${sample_id}_chr${chr}_${server}_magicalrsq_stats.txt', append = TRUE)
    
    # Create comparison plot
    library(ggplot2)
    pdf('${sample_id}_chr${chr}_${server}_rsq_comparison.pdf', width = 10, height = 6)
    
    # Plot 1: RSQ vs MagicalRsq-X
    p1 <- ggplot(merged[!is.na(RSQ) & !is.na(MagicalRsqX)], 
                 aes(x = RSQ, y = MagicalRsqX)) +
        geom_point(alpha = 0.1, size = 0.5) +
        geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 'dashed') +
        geom_hline(yintercept = 0.3, color = 'blue', linetype = 'dashed') +
        theme_bw() +
        labs(title = 'MagicalRsq-X vs Original Rsq',
             subtitle = paste('${sample_id} chr${chr} (${server})'),
             x = 'Original Rsq', y = 'MagicalRsq-X')
    print(p1)
    
    # Plot 2: Distribution by keep status
    p2 <- ggplot(merged, aes(x = MagicalRsqX, fill = keep)) +
        geom_histogram(bins = 50, alpha = 0.7) +
        geom_vline(xintercept = 0.3, color = 'red', linetype = 'dashed') +
        theme_bw() +
        labs(title = 'MagicalRsq-X Distribution',
             subtitle = 'Blue dashed line = 0.3 threshold',
             x = 'MagicalRsq-X', y = 'Count')
    print(p2)
    
    dev.off()
    
    # Write variant IDs to keep
    write.table(merged[keep == TRUE, .(ID)], 
                file = 'variants_to_keep.txt', 
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    " 
    
    # Step 3: Filter VCF
    bcftools view -i 'ID=@variants_to_keep.txt' ${vcf} \\
        -Oz -o ${sample_id}_chr${chr}_${server}_filtered.vcf.gz \\
        --threads ${task.cpus}
    
    bcftools index ${sample_id}_chr${chr}_${server}_filtered.vcf.gz
    
    echo "=== MagicalRsq-X filtering complete ==="
    """
}

/*
 * Process 3: Merge Chromosomes
 * Concatenate filtered chromosomes back together
 */
process CONCAT_CHROMOSOMES {
    tag "${sample_id}_${server}"
    label 'process_medium'
    publishDir "${params.outdir}/module6/02_concat/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), path(vcfs), path(indices)
    
    output:
    tuple val(sample_id), val(server), 
          path("${sample_id}_${server}_allchr.vcf.gz"),
          path("${sample_id}_${server}_allchr.vcf.gz.csi"), emit: concat_vcf
    
    script:
    """
    echo "=== Concatenating chromosomes for ${sample_id} (${server}) ==="
    
    # Create file list
    ls ${vcfs} | sort -V > vcf_list.txt
    
    # Concatenate
    bcftools concat -f vcf_list.txt \\
        -Oz -o ${sample_id}_${server}_allchr.vcf.gz \\
        --threads ${task.cpus}
    
    bcftools index ${sample_id}_${server}_allchr.vcf.gz
    
    echo "=== Concatenation complete ==="
    """
}

/*
 * Process 4: Convert to PLINK format for QC
 */
process VCF_TO_PLINK {
    tag "${sample_id}_${server}"
    label 'process_medium'
    publishDir "${params.outdir}/module6/03_plink/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), path(vcf), path(index)
    
    output:
    tuple val(sample_id), val(server),
          path("${sample_id}_${server}.bed"),
          path("${sample_id}_${server}.bim"),
          path("${sample_id}_${server}.fam"), emit: plink_files
    
    script:
    """
    echo "=== Converting VCF to PLINK for ${sample_id} (${server}) ==="
    
    plink2 --vcf ${vcf} dosage=DS \\
        --make-bed \\
        --out ${sample_id}_${server} \\
        --threads ${task.cpus} \\
        --memory ${task.memory.toMega()}
    
    echo "=== Conversion complete ==="
    """
}

/*
 * Process 5: Calculate Sample and Variant Call Rates
 */
process CALCULATE_CALL_RATES {
    tag "${sample_id}_${server}"
    label 'process_low'
    publishDir "${params.outdir}/module6/04_callrates/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), path(bed), path(bim), path(fam)
    
    output:
    tuple val(sample_id), val(server), path(bed), path(bim), path(fam),
          path("${sample_id}_${server}.imiss"),
          path("${sample_id}_${server}.lmiss"), emit: with_missingness
    path("${sample_id}_${server}_callrate_summary.txt"), emit: summary
    
    script:
    """
    echo "=== Calculating call rates for ${sample_id} (${server}) ==="
    
    # Calculate missingness
    plink --bfile ${sample_id}_${server} \\
        --missing \\
        --out ${sample_id}_${server} \\
        --threads ${task.cpus}
    
    # Summarize
    echo "Sample Call Rate Summary:" > ${sample_id}_${server}_callrate_summary.txt
    awk 'NR>1 {print \$6}' ${sample_id}_${server}.imiss | \\
        awk '{sum+=\$1; sumsq+=\$1*\$1} END {
            print "Mean sample missingness: " sum/NR;
            print "SD sample missingness: " sqrt(sumsq/NR - (sum/NR)^2);
            print "Min sample missingness: " (NR>0 ? min : 0);
            print "Max sample missingness: " max
        }' BEGIN{min=1;max=0} {if(\$1<min)min=\$1; if(\$1>max)max=\$1} >> ${sample_id}_${server}_callrate_summary.txt
    
    echo "" >> ${sample_id}_${server}_callrate_summary.txt
    echo "Variant Call Rate Summary:" >> ${sample_id}_${server}_callrate_summary.txt
    awk 'NR>1 {print \$5}' ${sample_id}_${server}.lmiss | \\
        awk '{sum+=\$1; sumsq+=\$1*\$1} END {
            print "Mean variant missingness: " sum/NR;
            print "SD variant missingness: " sqrt(sumsq/NR - (sum/NR)^2)
        }' >> ${sample_id}_${server}_callrate_summary.txt
    
    echo "=== Call rate calculation complete ==="
    """
}

/*
 * Process 6: Filter by Call Rate (95%)
 * Remove samples with <95% call rate
 * Remove variants with <95% call rate
 */
process FILTER_CALL_RATE {
    tag "${sample_id}_${server}"
    label 'process_medium'
    publishDir "${params.outdir}/module6/05_callrate_filtered/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), path(bed), path(bim), path(fam),
          path(imiss), path(lmiss)
    
    output:
    tuple val(sample_id), val(server),
          path("${sample_id}_${server}_cr95.bed"),
          path("${sample_id}_${server}_cr95.bim"),
          path("${sample_id}_${server}_cr95.fam"), emit: filtered_files
    path("${sample_id}_${server}_callrate_filter_log.txt"), emit: filter_log
    
    script:
    def mind = params.sample_call_rate ?: 0.05  // 95% call rate = 5% missingness
    def geno = params.variant_call_rate ?: 0.05
    """
    echo "=== Filtering by call rate (${server}) ==="
    echo "Sample threshold: ${mind * 100}% missingness (${(1-mind)*100}% call rate)"
    echo "Variant threshold: ${geno * 100}% missingness (${(1-geno)*100}% call rate)"
    
    # Count before filtering
    n_samples_before=\$(wc -l < ${fam})
    n_variants_before=\$(wc -l < ${bim})
    
    # Filter
    plink --bfile ${sample_id}_${server} \\
        --mind ${mind} \\
        --geno ${geno} \\
        --make-bed \\
        --out ${sample_id}_${server}_cr95 \\
        --threads ${task.cpus}
    
    # Count after filtering
    n_samples_after=\$(wc -l < ${sample_id}_${server}_cr95.fam)
    n_variants_after=\$(wc -l < ${sample_id}_${server}_cr95.bim)
    
    # Log results
    {
        echo "Call Rate Filtering Results (${server})"
        echo "========================================"
        echo "Samples before: \$n_samples_before"
        echo "Samples after: \$n_samples_after"
        echo "Samples removed: \$((n_samples_before - n_samples_after))"
        echo ""
        echo "Variants before: \$n_variants_before"
        echo "Variants after: \$n_variants_after"
        echo "Variants removed: \$((n_variants_before - n_variants_after))"
    } > ${sample_id}_${server}_callrate_filter_log.txt
    
    echo "=== Call rate filtering complete ==="
    """
}

/*
 * Process 7: Filter by Hardy-Weinberg Equilibrium (HWE)
 * Remove variants with HWE p-value < threshold
 * NOTE: Skipped by default (params.skip_hwe=true) for admixed populations
 *       HWE assumptions are violated in admixed/structured populations
 *       Enable with --skip_hwe false if working with homogeneous populations
 */
process FILTER_HWE {
    tag "${sample_id}_${server}"
    label 'process_medium'
    publishDir "${params.outdir}/module6/06_hwe_filtered/${server}", mode: 'copy'

    input:
    tuple val(sample_id), val(server), path(bed), path(bim), path(fam)

    when:
    !params.skip_hwe

    output:
    tuple val(sample_id), val(server),
          path("${sample_id}_${server}_hwe.bed"),
          path("${sample_id}_${server}_hwe.bim"),
          path("${sample_id}_${server}_hwe.fam"), emit: filtered_files
    path("${sample_id}_${server}.hwe"), emit: hwe_stats
    path("${sample_id}_${server}_hwe_filter_log.txt"), emit: filter_log

    script:
    def hwe_threshold = params.hwe_pvalue ?: 1e-6
    """
    echo "=== HWE filtering for ${sample_id} (${server}) ==="
    echo "HWE threshold: ${hwe_threshold}"

    # Count before
    n_variants_before=\$(wc -l < ${bim})

    # Calculate HWE
    plink --bfile ${sample_id}_${server} \\
        --hardy \\
        --out ${sample_id}_${server} \\
        --threads ${task.cpus}

    # Filter variants failing HWE
    plink --bfile ${sample_id}_${server} \\
        --hwe ${hwe_threshold} \\
        --make-bed \\
        --out ${sample_id}_${server}_hwe \\
        --threads ${task.cpus}

    # Count after
    n_variants_after=\$(wc -l < ${sample_id}_${server}_hwe.bim)

    # Log
    {
        echo "HWE Filtering Results (${server})"
        echo "========================================"
        echo "HWE threshold: ${hwe_threshold}"
        echo "Variants before: \$n_variants_before"
        echo "Variants after: \$n_variants_after"
        echo "Variants removed: \$((n_variants_before - n_variants_after))"
    } > ${sample_id}_${server}_hwe_filter_log.txt

    echo "=== HWE filtering complete ==="
    """
}

/*
 * Process 8: Calculate and Filter Extreme Heterozygosity
 * Remove samples with heterozygosity >3 SD from mean
 */
process FILTER_HETEROZYGOSITY {
    tag "${sample_id}_${server}"
    label 'process_medium'
    publishDir "${params.outdir}/module6/07_het_filtered/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), path(bed), path(bim), path(fam)
    
    output:
    tuple val(sample_id), val(server),
          path("${sample_id}_${server}_het.bed"),
          path("${sample_id}_${server}_het.bim"),
          path("${sample_id}_${server}_het.fam"), emit: filtered_files
    path("${sample_id}_${server}.het"), emit: het_stats
    path("${sample_id}_${server}_het_outliers.txt"), emit: outliers
    path("${sample_id}_${server}_het_plot.pdf"), emit: het_plot
    
    script:
    def sd_threshold = params.het_sd_threshold ?: 3
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    echo "=== Heterozygosity filtering for ${sample_id} (${server}) ==="
    
    # Calculate heterozygosity
    plink --bfile ${sample_id}_${server} \\
        --het \\
        --out ${sample_id}_${server} \\
        --threads ${task.cpus}
    
    # Identify outliers in R
    Rscript - <<'EOF'
    library(data.table)
    library(ggplot2)
    
    # Read heterozygosity data
    het <- fread("${sample_id}_${server}.het")
    
    # Calculate heterozygosity rate: (N(NM) - O(HOM)) / N(NM)
    het[, het_rate := (N.NM. - O.HOM.) / N.NM.]
    
    # Calculate mean and SD
    mean_het <- mean(het\$het_rate, na.rm = TRUE)
    sd_het <- sd(het\$het_rate, na.rm = TRUE)
    
    # Define outliers (>3 SD from mean)
    het[, outlier := abs(het_rate - mean_het) > ${sd_threshold} * sd_het]
    
    # Save outliers
    outliers <- het[outlier == TRUE, .(FID, IID)]
    fwrite(outliers, "${sample_id}_${server}_het_outliers.txt", 
           sep = "\\t", col.names = FALSE)
    
    # Summary stats
    cat(sprintf("Mean heterozygosity: %.4f\\n", mean_het))
    cat(sprintf("SD heterozygosity: %.4f\\n", sd_het))
    cat(sprintf("Lower bound (mean - %d SD): %.4f\\n", 
                ${sd_threshold}, mean_het - ${sd_threshold} * sd_het))
    cat(sprintf("Upper bound (mean + %d SD): %.4f\\n", 
                ${sd_threshold}, mean_het + ${sd_threshold} * sd_het))
    cat(sprintf("Number of outliers: %d\\n", sum(het\$outlier)))
    
    # Create plot
    pdf("${sample_id}_${server}_het_plot.pdf", width = 10, height = 6)
    
    p <- ggplot(het, aes(x = het_rate, fill = outlier)) +
        geom_histogram(bins = 50, alpha = 0.7) +
        geom_vline(xintercept = mean_het, color = "blue", linewidth = 1) +
        geom_vline(xintercept = c(mean_het - ${sd_threshold} * sd_het,
                                   mean_het + ${sd_threshold} * sd_het),
                   color = "red", linetype = "dashed", linewidth = 1) +
        theme_bw() +
        labs(title = "Heterozygosity Distribution",
             subtitle = paste("${sample_id} (${server})",
                            sprintf("Mean = %.4f, SD = %.4f", mean_het, sd_het)),
             x = "Heterozygosity Rate",
             y = "Count",
             fill = "Outlier") +
        scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "red"))
    
    print(p)
    dev.off()
EOF
    
    # Remove outliers if any exist
    if [ -s ${sample_id}_${server}_het_outliers.txt ]; then
        echo "Removing heterozygosity outliers..."
        plink --bfile ${sample_id}_${server} \\
            --remove ${sample_id}_${server}_het_outliers.txt \\
            --make-bed \\
            --out ${sample_id}_${server}_het \\
            --threads ${task.cpus}
    else
        echo "No heterozygosity outliers detected. Copying files..."
        cp ${sample_id}_${server}.bed ${sample_id}_${server}_het.bed
        cp ${sample_id}_${server}.bim ${sample_id}_${server}_het.bim
        cp ${sample_id}_${server}.fam ${sample_id}_${server}_het.fam
    fi
    
    echo "=== Heterozygosity filtering complete ==="
    """
}

/*
 * Process 9: Calculate KING Robust Kinship
 * Identify related individuals (1st-degree: kinship > 0.177)
 */
process CALCULATE_KINSHIP {
    tag "${sample_id}_${server}"
    label 'process_high'
    publishDir "${params.outdir}/module6/08_kinship/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), path(bed), path(bim), path(fam)
    
    output:
    tuple val(sample_id), val(server), path(bed), path(bim), path(fam),
          path("${sample_id}_${server}.kin0"), emit: with_kinship
    path("${sample_id}_${server}_kinship_summary.txt"), emit: summary
    
    script:
    """
    echo "=== Calculating kinship for ${sample_id} (${server}) ==="
    
    # Prune to independent SNPs for kinship calculation
    plink --bfile ${sample_id}_${server} \\
        --indep-pairwise 50 5 0.2 \\
        --out ${sample_id}_${server}_pruned \\
        --threads ${task.cpus}
    
    # Calculate KING robust kinship
    plink --bfile ${sample_id}_${server} \\
        --extract ${sample_id}_${server}_pruned.prune.in \\
        --make-king-table \\
        --king-table-filter ${params.kinship_threshold ?: 0.177} \\
        --out ${sample_id}_${server} \\
        --threads ${task.cpus}
    
    # Summarize related pairs
    if [ -f ${sample_id}_${server}.kin0 ]; then
        n_related_pairs=\$(tail -n +2 ${sample_id}_${server}.kin0 | wc -l)
        echo "Found \$n_related_pairs pairs of 1st-degree relatives" > ${sample_id}_${server}_kinship_summary.txt
        
        # Count unique individuals in related pairs
        tail -n +2 ${sample_id}_${server}.kin0 | \\
            awk '{print \$1"\\t"\$2"\\n"\$3"\\t"\$4}' | \\
            sort -u | wc -l | \\
            awk '{print "Unique individuals in related pairs: " \$1}' >> ${sample_id}_${server}_kinship_summary.txt
    else
        echo "No 1st-degree relatives found" > ${sample_id}_${server}_kinship_summary.txt
    fi
    
    echo "=== Kinship calculation complete ==="
    """
}

/*
 * Process 10: Remove Related Individuals
 * Prioritize samples with higher call rates
 */
process REMOVE_RELATEDS {
    tag "${sample_id}_${server}"
    label 'process_medium'
    publishDir "${params.outdir}/module6/09_unrelated/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), path(bed), path(bim), path(fam), path(kin0)
    path(imiss)  // From earlier call rate calculation
    
    output:
    tuple val(sample_id), val(server),
          path("${sample_id}_${server}_unrelated.bed"),
          path("${sample_id}_${server}_unrelated.bim"),
          path("${sample_id}_${server}_unrelated.fam"), emit: unrelated_files
    path("${sample_id}_${server}_removed_related.txt"), emit: removed_list
    path("${sample_id}_${server}_relatedness_log.txt"), emit: log
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    echo "=== Removing related individuals for ${sample_id} (${server}) ==="
    
    # Check if there are related pairs
    if [ ! -f ${kin0} ] || [ ! -s ${kin0} ]; then
        echo "No related individuals to remove"
        cp ${bed} ${sample_id}_${server}_unrelated.bed
        cp ${bim} ${sample_id}_${server}_unrelated.bim
        cp ${fam} ${sample_id}_${server}_unrelated.fam
        touch ${sample_id}_${server}_removed_related.txt
        echo "No related individuals found" > ${sample_id}_${server}_relatedness_log.txt
        exit 0
    fi
    
    # Prioritize samples with higher call rates
    Rscript - <<'EOF'
    library(data.table)
    
    # Read kinship file
    kin <- fread("${kin0}")
    setnames(kin, c("FID1", "IID1", "FID2", "IID2", "NSNP", "KINSHIP"))
    
    # Read missingness data
    miss <- fread("${imiss}")
    setnames(miss, c("FID", "IID", "MISS_PHENO", "N_MISS", "N_GENO", "F_MISS"))
    
    # For each related pair, remove the sample with higher missingness
    to_remove <- data.table()
    
    for (i in 1:nrow(kin)) {
        id1_miss <- miss[FID == kin\$FID1[i] & IID == kin\$IID1[i], F_MISS]
        id2_miss <- miss[FID == kin\$FID2[i] & IID == kin\$IID2[i], F_MISS]
        
        if (length(id1_miss) == 0 || length(id2_miss) == 0) next
        
        # Remove sample with higher missingness
        if (id1_miss > id2_miss) {
            to_remove <- rbind(to_remove, 
                             data.table(FID = kin\$FID1[i], IID = kin\$IID1[i]))
        } else {
            to_remove <- rbind(to_remove, 
                             data.table(FID = kin\$FID2[i], IID = kin\$IID2[i]))
        }
    }
    
    # Remove duplicates
    to_remove <- unique(to_remove)
    
    # Write removal list
    if (nrow(to_remove) > 0) {
        fwrite(to_remove, "${sample_id}_${server}_removed_related.txt",
               sep = "\\t", col.names = FALSE)
        cat(sprintf("Removing %d related individuals\\n", nrow(to_remove)))
    } else {
        file.create("${sample_id}_${server}_removed_related.txt")
        cat("No individuals to remove\\n")
    }
    
    # Log
    cat(sprintf("Total related pairs: %d\\n", nrow(kin)),
        file = "${sample_id}_${server}_relatedness_log.txt")
    cat(sprintf("Individuals removed: %d\\n", nrow(to_remove)),
        file = "${sample_id}_${server}_relatedness_log.txt", append = TRUE)
EOF
    
    # Remove related individuals
    if [ -s ${sample_id}_${server}_removed_related.txt ]; then
        plink --bfile ${sample_id}_${server} \\
            --remove ${sample_id}_${server}_removed_related.txt \\
            --make-bed \\
            --out ${sample_id}_${server}_unrelated \\
            --threads ${task.cpus}
    else
        cp ${bed} ${sample_id}_${server}_unrelated.bed
        cp ${bim} ${sample_id}_${server}_unrelated.bim
        cp ${fam} ${sample_id}_${server}_unrelated.fam
    fi
    
    echo "=== Relatedness filtering complete ==="
    """
}

/*
 * Process 11: Calculate Principal Components (Top 15 PCs)
 */
process CALCULATE_PCS {
    tag "${sample_id}_${server}"
    label 'process_high'
    publishDir "${params.outdir}/module6/10_pca/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), path(bed), path(bim), path(fam)
    
    output:
    tuple val(sample_id), val(server), 
          path("${sample_id}_${server}_pcs.eigenvec"),
          path("${sample_id}_${server}_pcs.eigenval"), emit: pcs
    path("${sample_id}_${server}_pc_plot.pdf"), emit: plot
    
    script:
    def n_pcs = params.n_pcs ?: 15
    """
    echo "=== Calculating ${n_pcs} PCs for ${sample_id} (${server}) ==="
    
    # LD pruning for PCA
    plink --bfile ${sample_id}_${server} \\
        --indep-pairwise 50 5 0.2 \\
        --out ${sample_id}_${server}_pca_pruned \\
        --threads ${task.cpus}
    
    # Calculate PCs
    plink2 --bfile ${sample_id}_${server} \\
        --extract ${sample_id}_${server}_pca_pruned.prune.in \\
        --pca ${n_pcs} \\
        --out ${sample_id}_${server}_pcs \\
        --threads ${task.cpus}
    
    # Create PCA plot
    Rscript - <<'EOF'
    library(data.table)
    library(ggplot2)
    
    # Read PCs
    pcs <- fread("${sample_id}_${server}_pcs.eigenvec")
    eigenval <- scan("${sample_id}_${server}_pcs.eigenval")
    
    # Calculate variance explained
    pve <- eigenval / sum(eigenval) * 100
    
    # Create plots
    pdf("${sample_id}_${server}_pc_plot.pdf", width = 12, height = 8)
    
    # Scree plot
    par(mfrow = c(2, 2))
    plot(1:length(pve), pve, type = "b", 
         xlab = "Principal Component", 
         ylab = "Variance Explained (%)",
         main = "Scree Plot")
    
    # PC1 vs PC2
    plot(pcs\$PC1, pcs\$PC2, 
         xlab = sprintf("PC1 (%.2f%%)", pve[1]),
         ylab = sprintf("PC2 (%.2f%%)", pve[2]),
         main = "PC1 vs PC2", pch = 16, cex = 0.5)
    
    # PC3 vs PC4
    plot(pcs\$PC3, pcs\$PC4,
         xlab = sprintf("PC3 (%.2f%%)", pve[3]),
         ylab = sprintf("PC4 (%.2f%%)", pve[4]),
         main = "PC3 vs PC4", pch = 16, cex = 0.5)
    
    # PC5 vs PC6
    plot(pcs\$PC5, pcs\$PC6,
         xlab = sprintf("PC5 (%.2f%%)", pve[5]),
         ylab = sprintf("PC6 (%.2f%%)", pve[6]),
         main = "PC5 vs PC6", pch = 16, cex = 0.5)
    
    dev.off()
EOF
    
    echo "=== PCA complete ==="
    """
}

/*
 * Process 12: MAF Filtering - Create Two Branches
 * Branch 1: MAF > 0.01 (for common variant GWAS)
 * Branch 2: No MAF filter (for rare variant studies)
 */
process MAF_FILTER_SPLIT {
    tag "${sample_id}_${server}"
    label 'process_medium'
    publishDir "${params.outdir}/module6/11_final_datasets/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), path(bed), path(bim), path(fam)
    
    output:
    tuple val(sample_id), val(server),
          path("${sample_id}_${server}_final_maf01.bed"),
          path("${sample_id}_${server}_final_maf01.bim"),
          path("${sample_id}_${server}_final_maf01.fam"), emit: maf_filtered
    tuple val(sample_id), val(server),
          path("${sample_id}_${server}_final_nomaf.bed"),
          path("${sample_id}_${server}_final_nomaf.bim"),
          path("${sample_id}_${server}_final_nomaf.fam"), emit: no_maf_filter
    path("${sample_id}_${server}_maf_summary.txt"), emit: summary
    
    script:
    def maf_threshold = params.maf_threshold ?: 0.01
    """
    echo "=== Creating MAF-filtered branches for ${sample_id} (${server}) ==="
    
    # Branch 1: MAF > 0.01
    echo "Creating MAF > ${maf_threshold} dataset..."
    plink --bfile ${sample_id}_${server} \\
        --maf ${maf_threshold} \\
        --make-bed \\
        --out ${sample_id}_${server}_final_maf01 \\
        --threads ${task.cpus}
    
    # Branch 2: No MAF filter (copy files)
    echo "Creating no-MAF-filter dataset..."
    cp ${bed} ${sample_id}_${server}_final_nomaf.bed
    cp ${bim} ${sample_id}_${server}_final_nomaf.bim
    cp ${fam} ${sample_id}_${server}_final_nomaf.fam
    
    # Summarize
    n_variants_nomaf=\$(wc -l < ${sample_id}_${server}_final_nomaf.bim)
    n_variants_maf=\$(wc -l < ${sample_id}_${server}_final_maf01.bim)
    n_samples=\$(wc -l < ${sample_id}_${server}_final_nomaf.fam)
    
    {
        echo "MAF Filtering Summary (${server})"
        echo "========================================"
        echo "Number of samples: \$n_samples"
        echo ""
        echo "No MAF filter dataset:"
        echo "  Total variants: \$n_variants_nomaf"
        echo ""
        echo "MAF > ${maf_threshold} dataset:"
        echo "  Total variants: \$n_variants_maf"
        echo "  Variants removed: \$((n_variants_nomaf - n_variants_maf))"
        echo "  Percent retained: \$(echo "scale=2; \$n_variants_maf * 100 / \$n_variants_nomaf" | bc)%"
    } > ${sample_id}_${server}_maf_summary.txt
    
    echo "=== MAF filtering complete ==="
    """
}

/*
 * Process 13: Generate QC Summary Report
 */
process GENERATE_QC_SUMMARY {
    tag "${sample_id}_${server}"
    label 'process_low'
    publishDir "${params.outdir}/module6/12_qc_reports/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), val(stage), path(report_files)
    
    output:
    path("${sample_id}_${server}_module6_qc_report.html"), emit: html_report
    path("${sample_id}_${server}_module6_qc_summary.txt"), emit: txt_summary
    
    script:
    """
    #!/usr/bin/env bash
    
    cat > ${sample_id}_${server}_module6_qc_summary.txt << EOF
================================================================================
MODULE 6: POST-MERGE QC SUMMARY
================================================================================
Sample ID: ${sample_id}
Imputation Server: ${server}
Date: \$(date)

QC Pipeline Steps:
--------------------------------------------------------------------------------
1. MagicalRsq-X Filtering (R² > 0.3 or genotyped)
2. Call Rate Filtering (95% for samples and variants)
3. Hardy-Weinberg Equilibrium (p > 1e-6)
4. Heterozygosity Outliers (±3 SD)
5. Relatedness Removal (1st-degree, kinship > 0.177)
6. Principal Components (Top 15 PCs)
7. MAF Filtering (Two branches: MAF > 0.01 and no filter)

See individual step reports in: ${params.outdir}/module6/

Final Datasets:
--------------------------------------------------------------------------------
- ${sample_id}_${server}_final_maf01.*  (Common variants, MAF > 0.01)
- ${sample_id}_${server}_final_nomaf.*  (All variants, no MAF filter)

================================================================================
EOF

    # Create simple HTML report
    cat > ${sample_id}_${server}_module6_qc_report.html << 'HTMLEOF'
<!DOCTYPE html>
<html>
<head>
    <title>Module 6 QC Report: ${sample_id} (${server})</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1 { color: #333; }
        h2 { color: #666; margin-top: 30px; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
        th { background-color: #4CAF50; color: white; }
        .summary-box { background-color: #f0f0f0; padding: 15px; margin: 20px 0; border-radius: 5px; }
    </style>
</head>
<body>
    <h1>Module 6: Post-Merge QC Report</h1>
    <div class="summary-box">
        <strong>Sample ID:</strong> ${sample_id}<br>
        <strong>Imputation Server:</strong> ${server}<br>
        <strong>Date:</strong> \$(date)
    </div>
    
    <h2>QC Pipeline Overview</h2>
    <ol>
        <li>MagicalRsq-X Filtering</li>
        <li>Call Rate Filtering</li>
        <li>Hardy-Weinberg Equilibrium</li>
        <li>Heterozygosity Filtering</li>
        <li>Relatedness Removal</li>
        <li>Principal Components Analysis</li>
        <li>MAF Filtering (Dual Output)</li>
    </ol>
    
    <h2>Final Datasets</h2>
    <p>Two analysis-ready datasets have been created:</p>
    <ul>
        <li><strong>MAF > 0.01:</strong> ${sample_id}_${server}_final_maf01.* (for common variant GWAS)</li>
        <li><strong>No MAF filter:</strong> ${sample_id}_${server}_final_nomaf.* (for rare variant studies)</li>
    </ul>
    
    <p>For detailed QC metrics, see the individual step reports in the output directory.</p>
</body>
</html>
HTMLEOF
    
    echo "=== QC summary report generated ==="
    """
}

/*
 * Process 14: Compare Imputation Performance
 * Compare TOPMed vs All of Us (if both available)
 */
process COMPARE_IMPUTATION_PERFORMANCE {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/module6/13_imputation_comparison", mode: 'copy'
    
    input:
    tuple val(sample_id), 
          path(topmed_bim), path(topmed_stats),
          path(allofus_bim), path(allofus_stats)
    
    output:
    path("${sample_id}_imputation_comparison.html"), emit: html_report
    path("${sample_id}_imputation_comparison.pdf"), emit: pdf_plots
    path("${sample_id}_imputation_metrics.txt"), emit: metrics
    
    when:
    topmed_bim.name != 'NO_FILE' && allofus_bim.name != 'NO_FILE'
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(data.table)
    library(ggplot2)
    library(cowplot)
    
    # Read variant counts
    topmed_vars <- fread("${topmed_bim}", header = FALSE)
    allofus_vars <- fread("${allofus_bim}", header = FALSE)
    
    # Read QC stats
    topmed_stats <- readLines("${topmed_stats}")
    allofus_stats <- readLines("${allofus_stats}")
    
    # Comparison metrics
    metrics <- list(
        topmed_n_variants = nrow(topmed_vars),
        allofus_n_variants = nrow(allofus_vars),
        shared_variants = length(intersect(topmed_vars\$V2, allofus_vars\$V2)),
        topmed_unique = length(setdiff(topmed_vars\$V2, allofus_vars\$V2)),
        allofus_unique = length(setdiff(allofus_vars\$V2, topmed_vars\$V2))
    )
    
    # Calculate MAF distributions
    topmed_maf <- as.numeric(sub(".*MAF=([0-9.]+).*", "\\\\1", topmed_vars\$V5))
    allofus_maf <- as.numeric(sub(".*MAF=([0-9.]+).*", "\\\\1", allofus_vars\$V5))
    
    # Create comparison plots
    pdf("${sample_id}_imputation_comparison.pdf", width = 14, height = 10)
    
    # Plot 1: Variant counts
    count_data <- data.frame(
        Server = c("TOPMed", "All of Us", "Shared", "TOPMed Unique", "All of Us Unique"),
        Count = c(metrics\$topmed_n_variants, metrics\$allofus_n_variants,
                 metrics\$shared_variants, metrics\$topmed_unique, metrics\$allofus_unique),
        Category = c("Total", "Total", "Shared", "Unique", "Unique")
    )
    
    p1 <- ggplot(count_data[1:2,], aes(x = Server, y = Count, fill = Server)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = scales::comma(Count)), vjust = -0.5) +
        theme_bw() +
        labs(title = "Total Variants by Imputation Server",
             y = "Number of Variants") +
        scale_fill_manual(values = c("TOPMed" = "#4CAF50", "All of Us" = "#2196F3"))
    
    # Plot 2: Venn diagram representation
    p2 <- ggplot(count_data[3:5,], aes(x = Server, y = Count, fill = Category)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = scales::comma(Count)), vjust = -0.5) +
        theme_bw() +
        labs(title = "Shared vs Unique Variants",
             y = "Number of Variants") +
        scale_fill_manual(values = c("Shared" = "#FFC107", "Unique" = "#E91E63"))
    
    # Plot 3: MAF distributions
    maf_data <- rbind(
        data.frame(MAF = topmed_maf[!is.na(topmed_maf)], Server = "TOPMed"),
        data.frame(MAF = allofus_maf[!is.na(allofus_maf)], Server = "All of Us")
    )
    
    p3 <- ggplot(maf_data, aes(x = MAF, fill = Server)) +
        geom_histogram(alpha = 0.5, position = "identity", bins = 50) +
        theme_bw() +
        labs(title = "MAF Distribution Comparison",
             x = "Minor Allele Frequency",
             y = "Count") +
        scale_fill_manual(values = c("TOPMed" = "#4CAF50", "All of Us" = "#2196F3"))
    
    # Combine plots
    combined <- plot_grid(p1, p2, p3, ncol = 2, nrow = 2)
    print(combined)
    
    dev.off()
    
    # Write metrics
    cat("Imputation Performance Comparison\\n", file = "${sample_id}_imputation_metrics.txt")
    cat("================================================================================\\n", 
        file = "${sample_id}_imputation_metrics.txt", append = TRUE)
    cat(sprintf("TOPMed variants: %s\\n", scales::comma(metrics\$topmed_n_variants)),
        file = "${sample_id}_imputation_metrics.txt", append = TRUE)
    cat(sprintf("All of Us variants: %s\\n", scales::comma(metrics\$allofus_n_variants)),
        file = "${sample_id}_imputation_metrics.txt", append = TRUE)
    cat(sprintf("Shared variants: %s (%.1f%%)\\n", 
                scales::comma(metrics\$shared_variants),
                100 * metrics\$shared_variants / min(metrics\$topmed_n_variants, metrics\$allofus_n_variants)),
        file = "${sample_id}_imputation_metrics.txt", append = TRUE)
    cat(sprintf("TOPMed unique: %s\\n", scales::comma(metrics\$topmed_unique)),
        file = "${sample_id}_imputation_metrics.txt", append = TRUE)
    cat(sprintf("All of Us unique: %s\\n", scales::comma(metrics\$allofus_unique)),
        file = "${sample_id}_imputation_metrics.txt", append = TRUE)
    
    # Create HTML report
    cat(paste0('
    <!DOCTYPE html>
    <html>
    <head>
        <title>Imputation Performance Comparison</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; }
            h1 { color: #333; }
            table { border-collapse: collapse; width: 80%; margin: 20px 0; }
            th, td { border: 1px solid #ddd; padding: 12px; text-align: right; }
            th { background-color: #4CAF50; color: white; text-align: left; }
            .metric { font-weight: bold; }
        </style>
    </head>
    <body>
        <h1>Imputation Performance Comparison</h1>
        <h2>Sample: ${sample_id}</h2>
        
        <h3>Variant Counts</h3>
        <table>
            <tr><th>Metric</th><th>Count</th><th>Percentage</th></tr>
            <tr><td class="metric">TOPMed Variants</td><td>', scales::comma(metrics\$topmed_n_variants), '</td><td>-</td></tr>
            <tr><td class="metric">All of Us Variants</td><td>', scales::comma(metrics\$allofus_n_variants), '</td><td>-</td></tr>
            <tr><td class="metric">Shared Variants</td><td>', scales::comma(metrics\$shared_variants), '</td><td>', 
            sprintf("%.1f%%", 100 * metrics\$shared_variants / min(metrics\$topmed_n_variants, metrics\$allofus_n_variants)), '</td></tr>
            <tr><td class="metric">TOPMed Unique</td><td>', scales::comma(metrics\$topmed_unique), '</td><td>-</td></tr>
            <tr><td class="metric">All of Us Unique</td><td>', scales::comma(metrics\$allofus_unique), '</td><td>-</td></tr>
        </table>
        
        <h3>Visualizations</h3>
        <p>See PDF report for detailed comparison plots.</p>
    </body>
    </html>
    '), file = "${sample_id}_imputation_comparison.html")
    
    cat("Imputation comparison complete\\n")
    """
}

// ============================================================================
// WORKFLOW DEFINITION
// ============================================================================

workflow MODULE6_POSTMERGE_QC {
    take:
    input_data  // Can be from Module 5 (re-imputed) or Module 4 (merged)
    
    main:
    // Step 1: Validate input source
    CHECK_INPUT_SOURCE(input_data)
    
    // Step 2: Split by chromosome for parallel MagicalRsq-X filtering
    chromosomes = Channel.of(1..22, 'X')
    
    chr_data = CHECK_INPUT_SOURCE.out.validated
        .combine(chromosomes)
        .map { sample_id, source, vcf, index, server, chr ->
            // Extract chromosome from VCF
            tuple(sample_id, source, vcf, index, server, chr)
        }
    
    // Get MagicalRsq-X scripts and models
    magicalrsqx_script = file("${params.magicalrsqx_dir}/scripts/calculate_magicalrsqx.R")
    magicalrsqx_models = file("${params.magicalrsqx_dir}/models")
    
    // Step 3: Apply MagicalRsq-X filtering per chromosome
    MAGICALRSQX_FILTER(
        chr_data,
        magicalrsqx_script,
        magicalrsqx_models
    )
    
    // Step 4: Group by sample and server, then concatenate chromosomes
    grouped_vcfs = MAGICALRSQX_FILTER.out.filtered_vcf
        .groupTuple(by: [0, 3])  // Group by sample_id and server
        .map { sample_id, vcfs, indices, server, chrs ->
            tuple(sample_id, server, vcfs.sort(), indices.sort())
        }
    
    CONCAT_CHROMOSOMES(grouped_vcfs)
    
    // Step 5: Convert to PLINK
    VCF_TO_PLINK(CONCAT_CHROMOSOMES.out.concat_vcf)
    
    // Step 6: Calculate call rates
    CALCULATE_CALL_RATES(VCF_TO_PLINK.out.plink_files)
    
    // Step 7: Filter by call rate
    FILTER_CALL_RATE(CALCULATE_CALL_RATES.out.with_missingness)
    
    // Step 8: Filter by HWE (skipped by default for admixed populations)
    // HWE assumptions are violated in admixed/structured populations
    // Enable with --skip_hwe false for homogeneous population studies
    FILTER_HWE(FILTER_CALL_RATE.out.filtered_files)

    // Conditionally route data: use HWE output if HWE ran, otherwise bypass
    hwe_output = params.skip_hwe
        ? FILTER_CALL_RATE.out.filtered_files
        : FILTER_HWE.out.filtered_files

    // Step 9: Filter heterozygosity
    FILTER_HETEROZYGOSITY(hwe_output)
    
    // Step 10: Calculate kinship
    CALCULATE_KINSHIP(FILTER_HETEROZYGOSITY.out.filtered_files)
    
    // Step 11: Remove related individuals (need imiss from Step 6)
    kinship_with_imiss = CALCULATE_KINSHIP.out.with_kinship
        .join(CALCULATE_CALL_RATES.out.with_missingness.map { 
            sample_id, server, bed, bim, fam, imiss, lmiss -> 
            tuple(sample_id, server, imiss) 
        }, by: [0, 1])
        .map { sample_id, server, bed, bim, fam, kin0, imiss ->
            tuple(sample_id, server, bed, bim, fam, kin0, imiss)
        }
    
    REMOVE_RELATEDS(
        kinship_with_imiss.map { it[0..5] },  // sample_id, server, bed, bim, fam, kin0
        kinship_with_imiss.map { it[6] }      // imiss
    )
    
    // Step 12: Calculate PCs
    CALCULATE_PCS(REMOVE_RELATEDS.out.unrelated_files)
    
    // Step 13: MAF filtering - create two branches
    MAF_FILTER_SPLIT(REMOVE_RELATEDS.out.unrelated_files)
    
    // Step 14: Generate QC summary
    qc_files = Channel.empty()  // Collect all QC reports
    
    GENERATE_QC_SUMMARY(
        MAF_FILTER_SPLIT.out.summary
            .map { sample_id, server, file -> tuple(sample_id, server, 'maf', file) }
    )
    
    // Step 15: Compare imputation performance (if both servers available)
    // Group by sample_id to compare TOPMed vs All of Us
    server_comparison = MAF_FILTER_SPLIT.out.maf_filtered
        .groupTuple(by: 0)  // Group by sample_id
        .filter { it[1].size() > 1 }  // Only if we have both servers
        .map { sample_id, servers, beds, bims, fams ->
            def topmed_idx = servers.findIndexOf { it == 'topmed' }
            def allofus_idx = servers.findIndexOf { it == 'allofus' }
            tuple(
                sample_id,
                topmed_idx >= 0 ? bims[topmed_idx] : file('NO_FILE'),
                topmed_idx >= 0 ? MAF_FILTER_SPLIT.out.summary[topmed_idx] : file('NO_FILE'),
                allofus_idx >= 0 ? bims[allofus_idx] : file('NO_FILE'),
                allofus_idx >= 0 ? MAF_FILTER_SPLIT.out.summary[allofus_idx] : file('NO_FILE')
            )
        }
    
    COMPARE_IMPUTATION_PERFORMANCE(server_comparison)
    
    emit:
    final_maf_filtered = MAF_FILTER_SPLIT.out.maf_filtered
    final_no_maf = MAF_FILTER_SPLIT.out.no_maf_filter
    pcs = CALCULATE_PCS.out.pcs
    qc_reports = GENERATE_QC_SUMMARY.out.html_report
    imputation_comparison = COMPARE_IMPUTATION_PERFORMANCE.out.html_report
}
