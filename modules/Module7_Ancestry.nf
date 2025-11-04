/*
================================================================================
    MODULE 7: ANCESTRY ANALYSIS 
================================================================================
    Comprehensive ancestry estimation including:
    - Global ancestry with GRAF-anc (genetic distance axes GD1-GD6)
      * Correctly parses 3-digit AncGroupID (major groups 1-8, subgroups 101-800)
    - ADMIXTURE analysis (K=2-12 with cross-validation)
      * Standard plots for all K values
      * SPECIAL: K=6 and K=9 organized by GRAF-anc major ancestry groups
      * Samples sorted by dominant component within each ancestry group
    - Local Ancestry Inference (LAI) - OPTIONAL:
      * RFMix v2 (default)
      * RFMix v1, FLARE, G-NOMIX (optional)
    - Principal Components Analysis colored by ancestry
    - QC metrics comparison by ancestry group
    - Imputation performance comparison by ancestry (major + subgroups)
    - Comprehensive visualizations
    
    Updates from original:
    - Fixed GRAF-anc to use correct executable and output format
    - Added major ancestry group parsing (100s: AFR, EUR, EAS, etc.)
    - Added subcontinental group parsing (101: Nigeria, 305: NE Europe, etc.)
    - New: ADMIXTURE plots organized and sorted by GRAF-anc ancestry
    - New: Imputation performance metrics stratified by ancestry
================================================================================
*/

// ============================================================================
// PROCESS DEFINITIONS
// ============================================================================

/*
 * Process 1: Prepare Data for Ancestry Analysis
 * LD prune and format for downstream tools
 */
process PREPARE_FOR_ANCESTRY {
    tag "${sample_id}_${server}"
    label 'process_medium'
    publishDir "${params.outdir}/module7/00_prep/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), path(bed), path(bim), path(fam)
    
    output:
    tuple val(sample_id), val(server),
          path("${sample_id}_${server}_pruned.bed"),
          path("${sample_id}_${server}_pruned.bim"),
          path("${sample_id}_${server}_pruned.fam"), emit: pruned_data
    tuple val(sample_id), val(server),
          path(bed), path(bim), path(fam), emit: original_data
    
    script:
    """
    echo "=== Preparing data for ancestry analysis (${server}) ==="
    
    # LD pruning for ancestry analysis
    # Use moderate pruning to retain ancestry-informative SNPs
    plink --bfile ${sample_id}_${server} \\
        --indep-pairwise 200 50 0.25 \\
        --out ${sample_id}_${server}_prune_list \\
        --threads ${task.cpus}
    
    plink --bfile ${sample_id}_${server} \\
        --extract ${sample_id}_${server}_prune_list.prune.in \\
        --make-bed \\
        --out ${sample_id}_${server}_pruned \\
        --threads ${task.cpus}
    
    echo "=== Data preparation complete ==="
    """
}

/*
 * Process 2: GRAF-anc Global Ancestry Estimation
 * CORRECTED: Uses grafanc executable with proper 3-digit AncGroupID parsing
 */
process GRAFANC_ANCESTRY {
    tag "${sample_id}_${server}"
    label 'process_high'
    publishDir "${params.outdir}/module7/01_grafanc/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), path(bed), path(bim), path(fam)
    
    output:
    tuple val(sample_id), val(server),
          path("${sample_id}_${server}_grafanc_results.txt"), emit: ancestry
    path("${sample_id}_${server}_grafanc_GD_plots.pdf"), emit: plots
    path("${sample_id}_${server}_grafanc_summary.txt"), emit: summary
    path("${sample_id}_${server}_grafanc_by_major_group.txt"), emit: major_groups
    path("${sample_id}_${server}_grafanc_by_subgroup.txt"), emit: subgroups
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    echo "=== Running GRAF-anc for ${sample_id} (${server}) ==="
    
    # Run GRAF-anc executable
    # grafanc takes: <input PLINK prefix> <output file> [options]
    grafanc ${sample_id}_${server} \\
        ${sample_id}_${server}_grafanc_results.txt \\
        --threads ${task.cpus} \\
        --maxmem ${task.memory.toMega()}
    
    echo "=== GRAF-anc complete, creating summaries and plots ==="
    
    # Create comprehensive summary and visualizations in R
    Rscript - <<'RSCRIPT'
    library(data.table)
    library(ggplot2)
    library(reshape2)
    
    # Read GRAF-anc output
    # Columns: Sample, #SNPs, GD1, GD2, GD3, EA1-EA4, AF1-AF3, EU1-EU3, 
    #          SA1-SA2, IC1-IC3, Pe, Pf, Pa, RawPe, RawPf, RawPa, AncGroupID
    anc <- fread("${sample_id}_${server}_grafanc_results.txt")
    
    # Define major ancestry group labels (first digit of 3-digit code)
    # 1XX = African, 2XX = Middle East/North Africa, 3XX = European, etc.
    major_group_labels <- c(
        "1" = "AFR - African",
        "2" = "MEN - Middle East/North Africa", 
        "3" = "EUR - European",
        "4" = "SAS - South Asian",
        "5" = "EAS - East Asian",
        "6" = "AMR - American",
        "7" = "OCN - Oceania",
        "8" = "MIX - Multi-ancestry"
    )
    
    # Define subcontinental group labels (full 3-digit codes)
    subgroup_labels <- c(
        "101" = "Nigeria", "102" = "Western Africa", "103" = "Central Africa",
        "104" = "Kenya", "105" = "Southern Africa", "106" = "Northeastern Africa",
        "107" = "African American", "108" = "Other Africa",
        "201" = "Northern Africa", "202" = "Middle East 1", "203" = "Middle East 2",
        "301" = "Finland", "302" = "Northern Europe", "303" = "Western Europe",
        "304" = "Southern Europe", "305" = "Northeastern Europe", 
        "306" = "Southeastern Europe", "307" = "Balkans", "308" = "Other Europe",
        "401" = "Asian India", "402" = "Gujarati India", 
        "403" = "Northern South Asia", "404" = "Sri Lanka", "405" = "Bangladesh",
        "501" = "Japan Ryukyu", "502" = "Japan Main Islands", "503" = "Korea",
        "504" = "Northern Asia", "505" = "Northern China 1", "506" = "Northern China 2",
        "507" = "Southern China 1", "508" = "Southern China 2", 
        "509" = "Southeast Asia", "510" = "Thailand", "511" = "Other East Asia",
        "601" = "Latin American 1", "602" = "Latin American 2", 
        "603" = "Native American",
        "700" = "Oceania",
        "800" = "Multiracial"
    )
    
    # Extract major ancestry group (first digit = 100s place)
    anc[, MajorGroup := floor(AncGroupID / 100)]
    anc[, MajorGroupLabel := major_group_labels[as.character(MajorGroup)]]
    
    # Add subgroup label
    anc[, SubgroupLabel := subgroup_labels[as.character(AncGroupID)]]
    # If subgroup not in our list, use the AncGroupID
    anc[is.na(SubgroupLabel), SubgroupLabel := paste0("Group_", AncGroupID)]
    
    # ========================================================================
    # SUMMARY 1: Overall statistics
    # ========================================================================
    cat("GRAF-anc Ancestry Analysis Summary\\n", 
        file = "${sample_id}_${server}_grafanc_summary.txt")
    cat("===================================\\n", 
        file = "${sample_id}_${server}_grafanc_summary.txt", append = TRUE)
    cat(sprintf("Sample: %s\\n", "${sample_id}"), 
        file = "${sample_id}_${server}_grafanc_summary.txt", append = TRUE)
    cat(sprintf("Server: %s\\n", "${server}"), 
        file = "${sample_id}_${server}_grafanc_summary.txt", append = TRUE)
    cat(sprintf("Total samples: %d\\n", nrow(anc)), 
        file = "${sample_id}_${server}_grafanc_summary.txt", append = TRUE)
    cat(sprintf("Mean ancestry SNPs genotyped: %.0f (SD: %.0f)\\n", 
               mean(anc[['#SNPs']], na.rm = TRUE), 
               sd(anc[['#SNPs']], na.rm = TRUE)),
        file = "${sample_id}_${server}_grafanc_summary.txt", append = TRUE)
    
    # ========================================================================
    # SUMMARY 2: Counts by Major Ancestry Group (Continental Level)
    # ========================================================================
    major_counts <- anc[, .N, by = .(MajorGroup, MajorGroupLabel)]
    setorder(major_counts, MajorGroup)
    major_counts[, Percentage := round(100 * N / sum(N), 2)]
    
    cat("\\n\\nMAJOR ANCESTRY GROUPS (Continental Level)\\n", 
        file = "${sample_id}_${server}_grafanc_summary.txt", append = TRUE)
    cat("==========================================\\n", 
        file = "${sample_id}_${server}_grafanc_summary.txt", append = TRUE)
    cat("1=AFR, 2=MEN, 3=EUR, 4=SAS, 5=EAS, 6=AMR, 7=OCN, 8=MIX\\n", 
        file = "${sample_id}_${server}_grafanc_summary.txt", append = TRUE)
    
    write.table(major_counts, 
               file = "${sample_id}_${server}_grafanc_summary.txt",
               append = TRUE, quote = FALSE, row.names = FALSE, sep = "\\t")
    
    # Save major groups table separately for downstream use
    fwrite(major_counts, "${sample_id}_${server}_grafanc_by_major_group.txt", 
           sep = "\\t")
    
    # ========================================================================
    # SUMMARY 3: Counts by Subcontinental Group
    # ========================================================================
    sub_counts <- anc[, .N, by = .(AncGroupID, SubgroupLabel, MajorGroupLabel)]
    setorder(sub_counts, AncGroupID)
    sub_counts[, Percentage := round(100 * N / sum(N), 2)]
    
    cat("\\n\\nSUBCONTINENTAL ANCESTRY GROUPS\\n", 
        file = "${sample_id}_${server}_grafanc_summary.txt", append = TRUE)
    cat("================================\\n", 
        file = "${sample_id}_${server}_grafanc_summary.txt", append = TRUE)
    
    write.table(sub_counts, 
               file = "${sample_id}_${server}_grafanc_summary.txt",
               append = TRUE, quote = FALSE, row.names = FALSE, sep = "\\t")
    
    # Save subgroups table separately for downstream use
    fwrite(sub_counts, "${sample_id}_${server}_grafanc_by_subgroup.txt", 
           sep = "\\t")
    
    # ========================================================================
    # SUMMARY 4: Average ancestry proportions by major group
    # ========================================================================
    prop_summary <- anc[, .(
        Mean_Pe = round(mean(Pe, na.rm = TRUE), 4),
        SD_Pe = round(sd(Pe, na.rm = TRUE), 4),
        Mean_Pf = round(mean(Pf, na.rm = TRUE), 4),
        SD_Pf = round(sd(Pf, na.rm = TRUE), 4),
        Mean_Pa = round(mean(Pa, na.rm = TRUE), 4),
        SD_Pa = round(sd(Pa, na.rm = TRUE), 4),
        N = .N
    ), by = MajorGroupLabel]
    
    cat("\\n\\nAVERAGE ANCESTRY PROPORTIONS BY MAJOR GROUP\\n", 
        file = "${sample_id}_${server}_grafanc_summary.txt", append = TRUE)
    cat("============================================\\n", 
        file = "${sample_id}_${server}_grafanc_summary.txt", append = TRUE)
    cat("Pe = European, Pf = African, Pa = East Asian\\n",
        file = "${sample_id}_${server}_grafanc_summary.txt", append = TRUE)
    
    write.table(prop_summary, 
               file = "${sample_id}_${server}_grafanc_summary.txt",
               append = TRUE, quote = FALSE, row.names = FALSE, sep = "\\t")
    
    # ========================================================================
    # VISUALIZATIONS
    # ========================================================================
    pdf("${sample_id}_${server}_grafanc_GD_plots.pdf", width = 16, height = 12)
    
    # Color palette for major groups
    major_colors <- c(
        "AFR - African" = "#E41A1C",
        "MEN - Middle East/North Africa" = "#FF7F00", 
        "EUR - European" = "#377EB8",
        "SAS - South Asian" = "#4DAF4A",
        "EAS - East Asian" = "#984EA3",
        "AMR - American" = "#FFFF33",
        "OCN - Oceania" = "#A65628",
        "MIX - Multi-ancestry" = "#999999"
    )
    
    # Plot 1: GD1 vs GD2 - Main population structure plot
    p1 <- ggplot(anc, aes(x = GD1, y = GD2, color = MajorGroupLabel)) +
        geom_point(alpha = 0.7, size = 2) +
        scale_color_manual(values = major_colors, na.value = "grey50") +
        theme_bw(base_size = 12) +
        labs(title = "GRAF-anc: GD1 vs GD2 (Global Population Structure)",
             subtitle = paste("${sample_id} (${server}) -", nrow(anc), "samples"),
             x = "GD1 (Genetic Distance 1)",
             y = "GD2 (Genetic Distance 2)",
             color = "Major\\nAncestry Group") +
        theme(legend.position = "right",
              plot.title = element_text(face = "bold"))
    print(p1)
    
    # Plot 2: GD1 vs GD3 - Useful for Latin American/Native American separation
    p2 <- ggplot(anc, aes(x = GD1, y = GD3, color = MajorGroupLabel)) +
        geom_point(alpha = 0.7, size = 2) +
        scale_color_manual(values = major_colors, na.value = "grey50") +
        theme_bw(base_size = 12) +
        labs(title = "GRAF-anc: GD1 vs GD3",
             subtitle = "Useful for American/Native American populations",
             x = "GD1 (Genetic Distance 1)",
             y = "GD3 (Genetic Distance 3)",
             color = "Major\\nAncestry Group") +
        theme(legend.position = "right")
    print(p2)
    
    # Plot 3: GD1 vs IC1 - Separates South Asians from Latin Americans
    if ("IC1" %in% names(anc)) {
        p3 <- ggplot(anc, aes(x = GD1, y = IC1, color = MajorGroupLabel)) +
            geom_point(alpha = 0.7, size = 2) +
            scale_color_manual(values = major_colors, na.value = "grey50") +
            theme_bw(base_size = 12) +
            labs(title = "GRAF-anc: GD1 vs IC1",
                 subtitle = "Separates South Asians from Latin Americans",
                 x = "GD1 (Genetic Distance 1)",
                 y = "IC1 (Inter-Continental 1)",
                 color = "Major\\nAncestry Group") +
            theme(legend.position = "right")
        print(p3)
    }
    
    # Plot 4: Ancestry proportions (Pe, Pf, Pa) by major group
    if (all(c("Pe", "Pf", "Pa") %in% names(anc))) {
        anc_props <- anc[, .(Sample, MajorGroupLabel, Pe, Pf, Pa)]
        anc_long <- melt(anc_props, 
                        id.vars = c("Sample", "MajorGroupLabel"),
                        variable.name = "Component",
                        value.name = "Proportion")
        
        # Create labels for components
        anc_long[, ComponentLabel := factor(Component, 
                                            levels = c("Pe", "Pf", "Pa"),
                                            labels = c("European", "African", "East Asian"))]
        
        p4 <- ggplot(anc_long, aes(x = MajorGroupLabel, y = Proportion, 
                                   fill = ComponentLabel)) +
            geom_boxplot(outlier.size = 0.5) +
            facet_wrap(~ ComponentLabel, ncol = 3) +
            theme_bw(base_size = 12) +
            labs(title = "Ancestry Proportions by Major Group",
                 subtitle = "Pe (European), Pf (African), Pa (East Asian)",
                 x = "Major Ancestry Group",
                 y = "Proportion",
                 fill = "Ancestry\\nComponent") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "bottom") +
            scale_fill_manual(values = c("European" = "#377EB8",
                                        "African" = "#E41A1C",
                                        "East Asian" = "#984EA3"))
        print(p4)
    }
    
    # Plot 5: Subcontinental structure (if multiple subgroups present)
    n_subgroups <- length(unique(anc\$SubgroupLabel))
    if (n_subgroups > 1 && n_subgroups <= 20) {
        
        p5 <- ggplot(anc, aes(x = GD1, y = GD2, color = SubgroupLabel)) +
            geom_point(alpha = 0.7, size = 2) +
            theme_bw(base_size = 12) +
            labs(title = "GRAF-anc: Subcontinental Population Structure",
                 subtitle = paste("${sample_id} (${server})"),
                 x = "GD1 (Genetic Distance 1)",
                 y = "GD2 (Genetic Distance 2)",
                 color = "Subgroup") +
            theme(legend.position = "right",
                  legend.text = element_text(size = 8))
        print(p5)
    }
    
    # Plot 6: Sample size barplot by major group
    p6 <- ggplot(major_counts, aes(x = reorder(MajorGroupLabel, -N), y = N,
                                   fill = MajorGroupLabel)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = paste0(N, "\\n(", Percentage, "%)")), 
                 vjust = -0.5, size = 3.5) +
        scale_fill_manual(values = major_colors, na.value = "grey50") +
        theme_bw(base_size = 12) +
        labs(title = "Sample Distribution by Major Ancestry Group",
             subtitle = "${sample_id} (${server})",
             x = "Major Ancestry Group",
             y = "Number of Samples") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none")
    print(p6)
    
    # Plot 7: Top 15 subcontinental groups
    if (nrow(sub_counts) > 1) {
        top_subs <- head(sub_counts[order(-N)], 15)
        
        p7 <- ggplot(top_subs, aes(x = reorder(SubgroupLabel, N), y = N,
                                   fill = MajorGroupLabel)) +
            geom_bar(stat = "identity") +
            geom_text(aes(label = N), hjust = -0.2, size = 3) +
            scale_fill_manual(values = major_colors, na.value = "grey50") +
            coord_flip() +
            theme_bw(base_size = 11) +
            labs(title = "Top 15 Subcontinental Ancestry Groups",
                 subtitle = "${sample_id} (${server})",
                 x = "Subcontinental Group",
                 y = "Number of Samples",
                 fill = "Major\\nGroup") +
            theme(legend.position = "right")
        print(p7)
    }
    
    dev.off()
    
    cat("\\n=== GRAF-anc analysis complete ===\\n")
    cat("Output files created:\\n")
    cat("  - ${sample_id}_${server}_grafanc_results.txt (full results)\\n")
    cat("  - ${sample_id}_${server}_grafanc_summary.txt (summary statistics)\\n")
    cat("  - ${sample_id}_${server}_grafanc_by_major_group.txt (counts by major group)\\n")
    cat("  - ${sample_id}_${server}_grafanc_by_subgroup.txt (counts by subgroup)\\n")
    cat("  - ${sample_id}_${server}_grafanc_GD_plots.pdf (visualizations)\\n")
RSCRIPT
    
    echo "=== GRAF-anc summary and visualizations complete ==="
    """
}

/*
 * Process 3: ADMIXTURE Analysis with Cross-Validation
 * Run for K=2 to K=12 (or user-specified range)
 */
process RUN_ADMIXTURE {
    tag "${sample_id}_${server}_K${k}"
    label 'process_high'
    publishDir "${params.outdir}/module7/02_admixture/${server}/K${k}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), path(bed), path(bim), path(fam), val(k)
    
    output:
    tuple val(sample_id), val(server), val(k),
          path("${sample_id}_${server}.${k}.Q"),
          path("${sample_id}_${server}.${k}.P"), emit: results
    path("${sample_id}_${server}_K${k}_CV.txt"), emit: cv_error
    
    script:
    def cv_folds = params.admixture_cv_folds ?: 5
    """
    echo "=== Running ADMIXTURE K=${k} for ${sample_id} (${server}) ==="
    
    # Run ADMIXTURE with cross-validation
    admixture --cv=${cv_folds} \\
        ${bed} \\
        ${k} \\
        -j${task.cpus} \\
        | tee ${sample_id}_${server}_K${k}_log.txt
    
    # Extract CV error
    grep "CV error" ${sample_id}_${server}_K${k}_log.txt > ${sample_id}_${server}_K${k}_CV.txt || true
    
    # Rename output files to include sample and server info
    mv ${sample_id}_${server}.${k}.Q ${sample_id}_${server}.${k}.Q
    mv ${sample_id}_${server}.${k}.P ${sample_id}_${server}.${k}.P
    
    echo "=== ADMIXTURE K=${k} complete ==="
    """
}

/*
 * Process 4: Aggregate ADMIXTURE Results and Create Standard Visualizations
 */
process SUMMARIZE_ADMIXTURE {
    tag "${sample_id}_${server}"
    label 'process_medium'
    publishDir "${params.outdir}/module7/02_admixture/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), val(k_values), path(q_files), path(p_files)
    path(cv_files)
    path(fam)
    
    output:
    path("${sample_id}_${server}_admixture_summary.txt"), emit: summary
    path("${sample_id}_${server}_admixture_cv_plot.pdf"), emit: cv_plot
    path("${sample_id}_${server}_admixture_K*.pdf"), emit: barplots
    tuple val(sample_id), val(server), path(q_files), path(fam), emit: for_grafanc_org
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(data.table)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    # Read sample IDs from FAM
    fam_data <- fread("${fam}", header = FALSE)
    sample_ids <- fam_data\$V2
    
    # Collect all K values
    k_values <- c(${k_values.join(',')})
    
    # Read CV errors
    cv_errors <- data.frame(K = integer(), CV_Error = numeric())
    cv_files <- list.files(pattern = "_CV.txt\$", full.names = TRUE)
    
    for (cv_file in cv_files) {
        k <- as.integer(gsub(".*_K([0-9]+)_CV.txt", "\\\\1", cv_file))
        cv_line <- readLines(cv_file)
        if (length(cv_line) > 0) {
            cv_error <- as.numeric(sub(".*\\\\(K=.*\\\\):\\\\s+([0-9.]+)", "\\\\1", cv_line[1]))
            cv_errors <- rbind(cv_errors, data.frame(K = k, CV_Error = cv_error))
        }
    }
    
    # Find optimal K and plot CV errors
    if (nrow(cv_errors) > 0) {
        optimal_k <- cv_errors[which.min(cv_errors\$CV_Error), "K"]
        
        # Plot CV errors
        pdf("${sample_id}_${server}_admixture_cv_plot.pdf", width = 10, height = 6)
        p_cv <- ggplot(cv_errors, aes(x = K, y = CV_Error)) +
            geom_line(size = 1) +
            geom_point(size = 3) +
            geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red") +
            theme_bw() +
            labs(title = "ADMIXTURE Cross-Validation Error",
                 subtitle = sprintf("Optimal K = %d (lowest CV error)", optimal_k),
                 x = "Number of Ancestral Populations (K)",
                 y = "Cross-Validation Error")
        print(p_cv)
        dev.off()
        
        # Write summary
        cat("ADMIXTURE Cross-Validation Summary\\n", 
            file = "${sample_id}_${server}_admixture_summary.txt")
        cat("===================================\\n", 
            file = "${sample_id}_${server}_admixture_summary.txt", append = TRUE)
        cat(sprintf("Sample: %s\\n", "${sample_id}"), 
            file = "${sample_id}_${server}_admixture_summary.txt", append = TRUE)
        cat(sprintf("Server: %s\\n", "${server}"), 
            file = "${sample_id}_${server}_admixture_summary.txt", append = TRUE)
        cat(sprintf("Optimal K: %d\\n", optimal_k), 
            file = "${sample_id}_${server}_admixture_summary.txt", append = TRUE)
        cat(sprintf("Minimum CV Error: %.6f\\n", min(cv_errors\$CV_Error)), 
            file = "${sample_id}_${server}_admixture_summary.txt", append = TRUE)
        cat("\\nCV Errors by K:\\n", 
            file = "${sample_id}_${server}_admixture_summary.txt", append = TRUE)
        write.table(cv_errors, 
                   file = "${sample_id}_${server}_admixture_summary.txt", 
                   append = TRUE, quote = FALSE, row.names = FALSE, sep = "\\t")
    }
    
    # Create barplots for each K
    q_files <- list.files(pattern = "\\\\.Q\$", full.names = TRUE)
    
    for (q_file in q_files) {
        k <- as.integer(sub(".*\\\\.([0-9]+)\\\\.Q", "\\\\1", q_file))
        
        # Read Q matrix
        q_matrix <- as.matrix(read.table(q_file))
        colnames(q_matrix) <- paste0("Pop", 1:k)
        rownames(q_matrix) <- sample_ids
        
        # Convert to long format
        q_long <- melt(q_matrix)
        colnames(q_long) <- c("Sample", "Population", "Proportion")
        q_long\$Sample <- factor(q_long\$Sample, levels = sample_ids)
        
        # Create stacked barplot
        pdf(sprintf("${sample_id}_${server}_admixture_K%d.pdf", k), 
            width = max(14, length(sample_ids) * 0.02), height = 6)
        
        # Get colors
        if (k <= 12) {
            colors <- brewer.pal(max(3, k), "Set3")[1:k]
        } else {
            colors <- rainbow(k)
        }
        
        p <- ggplot(q_long, aes(x = Sample, y = Proportion, fill = Population)) +
            geom_bar(stat = "identity", width = 1) +
            theme_minimal() +
            labs(title = sprintf("ADMIXTURE Analysis (K=%d)", k),
                 subtitle = "${sample_id} (${server})",
                 y = "Ancestry Proportion") +
            scale_fill_manual(values = colors) +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  panel.grid = element_blank())
        print(p)
        dev.off()
    }
    
    cat("\\nADMIXTURE analysis complete for all K values\\n")
    """
}

/*
 * Process 4B: Create ADMIXTURE Plots Organized by GRAF-anc Major Ancestry Groups
 * For K=6 and K=9, split by major ancestry, sorted by dominant component
 * NEW PROCESS
 */
process ADMIXTURE_BY_GRAFANC {
    tag "${sample_id}_${server}"
    label 'process_medium'
    publishDir "${params.outdir}/module7/02_admixture/${server}/by_ancestry", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), path(q_files), path(fam)
    path(grafanc_results)  // GRAF-anc results file
    
    output:
    path("${sample_id}_${server}_admixture_by_ancestry_K*.pdf"), emit: plots
    path("${sample_id}_${server}_admixture_ancestry_summary.txt"), emit: summary
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(data.table)
    library(ggplot2)
    library(RColorBrewer)
    library(cowplot)
    
    # Read sample IDs from FAM
    fam_data <- fread("${fam}", header = FALSE)
    setnames(fam_data, c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))
    
    # Read GRAF-anc results
    grafanc <- fread("${grafanc_results}")
    
    # Extract major ancestry group (first digit of AncGroupID)
    grafanc[, MajorGroup := floor(AncGroupID / 100)]
    
    # Define major group labels
    major_labels <- c(
        "1" = "AFR", "2" = "MEN", "3" = "EUR", "4" = "SAS",
        "5" = "EAS", "6" = "AMR", "7" = "OCN", "8" = "MIX"
    )
    grafanc[, MajorGroupLabel := major_labels[as.character(MajorGroup)]]
    
    # Match GRAF-anc Sample column to FAM IID
    setnames(grafanc, "Sample", "IID", skip_absent = TRUE)
    
    # Merge FAM with GRAF-anc ancestry
    merged <- merge(fam_data, grafanc[, .(IID, MajorGroup, MajorGroupLabel, 
                                          AncGroupID, Pe, Pf, Pa)], 
                   by = "IID", all.x = TRUE)
    
    # Process Q files for K=6 and K=9
    q_files <- list.files(pattern = "\\\\.Q\$", full.names = TRUE)
    
    # Function to create organized ADMIXTURE plot
    create_ancestry_organized_plot <- function(k, q_file) {
        cat(sprintf("\\nProcessing K=%d...\\n", k))
        
        # Read Q matrix
        q_matrix <- as.matrix(read.table(q_file))
        colnames(q_matrix) <- paste0("Pop", 1:k)
        
        # Combine with sample info and ancestry
        q_df <- data.table(
            IID = merged\$IID,
            MajorGroupLabel = merged\$MajorGroupLabel,
            q_matrix
        )
        
        # Remove samples with missing ancestry (if any)
        q_df <- q_df[!is.na(MajorGroupLabel)]
        
        # For each sample, find dominant ancestry component
        pop_cols <- paste0("Pop", 1:k)
        q_df[, DominantPop := pop_cols[apply(.SD, 1, which.max)], .SDcols = pop_cols]
        q_df[, DominantValue := apply(.SD, 1, max), .SDcols = pop_cols]
        
        # Within each major ancestry group, sort by:
        # 1. Dominant population (alphabetically)
        # 2. Dominant value (descending - highest to lowest)
        setorder(q_df, MajorGroupLabel, DominantPop, -DominantValue)
        
        # Add order index within each major group
        q_df[, OrderIndex := 1:.N, by = MajorGroupLabel]
        q_df[, SampleOrder := paste(MajorGroupLabel, sprintf("%05d", OrderIndex), sep = "_")]
        
        # Convert to long format for plotting
        q_long <- melt(q_df, 
                      id.vars = c("IID", "MajorGroupLabel", "SampleOrder", 
                                 "DominantPop", "DominantValue", "OrderIndex"),
                      variable.name = "Population", 
                      value.name = "Proportion")
        
        # Create ordered factor for proper plotting
        q_long[, SampleOrder := factor(SampleOrder, levels = unique(q_df\$SampleOrder))]
        
        # Get sample counts by major group
        group_counts <- q_df[, .N, by = MajorGroupLabel]
        setorder(group_counts, MajorGroupLabel)
        
        # Define colors for K populations
        if (k <= 9) {
            colors <- brewer.pal(max(3, k), "Set1")[1:k]
        } else {
            colors <- rainbow(k)
        }
        names(colors) <- paste0("Pop", 1:k)
        
        # Create faceted plot by major ancestry group
        # Use free_x scale so each facet only shows its samples
        p <- ggplot(q_long, aes(x = SampleOrder, y = Proportion, fill = Population)) +
            geom_bar(stat = "identity", width = 1) +
            facet_grid(. ~ MajorGroupLabel, scales = "free_x", space = "free_x") +
            scale_fill_manual(values = colors) +
            theme_minimal(base_size = 12) +
            labs(title = sprintf("ADMIXTURE (K=%d) Organized by GRAF-anc Ancestry", k),
                 subtitle = sprintf("%s (%s) - Sorted by dominant component within ancestry group", 
                                   "${sample_id}", "${server}"),
                 y = "Ancestry Proportion",
                 x = "Samples") +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  panel.grid = element_blank(),
                  strip.text = element_text(face = "bold", size = 10),
                  strip.background = element_rect(fill = "grey90", color = "black"),
                  panel.spacing = unit(0.5, "lines"),
                  legend.position = "bottom",
                  plot.title = element_text(face = "bold", size = 14))
        
        return(list(plot = p, data = q_df, counts = group_counts))
    }
    
    # Create summary file
    cat("ADMIXTURE Results Organized by GRAF-anc Ancestry\\n", 
        file = "${sample_id}_${server}_admixture_ancestry_summary.txt")
    cat("================================================\\n", 
        file = "${sample_id}_${server}_admixture_ancestry_summary.txt", append = TRUE)
    cat(sprintf("Sample: %s\\n", "${sample_id}"), 
        file = "${sample_id}_${server}_admixture_ancestry_summary.txt", append = TRUE)
    cat(sprintf("Server: %s\\n\\n", "${server}"), 
        file = "${sample_id}_${server}_admixture_ancestry_summary.txt", append = TRUE)
    
    # Process K=6
    k6_file <- q_files[grep("\\\\.6\\\\.Q\$", q_files)]
    if (length(k6_file) > 0) {
        result_k6 <- create_ancestry_organized_plot(6, k6_file)
        
        pdf("${sample_id}_${server}_admixture_by_ancestry_K6.pdf", 
            width = 20, height = 8)
        print(result_k6\$plot)
        dev.off()
        
        cat("K=6 Sample Distribution by Major Ancestry Group:\\n", 
            file = "${sample_id}_${server}_admixture_ancestry_summary.txt", append = TRUE)
        write.table(result_k6\$counts, 
                   file = "${sample_id}_${server}_admixture_ancestry_summary.txt",
                   append = TRUE, quote = FALSE, row.names = FALSE, sep = "\\t")
        
        # Calculate average admixture proportions by ancestry group
        avg_props <- result_k6\$data[, lapply(.SD, mean), 
                                    .SDcols = paste0("Pop", 1:6), 
                                    by = MajorGroupLabel]
        cat("\\nK=6 Average Admixture Proportions by Ancestry Group:\\n", 
            file = "${sample_id}_${server}_admixture_ancestry_summary.txt", append = TRUE)
        write.table(avg_props, 
                   file = "${sample_id}_${server}_admixture_ancestry_summary.txt",
                   append = TRUE, quote = FALSE, row.names = FALSE, sep = "\\t")
    }
    
    # Process K=9
    k9_file <- q_files[grep("\\\\.9\\\\.Q\$", q_files)]
    if (length(k9_file) > 0) {
        result_k9 <- create_ancestry_organized_plot(9, k9_file)
        
        pdf("${sample_id}_${server}_admixture_by_ancestry_K9.pdf", 
            width = 20, height = 8)
        print(result_k9\$plot)
        dev.off()
        
        cat("\\n\\nK=9 Sample Distribution by Major Ancestry Group:\\n", 
            file = "${sample_id}_${server}_admixture_ancestry_summary.txt", append = TRUE)
        write.table(result_k9\$counts, 
                   file = "${sample_id}_${server}_admixture_ancestry_summary.txt",
                   append = TRUE, quote = FALSE, row.names = FALSE, sep = "\\t")
        
        # Calculate average admixture proportions by ancestry group
        avg_props <- result_k9\$data[, lapply(.SD, mean), 
                                    .SDcols = paste0("Pop", 1:9), 
                                    by = MajorGroupLabel]
        cat("\\nK=9 Average Admixture Proportions by Ancestry Group:\\n", 
            file = "${sample_id}_${server}_admixture_ancestry_summary.txt", append = TRUE)
        write.table(avg_props, 
                   file = "${sample_id}_${server}_admixture_ancestry_summary.txt",
                   append = TRUE, quote = FALSE, row.names = FALSE, sep = "\\t")
    }
    
    cat("\\n\\nADMIXTURE plots organized by GRAF-anc ancestry complete\\n")
    cat("Samples are sorted within each ancestry group by their dominant component\\n")
    cat("This creates a 'waterfall' effect showing the gradient from highest to lowest\\n")
    cat("for the most common component in each ancestry group.\\n")
    """
}

/*
 * Process 5: Prepare for Local Ancestry Inference (OPTIONAL)
 * Phase data if not already phased, format for LAI tools
 */
process PREPARE_FOR_LAI {
    tag "${sample_id}_${server}"
    label 'process_high'
    publishDir "${params.outdir}/module7/03_local_ancestry/00_prep/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), path(vcf), path(index)
    
    output:
    tuple val(sample_id), val(server),
          path("${sample_id}_${server}_phased.vcf.gz"),
          path("${sample_id}_${server}_phased.vcf.gz.csi"), emit: phased_vcf
    
    when:
    params.run_lai
    
    script:
    """
    echo "=== Preparing for LAI: ${sample_id} (${server}) ==="
    
    # Check if already phased (from imputation servers)
    bcftools query -l ${vcf} | head -5
    
    # Most imputation servers return phased data
    # If not phased, we would need to phase here with SHAPEIT4 or Eagle
    
    # For now, assume imputation servers provide phased data
    # Just create symbolic link
    ln -s ${vcf} ${sample_id}_${server}_phased.vcf.gz
    ln -s ${index} ${sample_id}_${server}_phased.vcf.gz.csi
    
    echo "=== LAI preparation complete ==="
    """
}

/*
 * Process 6: RFMix v2 Local Ancestry Inference (DEFAULT - OPTIONAL)
 */
process RFMIX_V2 {
    tag "${sample_id}_${server}_chr${chr}"
    label 'process_high'
    publishDir "${params.outdir}/module7/03_local_ancestry/rfmix_v2/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), path(vcf), path(index), val(chr)
    path(reference_vcf)
    path(reference_sample_map)
    path(genetic_map)
    
    when:
    params.lai_tool == 'rfmix_v2' || params.lai_tool == 'all'
    
    output:
    tuple val(sample_id), val(server), val(chr),
          path("${sample_id}_${server}_chr${chr}.rfmix.Q"),
          path("${sample_id}_${server}_chr${chr}.msp.tsv"), emit: results
    path("${sample_id}_${server}_chr${chr}.rfmix.log"), emit: log
    
    script:
    """
    echo "=== Running RFMix v2 for chr${chr} (${server}) ==="
    
    # RFMix v2 command
    rfmix \\
        -f ${vcf} \\
        -r ${reference_vcf} \\
        -m ${reference_sample_map} \\
        -g ${genetic_map} \\
        -o ${sample_id}_${server}_chr${chr} \\
        --chromosome=${chr} \\
        -e ${params.rfmix_em_iterations ?: 5} \\
        --reanalyze-reference \\
        --n-threads=${task.cpus} \\
        2>&1 | tee ${sample_id}_${server}_chr${chr}.rfmix.log
    
    echo "=== RFMix v2 complete for chr${chr} ==="
    """
}

/*
 * Process 7: Summarize LAI Results (OPTIONAL)
 */
process SUMMARIZE_LAI {
    tag "${sample_id}_${server}"
    label 'process_medium'
    publishDir "${params.outdir}/module7/03_local_ancestry/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), path(lai_files)
    path(fam)
    
    when:
    params.run_lai
    
    output:
    path("${sample_id}_${server}_lai_summary.txt"), emit: summary
    path("${sample_id}_${server}_lai_genome_wide.pdf"), emit: plot optional true
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(data.table)
    library(ggplot2)
    
    # Read sample IDs
    fam_data <- fread("${fam}", header = FALSE)
    
    # Process LAI results (RFMix v2 format)
    lai_files <- list.files(pattern = ".msp.tsv\$", full.names = TRUE)
    
    if (length(lai_files) > 0) {
        # Read and combine all chromosomes
        lai_data <- rbindlist(lapply(lai_files, fread), fill = TRUE)
        
        cat("Local Ancestry Inference Summary\\n", 
            file = "${sample_id}_${server}_lai_summary.txt")
        cat("================================\\n", 
            file = "${sample_id}_${server}_lai_summary.txt", append = TRUE)
        cat(sprintf("Total samples: %d\\n", nrow(fam_data)), 
            file = "${sample_id}_${server}_lai_summary.txt", append = TRUE)
        cat(sprintf("LAI files processed: %d\\n", length(lai_files)), 
            file = "${sample_id}_${server}_lai_summary.txt", append = TRUE)
        
        # Create genome-wide LAI visualization if possible
        pdf("${sample_id}_${server}_lai_genome_wide.pdf", width = 16, height = 10)
        plot(1:10, main = "Local Ancestry Inference Results\\n(See individual chromosome files)")
        dev.off()
    } else {
        cat("No LAI results found\\n", 
            file = "${sample_id}_${server}_lai_summary.txt")
    }
    """
}

/*
 * Process 8: Create PCA Plots Colored by Ancestry (UPDATED)
 */
process PCA_BY_ANCESTRY {
    tag "${sample_id}_${server}"
    label 'process_low'
    publishDir "${params.outdir}/module7/04_pca_colored/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), path(eigenvec), path(eigenval)
    path(ancestry_file)  // From GRAF-anc (corrected format)
    
    output:
    path("${sample_id}_${server}_pca_by_ancestry.pdf"), emit: plot
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(data.table)
    library(ggplot2)
    library(cowplot)
    
    # Read PCs
    pcs <- fread("${eigenvec}")
    eigenval <- scan("${eigenval}")
    pve <- eigenval / sum(eigenval) * 100
    
    # Read GRAF-anc results (correct format)
    anc <- fread("${ancestry_file}")
    
    # Extract major ancestry group
    anc[, MajorGroup := floor(AncGroupID / 100)]
    
    # Define labels
    major_labels <- c(
        "1" = "AFR", "2" = "MEN", "3" = "EUR", "4" = "SAS",
        "5" = "EAS", "6" = "AMR", "7" = "OCN", "8" = "MIX"
    )
    anc[, MajorGroupLabel := major_labels[as.character(MajorGroup)]]
    
    # Merge with PCs (GRAF-anc uses "Sample" column)
    # Match to IID from PCs
    if ("#FID" %in% names(pcs)) {
        setnames(pcs, "#FID", "FID")
    }
    
    # GRAF-anc Sample column should match IID
    setnames(anc, "Sample", "IID", skip_absent = TRUE)
    
    merged <- merge(pcs, anc, by = "IID", all.x = TRUE)
    
    # Create plots colored by ancestry
    pdf("${sample_id}_${server}_pca_by_ancestry.pdf", width = 16, height = 12)
    
    # Define colors
    anc_colors <- c(
        "AFR" = "#E41A1C", "MEN" = "#FF7F00", "EUR" = "#377EB8",
        "SAS" = "#4DAF4A", "EAS" = "#984EA3", "AMR" = "#FFFF33",
        "OCN" = "#A65628", "MIX" = "#999999"
    )
    
    # PC1 vs PC2
    p1 <- ggplot(merged, aes(x = PC1, y = PC2, color = MajorGroupLabel)) +
        geom_point(alpha = 0.7, size = 2) +
        scale_color_manual(values = anc_colors, na.value = "grey50") +
        theme_bw(base_size = 14) +
        labs(title = "PC1 vs PC2 (Colored by GRAF-anc Major Ancestry)",
             subtitle = "${sample_id} (${server})",
             x = sprintf("PC1 (%.2f%% variance)", pve[1]),
             y = sprintf("PC2 (%.2f%% variance)", pve[2]),
             color = "Ancestry") +
        theme(legend.position = "right",
              plot.title = element_text(face = "bold", size = 16))
    
    # PC3 vs PC4
    p2 <- ggplot(merged, aes(x = PC3, y = PC4, color = MajorGroupLabel)) +
        geom_point(alpha = 0.7, size = 2) +
        scale_color_manual(values = anc_colors, na.value = "grey50") +
        theme_bw(base_size = 14) +
        labs(title = "PC3 vs PC4 (Colored by GRAF-anc Major Ancestry)",
             subtitle = "${sample_id} (${server})",
             x = sprintf("PC3 (%.2f%% variance)", pve[3]),
             y = sprintf("PC4 (%.2f%% variance)", pve[4]),
             color = "Ancestry") +
        theme(legend.position = "right",
              plot.title = element_text(face = "bold", size = 16))
    
    # Combine plots
    print(p1)
    print(p2)
    
    # Additional plot: PC1 vs PC2 colored by GD1 (continuous)
    if ("GD1" %in% names(merged)) {
        p3 <- ggplot(merged, aes(x = PC1, y = PC2, color = GD1)) +
            geom_point(alpha = 0.7, size = 2) +
            scale_color_viridis_c() +
            theme_bw(base_size = 14) +
            labs(title = "PC1 vs PC2 (Colored by GRAF-anc GD1)",
                 subtitle = "Continuous genetic distance measure",
                 x = sprintf("PC1 (%.2f%% variance)", pve[1]),
                 y = sprintf("PC2 (%.2f%% variance)", pve[2]),
                 color = "GD1") +
            theme(legend.position = "right")
        print(p3)
    }
    
    dev.off()
    
    cat("PCA plots colored by GRAF-anc ancestry created successfully\\n")
    """
}

/*
 * Process 9: Compare QC Metrics by Ancestry Group
 */
process QC_BY_ANCESTRY {
    tag "${sample_id}_${server}"
    label 'process_medium'
    publishDir "${params.outdir}/module7/05_qc_by_ancestry/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), path(bed), path(bim), path(fam)
    path(ancestry_file)
    
    output:
    path("${sample_id}_${server}_qc_by_ancestry.txt"), emit: summary
    path("${sample_id}_${server}_qc_by_ancestry.pdf"), emit: plot
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(data.table)
    library(ggplot2)
    library(cowplot)
    
    # Read ancestry
    anc <- fread("${ancestry_file}")
    anc[, MajorGroup := floor(AncGroupID / 100)]
    
    major_labels <- c("1" = "AFR", "2" = "MEN", "3" = "EUR", "4" = "SAS",
                     "5" = "EAS", "6" = "AMR", "7" = "OCN", "8" = "MIX")
    anc[, MajorGroupLabel := major_labels[as.character(MajorGroup)]]
    setnames(anc, "Sample", "IID", skip_absent = TRUE)
    
    # Read FAM file
    fam <- fread("${fam}", header = FALSE)
    setnames(fam, c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))
    
    # Merge ancestry with samples
    merged <- merge(fam, anc, by = "IID", all.x = TRUE)
    
    # Calculate missingness by ancestry
    system("plink --bfile ${sample_id}_${server} --missing --out ${sample_id}_${server}_temp")
    
    miss <- fread("${sample_id}_${server}_temp.imiss")
    miss_by_anc <- merge(miss, merged[, .(FID, IID, MajorGroupLabel)], 
                        by = c("FID", "IID"))
    
    # Summarize by ancestry
    summary_stats <- miss_by_anc[!is.na(MajorGroupLabel), .(
        n_samples = .N,
        mean_call_rate = mean(1 - F_MISS, na.rm = TRUE),
        sd_call_rate = sd(1 - F_MISS, na.rm = TRUE),
        min_call_rate = min(1 - F_MISS, na.rm = TRUE),
        max_call_rate = max(1 - F_MISS, na.rm = TRUE)
    ), by = MajorGroupLabel]
    
    # Write summary
    cat("QC Metrics by Ancestry Group\\n", 
        file = "${sample_id}_${server}_qc_by_ancestry.txt")
    cat("============================\\n", 
        file = "${sample_id}_${server}_qc_by_ancestry.txt", append = TRUE)
    write.table(summary_stats, 
               file = "${sample_id}_${server}_qc_by_ancestry.txt",
               append = TRUE, quote = FALSE, row.names = FALSE, sep = "\\t")
    
    # Create plots
    pdf("${sample_id}_${server}_qc_by_ancestry.pdf", width = 14, height = 10)
    
    anc_colors <- c(
        "AFR" = "#E41A1C", "MEN" = "#FF7F00", "EUR" = "#377EB8",
        "SAS" = "#4DAF4A", "EAS" = "#984EA3", "AMR" = "#FFFF33",
        "OCN" = "#A65628", "MIX" = "#999999"
    )
    
    # Sample size by ancestry
    p1 <- ggplot(summary_stats, aes(x = reorder(MajorGroupLabel, -n_samples), 
                                    y = n_samples, fill = MajorGroupLabel)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = n_samples), vjust = -0.5) +
        scale_fill_manual(values = anc_colors) +
        theme_bw(base_size = 12) +
        labs(title = "Sample Size by Ancestry",
             x = "Ancestry Group", y = "Number of Samples") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none")
    
    # Call rate by ancestry
    p2 <- ggplot(miss_by_anc[!is.na(MajorGroupLabel)], 
                aes(x = MajorGroupLabel, y = 1 - F_MISS, fill = MajorGroupLabel)) +
        geom_boxplot() +
        scale_fill_manual(values = anc_colors) +
        theme_bw(base_size = 12) +
        labs(title = "Call Rate Distribution by Ancestry",
             x = "Ancestry Group", y = "Call Rate") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none")
    
    combined <- plot_grid(p1, p2, ncol = 1)
    print(combined)
    
    dev.off()
    
    cat("QC by ancestry analysis complete\\n")
    """
}

/*
 * Process 10: Imputation Performance by Ancestry
 * Analyze imputation quality (INFO scores, R²) stratified by GRAF-anc ancestry
 * NEW PROCESS
 */
process IMPUTATION_PERFORMANCE_BY_ANCESTRY {
    tag "${sample_id}_${server}"
    label 'process_medium'
    publishDir "${params.outdir}/module7/06_imputation_performance/${server}", mode: 'copy'
    
    input:
    tuple val(sample_id), val(server), path(vcf), path(index)
    path(grafanc_results)
    path(grafanc_major_groups)
    path(grafanc_subgroups)
    
    output:
    path("${sample_id}_${server}_imputation_by_ancestry_summary.txt"), emit: summary
    path("${sample_id}_${server}_imputation_by_ancestry_plots.pdf"), emit: plots
    path("${sample_id}_${server}_imputation_metrics_by_major_group.txt"), emit: major_metrics
    path("${sample_id}_${server}_imputation_metrics_by_subgroup.txt"), emit: sub_metrics
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    echo "=== Analyzing imputation performance by GRAF-anc ancestry ==="
    
    # Extract imputation quality metrics from VCF
    # Most imputation servers provide INFO or R2 scores in the INFO field
    bcftools query -f '%CHROM\\t%POS\\t%ID\\t%INFO/INFO\\t%INFO/R2\\n' ${vcf} | \\
        head -100000 | \\
        gzip > ${sample_id}_${server}_imputation_info.txt.gz
    
    # Get sample list from VCF
    bcftools query -l ${vcf} > ${sample_id}_${server}_samples.txt
    
    echo "=== Running R analysis ==="
    
    Rscript - <<'RSCRIPT'
    library(data.table)
    library(ggplot2)
    library(cowplot)
    
    # Read GRAF-anc results
    grafanc <- fread("${grafanc_results}")
    
    # Extract ancestry groups
    grafanc[, MajorGroup := floor(AncGroupID / 100)]
    major_labels <- c(
        "1" = "AFR", "2" = "MEN", "3" = "EUR", "4" = "SAS",
        "5" = "EAS", "6" = "AMR", "7" = "OCN", "8" = "MIX"
    )
    grafanc[, MajorGroupLabel := major_labels[as.character(MajorGroup)]]
    
    # Define subgroup labels (abbreviated)
    subgroup_labels <- c(
        "101" = "Nigeria", "102" = "W_Africa", "103" = "C_Africa",
        "104" = "Kenya", "105" = "S_Africa", "106" = "NE_Africa",
        "107" = "Afr_American", "108" = "Other_Africa",
        "201" = "N_Africa", "202" = "Middle_East_1", "203" = "Middle_East_2",
        "301" = "Finland", "302" = "N_Europe", "303" = "W_Europe",
        "304" = "S_Europe", "305" = "NE_Europe", "306" = "SE_Europe",
        "307" = "Balkans", "308" = "Other_Europe",
        "401" = "Asian_India", "402" = "Gujarati", "403" = "N_S_Asia",
        "404" = "Sri_Lanka", "405" = "Bangladesh",
        "501" = "Japan_Ryukyu", "502" = "Japan", "503" = "Korea",
        "504" = "N_Asia", "505" = "N_China_1", "506" = "N_China_2",
        "507" = "S_China_1", "508" = "S_China_2", "509" = "SE_Asia",
        "510" = "Thailand", "511" = "Other_E_Asia",
        "601" = "Latin_Am_1", "602" = "Latin_Am_2", "603" = "Native_Am",
        "700" = "Oceania",
        "800" = "Multiracial"
    )
    grafanc[, SubgroupLabel := subgroup_labels[as.character(AncGroupID)]]
    grafanc[is.na(SubgroupLabel), SubgroupLabel := paste0("Group_", AncGroupID)]
    
    # Match sample column
    setnames(grafanc, "Sample", "IID", skip_absent = TRUE)
    
    # Read imputation INFO scores
    cat("Reading imputation metrics...\\n")
    info_data <- fread("${sample_id}_${server}_imputation_info.txt.gz")
    setnames(info_data, c("CHROM", "POS", "ID", "INFO", "R2"))
    
    # Handle missing values
    info_data[INFO == ".", INFO := NA]
    info_data[R2 == ".", R2 := NA]
    info_data[, INFO := as.numeric(INFO)]
    info_data[, R2 := as.numeric(R2)]
    
    # Determine which metric to use (INFO or R2)
    has_info <- sum(!is.na(info_data\$INFO)) > 0
    has_r2 <- sum(!is.na(info_data\$R2)) > 0
    
    if (has_info) {
        metric_col <- "INFO"
        metric_name <- "INFO Score"
    } else if (has_r2) {
        metric_col <- "R2"
        metric_name <- "R² Score"
    } else {
        stop("No imputation quality metrics found in VCF")
    }
    
    cat(sprintf("Using %s as quality metric\\n", metric_name))
    
    # Calculate overall statistics
    overall_stats <- data.table(
        Metric = metric_name,
        Mean = mean(info_data[[metric_col]], na.rm = TRUE),
        Median = median(info_data[[metric_col]], na.rm = TRUE),
        SD = sd(info_data[[metric_col]], na.rm = TRUE),
        Q25 = quantile(info_data[[metric_col]], 0.25, na.rm = TRUE),
        Q75 = quantile(info_data[[metric_col]], 0.75, na.rm = TRUE),
        N_variants = sum(!is.na(info_data[[metric_col]]))
    )
    
    # Statistics by chromosome
    chr_stats <- info_data[!is.na(get(metric_col)), .(
        Mean = mean(get(metric_col)),
        Median = median(get(metric_col)),
        SD = sd(get(metric_col)),
        N = .N
    ), by = CHROM]
    setorder(chr_stats, CHROM)
    
    # Read sample list and merge with ancestry
    samples <- fread("${sample_id}_${server}_samples.txt", header = FALSE)
    setnames(samples, "IID")
    
    sample_ancestry <- merge(samples, grafanc[, .(IID, MajorGroupLabel, SubgroupLabel, 
                                                   AncGroupID, MajorGroup)],
                            by = "IID", all.x = TRUE)
    
    # Count samples by ancestry
    major_counts <- sample_ancestry[!is.na(MajorGroupLabel), .N, by = MajorGroupLabel]
    setorder(major_counts, MajorGroupLabel)
    
    sub_counts <- sample_ancestry[!is.na(SubgroupLabel), .N, 
                                 by = .(SubgroupLabel, MajorGroupLabel)]
    setorder(sub_counts, SubgroupLabel)
    
    # ========================================================================
    # WRITE SUMMARY FILE
    # ========================================================================
    cat("Imputation Performance by GRAF-anc Ancestry\\n", 
        file = "${sample_id}_${server}_imputation_by_ancestry_summary.txt")
    cat("============================================\\n", 
        file = "${sample_id}_${server}_imputation_by_ancestry_summary.txt", append = TRUE)
    cat(sprintf("Sample: %s\\n", "${sample_id}"), 
        file = "${sample_id}_${server}_imputation_by_ancestry_summary.txt", append = TRUE)
    cat(sprintf("Server: %s\\n", "${server}"), 
        file = "${sample_id}_${server}_imputation_by_ancestry_summary.txt", append = TRUE)
    cat(sprintf("Quality Metric: %s\\n\\n", metric_name), 
        file = "${sample_id}_${server}_imputation_by_ancestry_summary.txt", append = TRUE)
    
    cat("Overall Imputation Quality:\\n", 
        file = "${sample_id}_${server}_imputation_by_ancestry_summary.txt", append = TRUE)
    write.table(overall_stats, 
               file = "${sample_id}_${server}_imputation_by_ancestry_summary.txt",
               append = TRUE, quote = FALSE, row.names = FALSE, sep = "\\t")
    
    cat("\\n\\nSample Distribution by Major Ancestry Group:\\n", 
        file = "${sample_id}_${server}_imputation_by_ancestry_summary.txt", append = TRUE)
    write.table(major_counts, 
               file = "${sample_id}_${server}_imputation_by_ancestry_summary.txt",
               append = TRUE, quote = FALSE, row.names = FALSE, sep = "\\t")
    
    cat("\\n\\nSample Distribution by Subcontinental Group:\\n", 
        file = "${sample_id}_${server}_imputation_by_ancestry_summary.txt", append = TRUE)
    write.table(sub_counts, 
               file = "${sample_id}_${server}_imputation_by_ancestry_summary.txt",
               append = TRUE, quote = FALSE, row.names = FALSE, sep = "\\t")
    
    cat("\\n\\nImputation Quality by Chromosome:\\n", 
        file = "${sample_id}_${server}_imputation_by_ancestry_summary.txt", append = TRUE)
    write.table(chr_stats, 
               file = "${sample_id}_${server}_imputation_by_ancestry_summary.txt",
               append = TRUE, quote = FALSE, row.names = FALSE, sep = "\\t")
    
    # Save detailed metrics
    fwrite(major_counts, "${sample_id}_${server}_imputation_metrics_by_major_group.txt", 
           sep = "\\t")
    fwrite(sub_counts, "${sample_id}_${server}_imputation_metrics_by_subgroup.txt", 
           sep = "\\t")
    
    # ========================================================================
    # CREATE VISUALIZATIONS
    # ========================================================================
    pdf("${sample_id}_${server}_imputation_by_ancestry_plots.pdf", 
        width = 16, height = 12)
    
    # Plot 1: Distribution of imputation quality scores
    p1 <- ggplot(info_data[!is.na(get(metric_col))], 
                aes(x = get(metric_col))) +
        geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
        geom_vline(xintercept = median(info_data[[metric_col]], na.rm = TRUE),
                  color = "red", linetype = "dashed", size = 1) +
        theme_bw(base_size = 14) +
        labs(title = sprintf("Distribution of Imputation Quality (%s)", metric_name),
             subtitle = sprintf("%s (%s) - %d variants", 
                               "${sample_id}", "${server}", 
                               sum(!is.na(info_data[[metric_col]]))),
             x = metric_name,
             y = "Number of Variants") +
        theme(plot.title = element_text(face = "bold"))
    print(p1)
    
    # Plot 2: Imputation quality by chromosome
    p2 <- ggplot(chr_stats, aes(x = factor(CHROM), y = Mean)) +
        geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
        geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), 
                     width = 0.3, color = "black") +
        theme_bw(base_size = 14) +
        labs(title = sprintf("Mean Imputation Quality by Chromosome (%s)", metric_name),
             subtitle = sprintf("%s (%s)", "${sample_id}", "${server}"),
             x = "Chromosome",
             y = sprintf("Mean %s ± SD", metric_name)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(face = "bold"))
    print(p2)
    
    # Plot 3: Sample distribution by major ancestry group
    anc_colors <- c(
        "AFR" = "#E41A1C", "MEN" = "#FF7F00", "EUR" = "#377EB8",
        "SAS" = "#4DAF4A", "EAS" = "#984EA3", "AMR" = "#FFFF33",
        "OCN" = "#A65628", "MIX" = "#999999"
    )
    
    p3 <- ggplot(major_counts, 
                aes(x = reorder(MajorGroupLabel, -N), y = N, 
                    fill = MajorGroupLabel)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = N), vjust = -0.5, size = 5) +
        scale_fill_manual(values = anc_colors, na.value = "grey50") +
        theme_bw(base_size = 14) +
        labs(title = "Sample Distribution by GRAF-anc Major Ancestry Group",
             subtitle = sprintf("%s (%s)", "${sample_id}", "${server}"),
             x = "Major Ancestry Group",
             y = "Number of Samples") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none",
              plot.title = element_text(face = "bold"))
    print(p3)
    
    # Plot 4: Sample distribution by subcontinental group (top 15)
    top_subgroups <- head(sub_counts[order(-N)], 15)
    
    p4 <- ggplot(top_subgroups, 
                aes(x = reorder(SubgroupLabel, N), y = N, 
                    fill = MajorGroupLabel)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = N), hjust = -0.2, size = 3) +
        scale_fill_manual(values = anc_colors, na.value = "grey50") +
        coord_flip() +
        theme_bw(base_size = 12) +
        labs(title = "Sample Distribution by Subcontinental Ancestry (Top 15)",
             subtitle = sprintf("%s (%s)", "${sample_id}", "${server}"),
             x = "Subcontinental Group",
             y = "Number of Samples",
             fill = "Major\\nGroup") +
        theme(plot.title = element_text(face = "bold"))
    print(p4)
    
    # Plot 5: Boxplot of INFO by chromosome (subset for clarity)
    if (nrow(info_data) > 50000) {
        info_subset <- info_data[sample(.N, 50000)]
    } else {
        info_subset <- info_data
    }
    
    p5 <- ggplot(info_subset[!is.na(get(metric_col))], 
                aes(x = factor(CHROM), y = get(metric_col))) +
        geom_boxplot(fill = "steelblue", alpha = 0.5, outlier.size = 0.5) +
        theme_bw(base_size = 14) +
        labs(title = sprintf("Distribution of %s by Chromosome", metric_name),
             subtitle = sprintf("%s (%s)", "${sample_id}", "${server}"),
             x = "Chromosome",
             y = metric_name) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(face = "bold"))
    print(p5)
    
    dev.off()
    
    cat("\\n=== Imputation performance analysis by ancestry complete ===\\n")
    cat("Note: Imputation quality is per-variant, not per-sample.\\n")
    cat("Sample counts show the distribution of ancestries in your dataset,\\n")
    cat("which is useful for understanding representation when stratifying analyses.\\n")
RSCRIPT
    
    echo "=== Analysis complete ==="
    """
}

// ============================================================================
// WORKFLOW DEFINITION
// ============================================================================

workflow MODULE7_ANCESTRY {
    take:
    qc_data      // From Module 6: (sample_id, server, bed, bim, fam)
    pcs          // From Module 6: (sample_id, server, eigenvec, eigenval)
    imputed_vcf  // From imputation: (sample_id, server, vcf, index)
    
    main:
    // Step 1: Prepare data for ancestry analysis
    PREPARE_FOR_ANCESTRY(qc_data)
    
    // Step 2: GRAF-anc global ancestry estimation (CORRECTED)
    GRAFANC_ANCESTRY(PREPARE_FOR_ANCESTRY.out.pruned_data)
    
    // Step 3: ADMIXTURE analysis for K=2 to K=12 (or user-specified)
    k_min = params.admixture_k_min ?: 2
    k_max = params.admixture_k_max ?: 12
    k_values = Channel.from(k_min..k_max)
    
    admixture_input = PREPARE_FOR_ANCESTRY.out.pruned_data.combine(k_values)
    RUN_ADMIXTURE(admixture_input)
    
    // Step 4: Aggregate ADMIXTURE results (standard plots)
    admixture_grouped = RUN_ADMIXTURE.out.results.groupTuple(by: [0, 1])
    
    SUMMARIZE_ADMIXTURE(
        admixture_grouped,
        RUN_ADMIXTURE.out.cv_error.collect(),
        PREPARE_FOR_ANCESTRY.out.pruned_data.map { it[4] }
    )
    
    // Step 4B: ADMIXTURE organized by GRAF-anc ancestry (K=6, K=9) - NEW
    ADMIXTURE_BY_GRAFANC(
        SUMMARIZE_ADMIXTURE.out.for_grafanc_org,
        GRAFANC_ANCESTRY.out.ancestry
    )
    
    // Step 5: Local Ancestry Inference (OPTIONAL - if enabled)
    if (params.run_lai) {
        // Get VCF data for LAI
        lai_input = imputed_vcf
        
        PREPARE_FOR_LAI(lai_input)
        
        // Split by chromosome
        chromosomes = Channel.of(1..22)
        lai_chr_data = PREPARE_FOR_LAI.out.phased_vcf.combine(chromosomes)
        
        // Run selected LAI tool
        if (params.lai_tool == 'rfmix_v2' || params.lai_tool == 'all') {
            reference_vcf = file(params.rfmix_reference_vcf)
            reference_map = file(params.rfmix_reference_map)
            genetic_maps = file(params.genetic_maps_dir)
            
            RFMIX_V2(
                lai_chr_data,
                reference_vcf,
                reference_map,
                genetic_maps
            )
            
            // Summarize LAI results
            lai_grouped = RFMIX_V2.out.results.groupTuple(by: [0, 1])
            
            SUMMARIZE_LAI(
                lai_grouped,
                PREPARE_FOR_ANCESTRY.out.pruned_data.map { it[4] }
            )
        }
    }
    
    // Step 6: PCA colored by ancestry
    pca_with_ancestry = pcs.join(GRAFANC_ANCESTRY.out.ancestry, by: [0, 1])
    PCA_BY_ANCESTRY(pca_with_ancestry)
    
    // Step 7: QC metrics by ancestry group
    qc_by_ancestry_input = PREPARE_FOR_ANCESTRY.out.original_data
        .join(GRAFANC_ANCESTRY.out.ancestry, by: [0, 1])
    QC_BY_ANCESTRY(qc_by_ancestry_input)
    
    // Step 8: Imputation performance by ancestry (NEW)
    imputation_ancestry_input = imputed_vcf
        .combine(GRAFANC_ANCESTRY.out.ancestry, by: [0, 1])
        .combine(GRAFANC_ANCESTRY.out.major_groups, by: [0, 1])
        .combine(GRAFANC_ANCESTRY.out.subgroups, by: [0, 1])
    
    IMPUTATION_PERFORMANCE_BY_ANCESTRY(imputation_ancestry_input)
    
    emit:
    // GRAF-anc outputs
    grafanc_ancestry = GRAFANC_ANCESTRY.out.ancestry
    grafanc_plots = GRAFANC_ANCESTRY.out.plots
    grafanc_summary = GRAFANC_ANCESTRY.out.summary
    grafanc_major_groups = GRAFANC_ANCESTRY.out.major_groups
    grafanc_subgroups = GRAFANC_ANCESTRY.out.subgroups
    
    // ADMIXTURE outputs
    admixture_summary = SUMMARIZE_ADMIXTURE.out.summary
    admixture_cv_plot = SUMMARIZE_ADMIXTURE.out.cv_plot
    admixture_barplots = SUMMARIZE_ADMIXTURE.out.barplots
    admixture_by_ancestry_plots = ADMIXTURE_BY_GRAFANC.out.plots
    admixture_by_ancestry_summary = ADMIXTURE_BY_GRAFANC.out.summary
    
    // LAI outputs (if enabled)
    lai_summary = params.run_lai ? SUMMARIZE_LAI.out.summary : Channel.empty()
    
    // Integration outputs
    pca_colored = PCA_BY_ANCESTRY.out.plot
    qc_by_ancestry = QC_BY_ANCESTRY.out.summary
    imputation_performance_summary = IMPUTATION_PERFORMANCE_BY_ANCESTRY.out.summary
    imputation_performance_plots = IMPUTATION_PERFORMANCE_BY_ANCESTRY.out.plots
    imputation_major_metrics = IMPUTATION_PERFORMANCE_BY_ANCESTRY.out.major_metrics
    imputation_sub_metrics = IMPUTATION_PERFORMANCE_BY_ANCESTRY.out.sub_metrics
}

// ============================================================================
// WORKFLOW END
// ============================================================================

