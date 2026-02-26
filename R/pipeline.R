#' Run Complete RNA-seq Analysis Pipeline
#'
#' Runs the complete RNA-seq analysis pipeline including preprocessing,
#' differential expression, pathway analysis, WGCNA, and visualization.
#'
#' @param counts_file Path to counts file (TSV format)
#' @param tpm_file Path to TPM file (TSV format)
#' @param metadata_file Path to metadata file (TSV format)
#' @param gtf_file Path to GTF annotation file
#' @param group1_condition Name of treatment/group1 condition
#' @param group2_condition Name of control/group2 condition
#' @param condition_column Column name in metadata for conditions (default: "condition")
#' @param sample_id_column Column name in metadata for sample IDs (default: "SampleID")
#' @param experiment_name Name for the experiment
#' @param output_dir Output directory
#' @param species Species ("mouse" or "human", default: "mouse")
#' @param run_wgcna Whether to run WGCNA analysis (default: TRUE)
#' @param run_tf_analysis Whether to run transcription factor analysis (default: TRUE)
#' @param run_immune_analysis Whether to run immune deconvolution (default: TRUE)
#' @param run_lincs_analysis Whether to run LINCS connectivity analysis (default: TRUE)
#' @param remove_samples Character vector of sample names to remove (optional)
#' @return List containing results from all analyses
#' @export
run_complete_pipeline <- function(counts_file, tpm_file, metadata_file, gtf_file,
                                 group1_condition, group2_condition,
                                 condition_column = "condition",
                                 sample_id_column = "SampleID",
                                 experiment_name, output_dir, species = "mouse",
                                 run_wgcna = TRUE, run_tf_analysis = TRUE,
                                 run_immune_analysis = TRUE, run_lincs_analysis = TRUE,
                                 remove_samples = NULL, lfc_threshold = 1, 
                                 color_volcano_up="#CA3433", color_volcano_down="#2B7CB6", genes_to_label=NULL, 
                                 point_size_volcano=4,label_size_volcano=5,n_labels_up=10,n_labels_down=10,
                                 covariates = NULL, design_formula = NULL, contrast_string = NULL,
                                 fdr_threshold = 0.05, gsea_custom_pathways=NULL, n_gsea_enrich_up=5, n_gsea_enrich_down=5, 
                                 color_gsea_down="#2B7CB6", color_gsea_up="#CA3433", color_gsea_ns="#C5C6C7",
                                 run_ssgsea=TRUE, ssgsea_extra_pathways=NULL, ssgsea_n_boxplot_pathways = 20,
                                 customed_pathways = NULL)
  {
  # !!! Added customed_pathways arg !!!
  
  # Start timing
  start_time <- Sys.time()


  message("Starting Complete RNA-seq Analysis Pipeline")

  message("Experiment: ", experiment_name)
  message("Output directory: ", output_dir)
  message("Species: ", species)
  message("Comparison: ", group1_condition, " vs ", group2_condition)
  message("")

  # Create output directory
  create_output_dir(output_dir)

  # Initialize results list
  results <- list()

  # ========== 1. Load and Preprocess Data ==========
  message("Step 1: Loading and preprocessing data...")

  # Load data
  gtf <- rtracklayer::import(gtf_file)
  counts <- read.table(counts_file, header = TRUE, sep = "\t")
  colnames(counts)<-c(colnames(counts)[1:2],clean_sample_names(colnames(counts)[3:length(colnames(counts))]))

  tpm <- read.table(tpm_file, header = TRUE, sep = "\t")
  colnames(tpm)<-c(colnames(tpm)[1:2],clean_sample_names(colnames(tpm)[3:length(colnames(tpm))]))
  metadata <- read.delim(metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  message("Loaded data:")
  message("  - Counts: ", nrow(counts), " genes x ", ncol(counts)-2, " samples")
  message("  - TPM: ", nrow(tpm), " genes x ", ncol(tpm)-2, " samples")
  message("  - Metadata: ", nrow(metadata), " samples")
  print(metadata)
  # Filter protein-coding genes
  pc_counts <- filter_protein_coding_genes(counts, gtf)
  pc_tpm <- filter_protein_coding_genes(tpm, gtf)

  # Prepare expression data
  pc_tpm_processed <- prepare_expression_data(
    pc_tpm,
    log_transform = TRUE,
    remove_samples = remove_samples
  )

  pc_counts_processed <- prepare_expression_data(
    pc_counts,
    log_transform = FALSE,
    remove_samples = remove_samples
  )

  # Match metadata
  metadata_matched <- match_metadata_to_expression(pc_counts_processed, metadata, sample_id_column)
  results$preprocessing <- list(
    pc_counts = pc_counts_processed,
    pc_tpm = pc_tpm_processed,
    metadata = metadata_matched
  )
  # ========== 2. Exploratory Data Analysis ==========
  message("\nStep 2: Exploratory data analysis...")
  # Define shared color mapping
  conditions <- unique(metadata_matched$condition)
  n_conditions <- length(conditions)
  # !!! Changed way palette is taken !!!
  pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9."Set1"))
  condition_colors <- setNames(
    pal(n_levels),
    levels_color
  )
  condition_colors <- setNames(
    RColorBrewer::brewer.pal(n_conditions, "Set1"),
    conditions
  )
  
  # PCA plot
  pca_result <- create_pca_plot(
    pc_tpm_processed, metadata_matched,
    color_by = condition_column,
    output_file = file.path(output_dir, paste0(experiment_name, "_PCA.pdf"))
  )

  # Expression heatmap
  color_map_list <- list(condition = condition_colors)

  create_expression_heatmap(
    pc_tpm_processed, metadata_matched,
    annotation_columns = condition_column,
    output_file = file.path(output_dir, paste0(experiment_name, "_expression_heatmap.pdf"))
  )

  results$exploratory <- list(
    pca = pca_result
  )

  # ========== 3. Differential Expression Analysis ==========
  message("\nStep 3: Differential expression analysis...")

  de_results <- run_differential_expression(
    counts_data = pc_counts_processed,
    metadata = metadata_matched,
    group1_condition = group1_condition,
    group2_condition = group2_condition,
    condition_column = condition_column,
    sample_id_column = sample_id_column,
    experiment_name = experiment_name,
    output_dir = output_dir,
    lfc_threshold = lfc_threshold,
    fdr_threshold = fdr_threshold,
    covariates = covariates,
    design_formula = design_formula,
    contrast_string = contrast_string
  )

  # Create visualizations
  create_volcano_plot(
    de_results$de_results$efit,
    fdr_threshold = fdr_threshold,              # keep your current choice
    lfc_threshold = lfc_threshold,
    output_file  = file.path(output_dir, paste0(experiment_name, "_volcano_plot.pdf")),
    color_up     = color_volcano_up,
    color_down   = color_volcano_down,
    point_size        = point_size_volcano,
    label_size        = label_size_volcano,
    n_labels_up       = n_labels_up,
    n_labels_down     = n_labels_down,
    highlight_genes = genes_to_label
    # color_ns keeps default "grey50", or you can add a color_volcano_ns arg as well
  )

  create_ma_plot(
    de_results$de_results$efit,
    output_file = file.path(output_dir, paste0(experiment_name, "_MA_plot.pdf"))
  )

  create_pie_chart(
    de_results$de_results$efit,
    output_file = file.path(output_dir, paste0(experiment_name, "_DE_pie_chart.pdf"))
  )

  results$differential_expression <- de_results

  # ========== 4. Pathway Analysis ==========
  message("\nStep 4: Pathway analysis...")

  # GSEA analysis
  gsea_results <- run_gsea_analysis(
    de_results = de_results$de_results$efit,
    species = ifelse(species == "mouse", "MM", "HS"),
    customed_pathways = customed_pathways
  )

  # Extract gene sets and ranks once here
  gsea_gene_sets  <- attr(gsea_results, "gene_sets")
  gsea_gene_ranks <- attr(gsea_results, "gene_ranks")
  print(gsea_gene_sets)
  print(gsea_gene_ranks)
  # Create GSEA visualizations
  plot_gsea_barplot(
    gsea_results,
    output_file = file.path(output_dir, paste0(experiment_name, "_GSEA_barplot.pdf")),
    color_up = color_gsea_up,
    color_down = color_gsea_down,
    color_ns=color_gsea_ns
  )

  plot_gsea_dotplot(
    gsea_results,
    output_file = file.path(output_dir, paste0(experiment_name, "_GSEA_dotplot.pdf")),
    color_up = color_gsea_up,
    color_down = color_gsea_down,
    color_ns = color_gsea_ns
  )
  create_gsea_table_plot(
    gsea_results,
    output_file = file.path(output_dir, paste0(experiment_name, "_GSEA_tableplot.pdf")),
    gene_sets = gsea_gene_sets,
    gene_ranks = gsea_gene_ranks
  )

  plot_gsea_enrichment(
    gsea_results,
    pathways    = gsea_custom_pathways,   # NULL ⇒ auto top up/down
    n_up        = n_gsea_enrich_up,
    n_down      = n_gsea_enrich_down,
    output_file = file.path(output_dir, paste0(experiment_name, "_GSEA_enrichment.pdf")),
    gene_sets   = gsea_gene_sets,
    gene_ranks  = gsea_gene_ranks
  )

  # Save pathway results
  save_pathway_results(gsea_results, NULL, experiment_name, output_dir)

  results$pathway_analysis <- list(
    gsea = gsea_results,
    gene_sets = gsea_gene_sets,
    gene_ranks = gsea_gene_ranks
  )

  if (run_ssgsea) {
    message("\nStep 5a: ssGSEA pathway activity analysis...")

    ssgsea_results <- run_ssgsea_analysis(
      expr_data        = pc_tpm_processed,
      metadata         = metadata_matched,
      condition_column = condition_column,
      sample_id_column = sample_id_column,
      group1           = group1_condition,
      group2           = group2_condition,
      species          = species,
      extra_pathways   = ssgsea_extra_pathways,
      output_dir       = output_dir,
      experiment_name  = experiment_name,
      n_boxplot_pathways = ssgsea_n_boxplot_pathways
    )

    results$ssgsea <- ssgsea_results
  }


  if (run_wgcna) {
    message("\nStep 5b: WGCNA enriched co-expression analysis...")

    enhanced_results <- run_enhanced_wgcna_analysis(
    expr_data = pc_tpm_processed,
    metadata = metadata_matched,
    trait_column = condition_column,
    output_dir = output_dir,
    experiment_name = experiment_name,
    organism = species, # Explicitly set to mouse
    species_id = 10090 # Explicitly set the mouse ID for STRING-db
  )
    results$wgcna_enhanced <- enhanced_results
  }

  # ========== 6. Transcription Factor Analysis ==========
  if (run_tf_analysis) {
    message("\nStep 6: Transcription factor activity analysis...")

    tf_results <- run_tf_analysis(
      expr_data = pc_tpm_processed,
      metadata = metadata_matched,
      condition_column = condition_column,
      group1 = group1_condition,
      group2 = group2_condition,
      species = species,
      output_dir = output_dir,
      experiment_name = experiment_name,
      fdr_threshold = 0.05, lfc_threshold = 1
    )

    results$transcription_factors <- tf_results
  }




  # ========== 7. Immune Cell Deconvolution ==========
  if (run_immune_analysis) {
    message("\nStep 7: Immune cell deconvolution analysis...")

    immune_results <- run_immune_deconvolution(
      expr_data = pc_tpm_processed,
      metadata = metadata_matched,
      condition_column = condition_column,
      species = species,
      output_dir = output_dir,
      experiment_name = experiment_name
    )

    results$immune_analysis <- immune_results
  }

  # ========== 8. LINCS Connectivity Analysis ==========
  if (run_lincs_analysis) {
    message("\nStep 8: LINCS connectivity analysis...")

    lincs_results <- run_lincs_analysis(
      de_results = de_results$de_results$efit,
      species = species,
      output_dir = output_dir,
      experiment_name = experiment_name
    )

    results$lincs <- lincs_results
  }

  # ========== 9. Create Summary Report ==========
  message("\nStep 9: Creating summary report...")

  create_analysis_summary_report(results, experiment_name, output_dir,
                                group1_condition, group2_condition)

  # Create comprehensive summary plot
  create_summary_plot(
    expr_data = pc_tpm_processed,
    de_results = de_results$de_results$efit,
    metadata = metadata_matched,
    condition_column = condition_column,
    output_file = file.path(output_dir, paste0(experiment_name, "_summary_plot.pdf"))
  )

  # Calculate total runtime
  end_time <- Sys.time()
  total_time <- end_time - start_time


  message("RNA-seq Analysis Pipeline Completed Successfully!")

  message("Total runtime: ", round(total_time, 2), " ", attr(total_time, "units"))
  message("Results saved in: ", output_dir)


  return(results)
}

#' Run Streamlined Pipeline
#'
#' Runs a streamlined version of the pipeline with essential analyses only.
#'
#' @param counts_file Path to counts file
#' @param tpm_file Path to TPM file
#' @param metadata_file Path to metadata file
#' @param gtf_file Path to GTF file
#' @param group1_condition Treatment condition name
#' @param group2_condition Control condition name
#' @param experiment_name Experiment name
#' @param output_dir Output directory
#' @param species Species (default: "mouse")
#' @return List with essential results
#' @export
run_streamlined_pipeline <- function(counts_file, tpm_file, metadata_file, gtf_file,
                                   group1_condition, group2_condition,
                                   experiment_name, output_dir, species = "mouse") {

  message("Running streamlined RNA-seq analysis pipeline...")

  # Run essential analyses only
  results <- run_complete_pipeline(
    counts_file = counts_file,
    tpm_file = tpm_file,
    metadata_file = metadata_file,
    gtf_file = gtf_file,
    group1_condition = group1_condition,
    group2_condition = group2_condition,
    experiment_name = experiment_name,
    output_dir = output_dir,
    species = species,
    run_wgcna = FALSE,
    run_tf_analysis = FALSE,
    run_immune_analysis = FALSE,
    run_lincs_analysis = FALSE
  )

  message("Streamlined analysis completed!")
  return(results)
}

#' Create Analysis Summary Report
#'
#' Creates a comprehensive text summary report of all analyses.
#'
#' @param results Results list from pipeline
#' @param experiment_name Experiment name
#' @param output_dir Output directory
#' @param group1_condition Group 1 condition name
#' @param group2_condition Group 2 condition name
#' @return None (creates report file)
#' @keywords internal
create_analysis_summary_report <- function(results, experiment_name, output_dir,
                                         group1_condition, group2_condition) {

  report_file <- file.path(output_dir, paste0(experiment_name, "_ANALYSIS_SUMMARY_REPORT.txt"))

  # Create report header
  cat("RNA-seq Analysis Summary Report\n", file = report_file)
  cat("===============================\n\n", file = report_file, append = TRUE)
  cat("Experiment:", experiment_name, "\n", file = report_file, append = TRUE)
  cat("Analysis Date:", Sys.Date(), "\n", file = report_file, append = TRUE)
  cat("Comparison:", group1_condition, "vs", group2_condition, "\n\n", file = report_file, append = TRUE)

  # Data preprocessing summary
  cat("DATA PREPROCESSING:\n", file = report_file, append = TRUE)
  cat("-------------------\n", file = report_file, append = TRUE)
  if (!is.null(results$preprocessing)) {
    cat("Protein-coding genes:", nrow(results$preprocessing$pc_counts), "\n",
        file = report_file, append = TRUE)
    cat("Samples analyzed:", ncol(results$preprocessing$pc_counts), "\n\n",
        file = report_file, append = TRUE)
  }

  # Differential expression summary
  cat("DIFFERENTIAL EXPRESSION:\n", file = report_file, append = TRUE)
  cat("------------------------\n", file = report_file, append = TRUE)
  if (!is.null(results$differential_expression$summary)) {
    summary_de <- results$differential_expression$summary
    cat("Total genes tested:", summary_de$total_genes, "\n", file = report_file, append = TRUE)
    cat("Significant genes:", summary_de$significant_genes, "\n", file = report_file, append = TRUE)
    cat("Upregulated:", summary_de$upregulated, "\n", file = report_file, append = TRUE)
    cat("Downregulated:", summary_de$downregulated, "\n", file = report_file, append = TRUE)
    cat("Percent significant:", summary_de$percent_significant, "%\n\n", file = report_file, append = TRUE)
  }

  # Pathway analysis summary
  cat("PATHWAY ANALYSIS:\n", file = report_file, append = TRUE)
  cat("-----------------\n", file = report_file, append = TRUE)
  if (!is.null(results$pathway_analysis$gsea)) {
    gsea_res <- results$pathway_analysis$gsea
    sig_pathways <- sum(gsea_res$padj < 0.05, na.rm = TRUE)
    cat("Total pathways tested:", nrow(gsea_res), "\n", file = report_file, append = TRUE)
    cat("Significant pathways (padj < 0.05):", sig_pathways, "\n", file = report_file, append = TRUE)

    # Top 5 enriched pathways
    top_pathways <- head(gsea_res[order(-gsea_res$NES), ], 5)
    cat("Top 5 enriched pathways:\n", file = report_file, append = TRUE)
    for (i in 1:nrow(top_pathways)) {
      cat("  ", i, ". ", gsub("HALLMARK_", "", top_pathways$pathway[i]),
          " (NES: ", round(top_pathways$NES[i], 2), ")\n",
          file = report_file, append = TRUE)
    }
    cat("\n", file = report_file, append = TRUE)
  }

  # WGCNA summary
  if (!is.null(results$wgcna)) {
    cat("WGCNA CO-EXPRESSION ANALYSIS:\n", file = report_file, append = TRUE)
    cat("-----------------------------\n", file = report_file, append = TRUE)
    n_modules <- length(unique(results$wgcna$moduleColors))
    cat("Number of modules identified:", n_modules, "\n", file = report_file, append = TRUE)

    # Module sizes
    module_sizes <- table(results$wgcna$moduleColors)
    largest_module <- names(module_sizes)[which.max(module_sizes)]
    cat("Largest module:", largest_module, "(", max(module_sizes), "genes )\n",
        file = report_file, append = TRUE)
    cat("\n", file = report_file, append = TRUE)
  }

  # Transcription factor analysis summary
  if (!is.null(results$transcription_factors)) {
    cat("TRANSCRIPTION FACTOR ANALYSIS:\n", file = report_file, append = TRUE)
    cat("-------------------------------\n", file = report_file, append = TRUE)
    tf_diff <- results$transcription_factors$differential_results
    n_activated <- sum(tf_diff$Regulation == "Activated", na.rm = TRUE)
    n_repressed <- sum(tf_diff$Regulation == "Repressed", na.rm = TRUE)

    cat("Total TFs analyzed:", nrow(tf_diff), "\n", file = report_file, append = TRUE)
    cat("Activated TFs:", n_activated, "\n", file = report_file, append = TRUE)
    cat("Repressed TFs:", n_repressed, "\n", file = report_file, append = TRUE)

    # Top activated TFs
    if (n_activated > 0) {
      top_activated <- head(tf_diff[tf_diff$Regulation == "Activated", ], 3)
      cat("Top activated TFs:", paste(top_activated$TF, collapse = ", "), "\n",
          file = report_file, append = TRUE)
    }

    if (n_repressed > 0) {
      top_repressed <- head(tf_diff[tf_diff$Regulation == "Repressed", ], 3)
      cat("Top repressed TFs:", paste(top_repressed$TF, collapse = ", "), "\n",
          file = report_file, append = TRUE)
    }
    cat("\n", file = report_file, append = TRUE)
  }

  # Immune analysis summary
  if (!is.null(results$immune_analysis)) {
    cat("IMMUNE CELL DECONVOLUTION:\n", file = report_file, append = TRUE)
    cat("---------------------------\n", file = report_file, append = TRUE)

    if (!is.null(results$immune_analysis$statistical_results)) {
      stats_results <- results$immune_analysis$statistical_results
      sig_immune <- stats_results[stats_results$Test_Method == "Wilcoxon" &
                                 stats_results$Adjusted_P_Value < 0.05, ]

      cat("Cell types with significant differences:", nrow(sig_immune), "\n",
          file = report_file, append = TRUE)

      if (nrow(sig_immune) > 0) {
        cat("Significantly different cell types:\n", file = report_file, append = TRUE)
        for (i in 1:min(5, nrow(sig_immune))) {
          cat("  ", sig_immune$Cell_Type[i], " (adj.p = ",
              signif(sig_immune$Adjusted_P_Value[i], 3), ")\n",
              file = report_file, append = TRUE)
        }
      }
    }
    cat("\n", file = report_file, append = TRUE)
  }

  # LINCS analysis summary
  if (!is.null(results$lincs)) {
    cat("LINCS CONNECTIVITY ANALYSIS:\n", file = report_file, append = TRUE)
    cat("-----------------------------\n", file = report_file, append = TRUE)

    if (!is.null(results$lincs$results_table)) {
      lincs_table <- results$lincs$results_table
      score_col <- ifelse("NCS" %in% colnames(lincs_table), "NCS", "score")

      cat("Total drug connections found:", nrow(lincs_table), "\n",
          file = report_file, append = TRUE)
      cat("Drugs with opposite effects:", sum(lincs_table[[score_col]] < -0.1), "\n",
          file = report_file, append = TRUE)
      cat("Drugs with similar effects:", sum(lincs_table[[score_col]] > 0.1), "\n",
          file = report_file, append = TRUE)

      # Top 3 opposite effect drugs (potential therapeutics)
      opposite_drugs <- head(lincs_table[order(lincs_table[[score_col]]), ], 3)
      if ("pert_iname" %in% colnames(opposite_drugs)) {
        cat("Top potential therapeutic compounds:\n", file = report_file, append = TRUE)
        for (i in 1:nrow(opposite_drugs)) {
          cat("  ", i, ". ", opposite_drugs$pert_iname[i],
              " (Score: ", round(opposite_drugs[[score_col]][i], 3), ")\n",
              file = report_file, append = TRUE)
        }
      }
    }
    cat("\n", file = report_file, append = TRUE)
  }

  # File outputs summary
  cat("OUTPUT FILES:\n", file = report_file, append = TRUE)
  cat("-------------\n", file = report_file, append = TRUE)
  cat("All analysis results, plots, and data files have been saved to:\n",
      file = report_file, append = TRUE)
  cat(output_dir, "\n\n", file = report_file, append = TRUE)

  cat("Key output files include:\n", file = report_file, append = TRUE)
  cat("- ", experiment_name, "_DE_results.tsv (Differential expression results)\n",
      file = report_file, append = TRUE)
  cat("- ", experiment_name, "_GSEA_results.tsv (Pathway analysis results)\n",
      file = report_file, append = TRUE)
  cat("- ", experiment_name, "_volcano_plot.pdf (Volcano plot)\n",
      file = report_file, append = TRUE)
  cat("- ", experiment_name, "_summary_plot.pdf (Summary visualization)\n",
      file = report_file, append = TRUE)

  if (!is.null(results$wgcna)) {
    cat("- ", experiment_name, "_module_assignments.tsv (WGCNA modules)\n",
        file = report_file, append = TRUE)
  }

  if (!is.null(results$transcription_factors)) {
    cat("- ", experiment_name, "_differential_TF_results.tsv (TF activity results)\n",
        file = report_file, append = TRUE)
  }

  if (!is.null(results$immune_analysis)) {
    cat("- ", experiment_name, "_mMCPcounter_scores.tsv (Immune cell scores)\n",
        file = report_file, append = TRUE)
  }

  if (!is.null(results$lincs)) {
    cat("- ", experiment_name, "_LINCS_results.tsv (Drug connectivity results)\n",
        file = report_file, append = TRUE)
  }

  cat("\nAnalysis completed successfully!\n", file = report_file, append = TRUE)
  cat("For questions or issues, please refer to the package documentation.\n",
      file = report_file, append = TRUE)

  message("Analysis summary report saved to: ", report_file)
}

#' Validate Pipeline Inputs
#'
#' Validates that all required input files exist and have the correct format.
#'
#' @param counts_file Path to counts file
#' @param tpm_file Path to TPM file
#' @param metadata_file Path to metadata file
#' @param gtf_file Path to GTF file
#' @param group1_condition Group 1 condition name
#' @param group2_condition Group 2 condition name
#' @param condition_column Condition column name
#' @param sample_id_column Sample ID column name
#' @return Logical indicating if all inputs are valid
#' @export
validate_pipeline_inputs <- function(counts_file, tpm_file, metadata_file, gtf_file,
                                    group1_condition, group2_condition,
                                    condition_column = "condition",
                                    sample_id_column = "SampleID") {

  message("Validating pipeline inputs...")

  # Check file existence
  files_to_check <- c(counts_file, tpm_file, metadata_file, gtf_file)
  file_names <- c("counts", "TPM", "metadata", "GTF")

  for (i in seq_along(files_to_check)) {
    if (!file.exists(files_to_check[i])) {
      stop("File not found: ", files_to_check[i], " (", file_names[i], " file)")
    }
    message("Found ", file_names[i], " file: ", files_to_check[i])
  }

  # Check metadata format
  tryCatch({
    metadata <- read.table(metadata_file, header = TRUE, sep = "\t")

    # Check required columns
    if (!condition_column %in% colnames(metadata)) {
      stop("Condition column '", condition_column, "' not found in metadata")
    }

    if (!sample_id_column %in% colnames(metadata)) {
      stop("Sample ID column '", sample_id_column, "' not found in metadata")
    }

    # Check conditions exist
    available_conditions <- unique(metadata[[condition_column]])
    if (!group1_condition %in% available_conditions) {
      stop("Group 1 condition '", group1_condition, "' not found in metadata. Available: ",
           paste(available_conditions, collapse = ", "))
    }

    if (!group2_condition %in% available_conditions) {
      stop("Group 2 condition '", group2_condition, "' not found in metadata. Available: ",
           paste(available_conditions, collapse = ", "))
    }

    message("Metadata validation passed")
    message("  - Samples: ", nrow(metadata))
    message("  - Conditions: ", paste(available_conditions, collapse = ", "))
    message("  - Comparison: ", group1_condition, " (n=",
            sum(metadata[[condition_column]] == group1_condition), ") vs ",
            group2_condition, " (n=",
            sum(metadata[[condition_column]] == group2_condition), ")")

  }, error = function(e) {
    stop("Error reading metadata file: ", e$message)
  })

  # Quick check of expression files format
  tryCatch({
    counts_sample <- read.table(counts_file, header = TRUE, sep = "\t", nrows = 5)
    if (!"gene_id" %in% colnames(counts_sample)) {
      warning("'gene_id' column not found in counts file. Expected format may be incorrect.")
    } else {
      message("Counts file format appears correct")
    }

    tpm_sample <- read.table(tpm_file, header = TRUE, sep = "\t", nrows = 5)
    if (!"gene_id" %in% colnames(tpm_sample)) {
      warning("'gene_id' column not found in TPM file. Expected format may be incorrect.")
    } else {
      message("TPM file format appears correct")
    }

  }, error = function(e) {
    warning("Could not validate expression file formats: ", e$message)
  })

  message("All input validation checks passed!")
  return(TRUE)
}

#' Get Pipeline Configuration
#'
#' Returns the default configuration settings for the pipeline.
#'
#' @return List with default pipeline settings
#' @export
get_pipeline_config <- function() {

  config <- list(
    # Differential expression settings
    fdr_threshold = 0.05,
    lfc_threshold = 1,

    # GSEA settings
    gsea_min_size = 15,
    gsea_max_size = 500,
    gsea_n_perm = 100000,

    # WGCNA settings
    wgcna_n_genes = 3000,
    wgcna_min_module_size = 30,
    wgcna_soft_power = NULL,  # Auto-detect

    # Transcription factor settings
    tf_confidence_levels = c("A", "B", "C"),

    # Immune deconvolution settings
    immune_method = "mMCPcounter",

    # LINCS settings
    lincs_n_up = 150,
    lincs_n_down = 150,

    # Visualization settings
    plot_width = 10,
    plot_height = 8,
    heatmap_n_genes = 500,

    # General settings
    n_threads = 2,
    seed = 123
  )

  return(config)
}

#' Update Pipeline Configuration
#'
#' Updates pipeline configuration with user-specified settings.
#'
#' @param config_updates List of configuration updates
#' @return Updated configuration list
#' @export
update_pipeline_config <- function(config_updates) {

  default_config <- get_pipeline_config()

  # Update with user settings
  for (setting in names(config_updates)) {
    if (setting %in% names(default_config)) {
      default_config[[setting]] <- config_updates[[setting]]
      message("Updated ", setting, " to: ", config_updates[[setting]])
    } else {
      warning("Unknown configuration setting: ", setting)
    }
  }

  return(default_config)
}
