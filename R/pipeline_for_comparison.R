search_mate <- function(vector, mate, pattern, split) {
  vector <- grep(pattern = pattern, unique(vector), value = TRUE)
  mate <- strsplit(mate, split = split)[[1]]
  matches <- sapply(vector, function(x) {
   one <- strsplit(vector, split = split)
   any(mapply(´==´, one, mate, SIMPLIFY = TRUE))
  }
                    
  return(vector[matches])
}


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
                                 stratify_by = NULL, design_formula = NULL, contrast_string = NULL,
                                 experiment_name, output_dir, species = "mouse",
                                 run_wgcna = TRUE, run_tf_analysis = TRUE, run_ssgsea=TRUE,
                                 run_immune_analysis = TRUE, run_lincs_analysis = TRUE,
                                 remove_samples = NULL, lfc_threshold = 1,
                                 color_volcano_up="#CA3433", color_volcano_down="#2B7CB6", genes_to_label=NULL,
                                 fdr_threshold = 0.05, point_size_volcano=4, label_size_volcano=5, n_labels_up=10, n_labels_down=10,
                                 gsea_custom_pathways=NULL, n_gsea_enrich_up=5, n_gsea_enrich_down=5, 
                                 color_gsea_down="#2B7CB6", color_gsea_up="#CA3433", color_gsea_ns="#C5C6C7", 
                                 ssgsea_extra_pathways=NULL, ssgsea_n_boxplot_pathways = 20, customed_pathways)
  {

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

  metadata_matched$merged_col <- apply(metadata_matched[ , c(condition,stratify_by), drop = FALSE], 1, paste, collapse = "--")
  
  # ========== 2. Exploratory Data Analysis ==========
  message("\nStep 2: Exploratory data analysis...")
  # Define shared color mapping
  conditions <- unique(metadata_matched$merged_col)
  n_conditions <- length(conditions)
  # !!! Changed way palette is taken !!!
  pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9."Set1"))
  condition_colors <- setNames(
    pal(n_levels),
    levels_color
  )
  
  # PCA plot
  pca_result <- create_pca_plot(
    pc_tpm_processed, metadata_matched,
    color_by = "merged_col",
    output_file = file.path(output_dir, paste0(experiment_name, "_PCA.pdf"))
  )

  # Expression heatmap
  color_map_list <- list(condition = condition_colors)

  create_expression_heatmap(
    pc_tpm_processed, metadata_matched,
    annotation_columns = "merged_col",
    output_file = file.path(output_dir, paste0(experiment_name, "_expression_heatmap.pdf"))
  )

  results$exploratory <- list(
    pca = pca_result
  )

  # ========== 3. Differential Expression Analysis ==========
  message("\nStep 3: Differential expression analysis...")

  de_results <- list()
  
  for (group in grep(pattern = group1_condition, unique(metadata_matched$merged_col), value = TRUE)) {
    samples <- c(metadata_matched$SampleID[metadata_matched$merged_col == search_mate(metadata_matched$merged_col, 
                                                                                     mate = group, pattern = group2_condition, 
                                                                                     split = "--")],
                metadata_matched$SampleID[metadata_matched$merged_col == group])

    de_results[[group]] <- run_differential_expression(
    counts_data = pc_counts_processed[, samples],
    metadata = metadata_matched[metadata_matched %in% samples, ],
    group1_condition = group1_condition,
    group2_condition = group2_condition,
    condition_column = condition_column,
    sample_id_column = sample_id_column,
    experiment_name = group,
    output_dir = file.path(output_dir, group),
    lfc_threshold = lfc_threshold,
    fdr_threshold = fdr_threshold,
    covariates = covariates,
    design_formula = design_formula,
    contrast_string = contrast_string
    )

    # Create visualizations
    create_volcano_plot(
      de_results[[group]]$de_results$efit,
      fdr_threshold = fdr_threshold,              # keep your current choice
      lfc_threshold = lfc_threshold,
      output_file  = file.path(output_dir, group, paste0(group, "_volcano_plot.pdf")),
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
      de_results[[group]]$de_results$efit,
      output_file = file.path(output_dir, group, paste0(group, "_MA_plot.pdf"))
    )

    create_pie_chart(
      de_results[[group]]$de_results$efit,
      output_file = file.path(output_dir, group, paste0(group, "_DE_pie_chart.pdf"))
    )

  }

  
  results$differential_expression <- de_results

  # ========== 4. Pathway Analysis ==========
  message("\nStep 4: Pathway analysis...")

  gsea_results <- gsea_gene_sets <- gsea_gene_ranks <- list()
  
  for (group in names(de_results)) {
      # GSEA analysis
    gsea_results[[group]] <- run_gsea_analysis(
      de_results = de_results[[group]]$de_results$efit,
      species = ifelse(species == "mouse", "MM", "HS"),
      customed_pathways = customed_pathways
    )

      # Extract gene sets and ranks once here
    gsea_gene_sets[[group]]  <- attr(gsea_results, "gene_sets")
    gsea_gene_ranks[[group]] <- attr(gsea_results, "gene_ranks")
    print(gsea_gene_sets[[group]])
    print(gsea_gene_ranks[[group]])
    
      # Create GSEA visualizations
    plot_gsea_barplot(
      gsea_results[[group]],
      output_file = file.path(output_dir, group, paste0(group, "_GSEA_barplot.pdf")),
      color_up = color_gsea_up,
      color_down = color_gsea_down,
      color_ns=color_gsea_ns
    )
    
    plot_gsea_dotplot(
      gsea_results[[group]],
      output_file = file.path(output_dir, group, paste0(group, "_GSEA_dotplot.pdf")),
      color_up = color_gsea_up,
      color_down = color_gsea_down,
      color_ns = color_gsea_ns
    )
    create_gsea_table_plot(
      gsea_results[[group]],
      output_file = file.path(output_dir, group, paste0(group, "_GSEA_tableplot.pdf")),
      gene_sets = gsea_gene_sets[[group]],
      gene_ranks = gsea_gene_ranks[[group]]
    )

    
    plot_gsea_enrichment(
      gsea_results[[group]],
      pathways    = gsea_custom_pathways,   # NULL ⇒ auto top up/down
      n_up        = n_gsea_enrich_up,
      n_down      = n_gsea_enrich_down,
      output_file = file.path(output_dir, group, paste0(group, "_GSEA_enrichment.pdf")),
      gene_sets   = gsea_gene_sets[[group]],
      gene_ranks  = gsea_gene_ranks[[group]]
    )

    # Save pathway results
    save_pathway_results(gsea_results[[group]], NULL, group, file.path(output_dir, group))
  }

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
