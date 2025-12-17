#' Run DoRothEA Transcription Factor Activity Analysis
#'
#' Performs transcription factor activity analysis using DoRothEA regulons and VIPER.
#'
#' @param expr_data Log-transformed expression data (genes x samples)
#' @param metadata Sample metadata
#' @param condition_column Column name for condition comparison
#' @param group1 Name of group 1 (treatment)
#' @param group2 Name of group 2 (control)
#' @param species Species ("mouse" or "human", default: "mouse")
#' @param confidence_levels DoRothEA confidence levels to include (default: c("A", "B", "C"))
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return List containing TF activity scores and differential TF analysis results
#' @export

run_tf_analysis <- function(expr_data, metadata, condition_column, group1, group2,
                           species = "mouse", confidence_levels = c("A", "B", "C"),
                           output_dir, experiment_name,fdr_threshold = 0.05, lfc_threshold = 1) {

  # Check required packages
  required_pkgs <- c("dorothea", "viper", "biomaRt", "limma")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed")
    }
  }

  create_output_dir(output_dir)
  message("Starting transcription factor activity analysis...")

  # Load DoRothEA regulons
  regulon <- load_dorothea_regulons(species, confidence_levels)

  # Convert gene symbols to ENSEMBL IDs
  regulon_ensembl <- map_regulons_to_ensembl(regulon, species)

  # Prepare expression data with gene symbols
  expr_symbols <- expr_data

  rownames(expr_symbols) <- sub(".*_", "", rownames(expr_symbols))  # Extract gene symbol
  expr_symbols <- as.matrix(expr_symbols)   # force numeric
  mode(expr_symbols) <- "numeric"

  # Optional: remove genes with zero variance
  expr_symbols <- expr_symbols[apply(expr_symbols, 1, sd) != 0, ]
  # Z-score normalize for VIPER
  expr_z <- t(scale(t(expr_symbols)))

  # Convert regulon to VIPER format
  regulon_viper <- convert_to_viper_format(regulon_ensembl)
  # regulon_viper <- lapply(regulon_viper, function(x) {
  #   x$tfmode <- as.numeric(x$tfmode)
  #   x$likelihood <- as.numeric(x$likelihood)
  #   return(x)
  # })

  # Run VIPER
  message("Running VIPER to estimate TF activities...")
  #print(expr_z)
  #print(regulon_viper)
  tf_activities <- viper::viper(expr_z, regulon_viper, verbose = TRUE)

  # Differential TF activity analysis
  message("Performing differential TF activity analysis...")
  diff_tf_results <- analyze_differential_tf_activity(tf_activities, metadata, condition_column,
                                                     group1, group2,fdr_threshold = 0.05, lfc_threshold = 1)

  # Create visualizations
  plot_tf_volcano <- function(diff_tf_results, output_dir, experiment_name,
                          fdr_threshold = 0.05, lfc_threshold = 1, n_labels = 20) {

  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("ggrepel", quietly = TRUE)) {
    message("Skipping TF volcano plot - required packages not available")
    return()
  }

  # Prepare data
  plot_data <- diff_tf_results
  plot_data$neg_log10_padj <- -log10(plot_data$adj.P.Val)

  # Select top TFs to label
  plot_data$label <- ""
  top_activated <- head(plot_data[plot_data$Regulation == "Activated", ], n_labels/2)
  top_repressed <- head(plot_data[plot_data$Regulation == "Repressed", ], n_labels/2)

  plot_data$label[plot_data$TF %in% c(top_activated$TF, top_repressed$TF)] <-
    plot_data$TF[plot_data$TF %in% c(top_activated$TF, top_repressed$TF)]

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = logFC, y = neg_log10_padj, color = Regulation)) +
    ggplot2::geom_point(alpha = 0.7, size = 2) +
    ggplot2::scale_color_manual(values = c(
      "Activated" = "#E31A1C",
      "Repressed" = "#1F78B4",
      "Not Significant" = "grey70"
    )) +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                       linetype = "dashed", color = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(fdr_threshold),
                       linetype = "dashed", color = "grey40") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = "Transcription Factor Activity Changes",
      x = "Log2 Fold Change (TF Activity)",
      y = "-log10(Adjusted P-value)",
      color = "Regulation"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )

  # Add labels if ggrepel is available
  if (any(plot_data$label != "")) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = label),
      max.overlaps = 20,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey50",
      size = 3
    )
  }

  # Save plot
  output_file <- file.path(output_dir, paste0(experiment_name, "_TF_volcano_plot.pdf"))
  ggplot2::ggsave(output_file, p, width = 10, height = 8)

  message("TF volcano plot saved to: ", output_file)
  return(p)
}

#' Plot TF Activity Barplot
#'
#' Creates a barplot of top differentially active transcription factors.
#'
#' @param diff_tf_results Differential TF analysis results
#' @param n_tfs Number of top TFs to show (default: 20)
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return ggplot object
#' @keywords internal
plot_tf_barplot <- function(diff_tf_results, n_tfs = 20, output_dir, experiment_name) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Skipping TF barplot - ggplot2 not available")
    return()
  }

  # Select top TFs by absolute log fold change and significance
  sig_tfs <- diff_tf_results[diff_tf_results$adj.P.Val < 0.05, ]

  if (nrow(sig_tfs) == 0) {
    message("No significant TFs found for barplot")
    return()
  }

  # Get top activated and repressed
  n_each <- min(n_tfs/2, nrow(sig_tfs))
  top_activated <- head(sig_tfs[sig_tfs$logFC > 0, ], n_each)
  top_repressed <- head(sig_tfs[sig_tfs$logFC < 0, ], n_each)

  plot_data <- rbind(top_activated, top_repressed)
  plot_data$TF <- factor(plot_data$TF, levels = plot_data$TF[order(plot_data$logFC)])

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = TF, y = logFC, fill = Regulation)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = c(
      "Activated" = "#E31A1C",
      "Repressed" = "#1F78B4"
    )) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(
      title = "Top Differentially Active Transcription Factors",
      x = "Transcription Factor",
      y = "Log2 Fold Change (TF Activity)",
      fill = "Regulation"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")

  # Save plot
  output_file <- file.path(output_dir, paste0(experiment_name, "_TF_barplot.pdf"))
  ggplot2::ggsave(output_file, p, width = 8, height = 6)

  message("TF barplot saved to: ", output_file)
  return(p)
}

#' Run Master Regulator Analysis
#'
#' Identifies master transcription factors using VIPER's msviper function.
#'
#' @param diff_tf_results Differential TF analysis results
#' @param regulon_viper VIPER-formatted regulon
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return Master regulator analysis results
#' @keywords internal
run_master_regulator_analysis <- function(diff_tf_results, regulon_viper, output_dir, experiment_name) {

  if (!requireNamespace("viper", quietly = TRUE)) {
    message("Skipping master regulator analysis - viper not available")
    return(NULL)
  }

  message("Running master regulator analysis...")

  # Create signature from t-statistics
  tf_stats <- diff_tf_results$t
  names(tf_stats) <- diff_tf_results$TF

  # Run msviper
  tryCatch({
    mra_results <- viper::msviper(
      signature = tf_stats,
      regulon = regulon_viper,
      minsize = 10,
      ges.filter = FALSE
    )

    # Create visualization
    output_file <- file.path(output_dir, paste0(experiment_name, "_master_regulators.pdf"))
    pdf(output_file, width = 15, height = 8)
    plot(mra_results, main = paste("Master Regulators -", experiment_name))
    dev.off()

    message("Master regulator analysis plot saved to: ", output_file)

    # Extract summary
    mra_summary <- summary(mra_results)

    return(list(results = mra_results, summary = mra_summary))

  }, error = function(e) {
    message("Master regulator analysis failed: ", e$message)
    return(NULL)
  })
}

#' Save TF Analysis Results
#'
#' Saves all transcription factor analysis results to files.
#'
#' @param tf_activities TF activity matrix
#' @param diff_tf_results Differential TF analysis results
#' @param master_regulators Master regulator analysis results
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return None (invisible)
#' @keywords internal
save_tf_results <- function(tf_activities, diff_tf_results, master_regulators, output_dir, experiment_name) {

  # Save TF activities
  tf_file <- file.path(output_dir, paste0(experiment_name, "_TF_activities.tsv"))
  write.table(tf_activities, file = tf_file, sep = "\t", quote = FALSE, row.names = TRUE)

  # Save differential results
  diff_file <- file.path(output_dir, paste0(experiment_name, "_differential_TF_results.tsv"))
  write.table(diff_tf_results, file = diff_file, sep = "\t", quote = FALSE, row.names = FALSE)

  # Save significant TFs
  sig_activated <- diff_tf_results[diff_tf_results$Regulation == "Activated", "TF"]
  sig_repressed <- diff_tf_results[diff_tf_results$Regulation == "Repressed", "TF"]

  if (length(sig_activated) > 0) {
    activated_file <- file.path(output_dir, paste0(experiment_name, "_activated_TFs.txt"))
    writeLines(sig_activated, activated_file)
  }

  if (length(sig_repressed) > 0) {
    repressed_file <- file.path(output_dir, paste0(experiment_name, "_repressed_TFs.txt"))
    writeLines(sig_repressed, repressed_file)
  }

  # Save master regulator results if available
  if (!is.null(master_regulators)) {
    mra_file <- file.path(output_dir, paste0(experiment_name, "_master_regulators.rds"))
    saveRDS(master_regulators, mra_file)
  }

  message("TF analysis results saved to: ", output_dir)

  invisible()
}

#' Create TF Network Visualization
#'
#' Creates a network visualization of top transcription factors and their targets.
#'
#' @param diff_tf_results Differential TF analysis results
#' @param regulon Original regulon data frame
#' @param n_tfs Number of top TFs to include (default: 10)
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return None (creates plot)
#' @export

create_tf_network_plot <- function(diff_tf_results, regulon, n_tfs = 10, output_dir, experiment_name) {

  if (!requireNamespace("igraph", quietly = TRUE) ||
      !requireNamespace("ggraph", quietly = TRUE)) {
    message("Skipping TF network plot - igraph or ggraph not available")
    return()
  }

  # Get top significant TFs
  sig_tfs <- diff_tf_results[diff_tf_results$adj.P.Val < 0.05, ]
  top_tfs <- head(sig_tfs[order(-abs(sig_tfs$logFC)), ], n_tfs)

  if (nrow(top_tfs) == 0) {
    message("No significant TFs found for network plot")
    return()
  }

  # Filter regulon for top TFs
  network_regulon <- regulon[regulon$tf %in% top_tfs$TF, ]

  # Limit targets per TF for visualization
  max_targets_per_tf <- 10
  network_regulon <- network_regulon %>%
    dplyr::group_by(tf) %>%
    dplyr::slice_head(n = max_targets_per_tf) %>%
    dplyr::ungroup()

  # Create igraph object
  edges <- network_regulon[, c("tf", "target", "mor")]
  graph <- igraph::graph_from_data_frame(edges, directed = TRUE)

  # Add node attributes
  igraph::V(graph)$type <- ifelse(igraph::V(graph)$name %in% top_tfs$TF, "TF", "Target")
  igraph::V(graph)$logFC <- NA
  igraph::V(graph)$logFC[match(top_tfs$TF, igraph::V(graph)$name)] <- top_tfs$logFC

  # Create plot
  output_file <- file.path(output_dir, paste0(experiment_name, "_TF_network.pdf"))

  pdf(output_file, width = 12, height = 10)

  p <- ggraph::ggraph(graph, layout = "stress") +
    ggraph::geom_edge_link(ggplot2::aes(color = factor(mor)),
                          arrow = ggplot2::arrow(length = ggplot2::unit(2, "mm")),
                          alpha = 0.6) +
    ggraph::geom_node_point(ggplot2::aes(color = type, size = type)) +
    ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 3) +
    ggplot2::scale_color_manual(
      name = "Node Type",
      values = c("TF" = "#E31A1C", "Target" = "#377EB8")
    ) +
    ggplot2::scale_size_manual(
      name = "Node Type",
      values = c("TF" = 4, "Target" = 2)
    ) +
    ggraph::scale_edge_color_manual(
      name = "Regulation",
      values = c("1" = "#E31A1C", "-1" = "#1F78B4"),
      labels = c("Activation", "Repression")
    ) +
    ggplot2::labs(title = paste("Transcription Factor Network -", experiment_name)) +
    ggraph::theme_graph() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )

  print(p)
  dev.off()

  message("TF network plot saved to: ", output_file)
}
  print(head(tf_activities))

  plot_tf_activity_heatmap(tf_activities, metadata, condition_column, output_dir, experiment_name)

  plot_tf_volcano(diff_tf_results, output_dir, experiment_name,fdr_threshold = 0.05, lfc_threshold = 1)

  plot_tf_barplot(diff_tf_results, n_tfs = 20, output_dir, experiment_name)


  # Master regulator analysis
  if (requireNamespace("viper", quietly = TRUE)) {
    master_regulators <- run_master_regulator_analysis(diff_tf_results, regulon_viper, output_dir, experiment_name)
  } else {
    master_regulators <- NULL
  }

  # Save results
  save_tf_results(tf_activities, diff_tf_results, master_regulators, output_dir, experiment_name)

  message("Transcription factor analysis completed!")

  return(list(
    tf_activities = tf_activities,
    differential_results = diff_tf_results,
    master_regulators = master_regulators,
    regulon = regulon_viper
  ))
}

#' Load DoRothEA Regulons
#'
#' Loads DoRothEA transcription factor regulons.
#'
#' @param species Species ("mouse" or "human")
#' @param confidence_levels Confidence levels to include
#' @return DoRothEA regulon data frame
#' @keywords internal
load_dorothea_regulons <- function(species, confidence_levels) {

  if (!requireNamespace("dorothea", quietly = TRUE)) {
    stop("Package 'dorothea' is required")
  }

  if (species == "mouse") {
    data("dorothea_mm", package = "dorothea")
    regulon <- dorothea_mm
  } else if (species == "human") {
    data("dorothea_hs", package = "dorothea")
    regulon <- dorothea_hs
  } else {
    stop("Species must be 'mouse' or 'human'")
  }

  # Filter by confidence levels
  regulon <- regulon[regulon$confidence %in% confidence_levels, ]

  message("Loaded ", nrow(regulon), " TF-target interactions from DoRothEA")
  message("Confidence levels: ", paste(confidence_levels, collapse = ", "))

  return(regulon)
}

#' Map Regulons to ENSEMBL IDs
#'
#' Maps gene symbols in regulons to ENSEMBL IDs using biomaRt.
#'
#' @param regulon DoRothEA regulon data frame
#' @param species Species ("mouse" or "human")
#' @return Regulon with ENSEMBL IDs
#' @keywords internal
map_regulons_to_ensembl <- function(regulon, species) {

  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    warning("biomaRt not available, using gene symbols directly")
    return(regulon)
  }

  # Set up biomaRt
  if (species == "mouse") {
    mart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    symbol_attr <- "mgi_symbol"
  } else {
    mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    symbol_attr <- "hgnc_symbol"
  }

  # Map TF symbols to ENSEMBL
  tf_symbols <- unique(regulon$tf)
  tf_mapping <- biomaRt::getBM(
    attributes = c(symbol_attr, "ensembl_gene_id"),
    filters = symbol_attr,
    values = tf_symbols,
    mart = mart
  )
  tf_mapping <- tf_mapping[tf_mapping[[symbol_attr]] != "" & tf_mapping$ensembl_gene_id != "", ]
  tf_mapping_vec <- setNames(tf_mapping$ensembl_gene_id, tf_mapping[[symbol_attr]])

  # Map target symbols to ENSEMBL
  target_symbols <- unique(regulon$target)
  target_mapping <- biomaRt::getBM(
    attributes = c(symbol_attr, "ensembl_gene_id"),
    filters = symbol_attr,
    values = target_symbols,
    mart = mart
  )
  target_mapping <- target_mapping[target_mapping[[symbol_attr]] != "" & target_mapping$ensembl_gene_id != "", ]
  target_mapping_vec <- setNames(target_mapping$ensembl_gene_id, target_mapping[[symbol_attr]])

  # Map regulon
  regulon$tf_ensembl <- tf_mapping_vec[regulon$tf]
  regulon$target_ensembl <- target_mapping_vec[regulon$target]

  # Remove unmapped entries
  regulon <- regulon[!is.na(regulon$tf_ensembl) & !is.na(regulon$target_ensembl), ]

  message("Mapped to ENSEMBL IDs: ", nrow(regulon), " interactions remain")

  return(regulon)
}

#' Convert to VIPER Format
#'
#' Converts DoRothEA regulons to VIPER format.
#'
#' @param regulon Regulon data frame with ENSEMBL IDs
#' @return VIPER-formatted regulon list
#' @keywords internal
convert_to_viper_format <- function(regulon) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required")
  }

  # Group by TF and create VIPER format using Ensembl IDs
  regulon_viper <- regulon %>%
    dplyr::group_by(tf) %>%
    dplyr::summarise(
      tfmode = list(setNames(mor, target_ensembl)),        # <-- use target_ensembl
      likelihood = list(setNames(rep(1, length(target_ensembl)), target_ensembl)), # <-- use target_ensembl
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      tf_regulon = purrr::map2(tfmode, likelihood, ~ list(tfmode = .x, likelihood = .y))
    ) %>%
    { setNames(.$tf_regulon, .$tf) }

  message("Created VIPER regulon with ", length(regulon_viper), " transcription factors")

  return(regulon_viper)
}

#' Analyze Differential TF Activity
#'
#' Performs differential analysis of transcription factor activities.
#'
#' @param tf_activities TF activity matrix from VIPER
#' @param metadata Sample metadata
#' @param condition_column Condition column name
#' @param group1 Group 1 name
#' @param group2 Group 2 name
#' @return Differential TF analysis results
#' @keywords internal
analyze_differential_tf_activity <- function(tf_activities, metadata, condition_column, group1, group2,fdr_threshold = 0.05, lfc_threshold = 1) {

  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required")
  }

  # Match metadata to TF activities
  sample_order <- match(colnames(tf_activities), metadata$SampleID)
  metadata_ordered <- metadata[sample_order, ]

  # Create condition factor
  condition <- factor(metadata_ordered[[condition_column]], levels = c(group2, group1))


  # Keep only samples with non-NA condition
  keep <- !is.na(condition)

  tf_activities_sub <- tf_activities[, keep]
  condition_sub <- droplevels(condition[keep])

  # Redefine design matrix
  design <- model.matrix(~ condition_sub)
  colnames(design) <- c("Intercept", paste0(group1, "_vs_", group2))

  # Fit model
  fit <- limma::lmFit(tf_activities_sub, design)
  fit <- limma::eBayes(fit)

  # Extract results
  results <- limma::topTable(fit, coef = 2, number = Inf)
  results$TF <- rownames(results)

  # Add significance labels
  results$Regulation <- "Not Significant"
  results$Regulation[results$adj.P.Val < fdr_threshold & results$logFC > lfc_threshold] <- "Activated"
  results$Regulation[results$adj.P.Val < fdr_threshold & results$logFC < lfc_threshold] <- "Repressed"

  message("Found ", sum(results$Regulation == "Activated"), " activated and ",
          sum(results$Regulation == "Repressed"), " repressed TFs (adj.P.Val < ",fdr_threshold,")")

  return(results)
}

#' Plot TF Activity Heatmap
#'
#' Creates a heatmap of transcription factor activities.
#'
#' @param tf_activities TF activity matrix
#' @param metadata Sample metadata
#' @param condition_column Condition column name
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @param n_top_tfs Number of most variable TFs to show (default: 50)
#' @return None (creates plot)
#' @keywords internal
plot_tf_activity_heatmap <- function(tf_activities,
                                     metadata,
                                     condition_column,
                                     output_dir,
                                     experiment_name,
                                     n_top_tfs = 50) {

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    message("Skipping TF heatmap - ComplexHeatmap not available")
    return()
  }

  # Select most variable TFs
  tf_vars <- apply(tf_activities, 1, var, na.rm = TRUE)

  top_tfs <- names(sort(tf_vars, decreasing = TRUE))[1:min(n_top_tfs, length(tf_vars))]

  tf_subset <- tf_activities[top_tfs, , drop = FALSE]


  # Z-score scaling
  tf_scaled <- t(scale(t(tf_subset)))


  # Match metadata to columns
  sample_order <- match(colnames(tf_scaled), metadata$SampleID)
  metadata_ordered <- metadata[sample_order, , drop = FALSE]

  # Extract condition vector
  cond_vec <- metadata_ordered[[condition_column]]

  # Unique non-NA conditions
  conds <- sort(unique(na.omit(cond_vec)))
  n_conds <- length(conds)

  if (n_conds == 0) {
    stop("No non-NA conditions found in metadata for column: ", condition_column)
  }

  if (n_conds == 1) {
    # Single condition – just pick one color
    cond_cols <- "#E41A1C"
  } else if (n_conds <= 9) {
    # Brewer needs n >= 3, so take first n_conds colors
    cond_cols <- RColorBrewer::brewer.pal(max(n_conds, 3), "Set1")[seq_len(n_conds)]
  } else {
    cond_cols <- circlize::rand_color(n_conds)
  }

  col_colors <- setNames(cond_cols, conds)
  print(col_colors)

  # Column annotation (NA conditions will just have no color)
  ha <- ComplexHeatmap::HeatmapAnnotation(
    Condition = factor(cond_vec, levels = conds),
    col = list(Condition = col_colors)
  )


  # Create heatmap
  output_file <- file.path(output_dir, paste0(experiment_name, "_TF_activity_heatmap.pdf"))

  pdf(output_file, width = 10, height = 12)
  ht <- ComplexHeatmap::Heatmap(
    tf_scaled,
    name = "TF Activity\n(Z-score)",
    col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    top_annotation = ha,
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_names_gp = grid::gpar(fontsize = 8),
    column_title = paste("Top", n_top_tfs, "Variable Transcription Factors")
  )
  print(ht)
  dev.off()

  message("TF activity heatmap saved to: ", output_file)
}


