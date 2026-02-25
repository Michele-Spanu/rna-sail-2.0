#' Create Volcano Plot
#'
#' Creates a volcano plot from limma differential expression results.
#'
#' @param de_results Differential expression results (efit object or data frame)
#' @param fdr_threshold FDR threshold for significance (default: 0.1)
#' @param lfc_threshold Log fold change threshold (default: 1)
#' @param output_file Output PDF file path (optional)
#' @param n_labels Number of top genes to label (default: 10)
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
#' @return ggplot object
#' @export
create_volcano_plot <- function(
    de_results,
    fdr_threshold      = 0.1,
    lfc_threshold      = 1,
    output_file        = NULL,
    n_labels           = 10,
    width              = 10,
    height             = 8,
    color_up           = "#E31A1C",   # default red
    color_down         = "#1F78B4",   # default blue
    color_ns           = "grey50",    # default non-significant
    point_size         = 2.5,         # size of dots
    label_size         = 3.5,         # size of gene labels
    n_labels_up        = NULL,        # how many up genes to label
    n_labels_down      = NULL,        # how many down genes to label
    gene_label_column  = NULL,        # column to use for labels (e.g. "gene_name")
    highlight_genes    = NULL         # NEW: vector of genes to force-label
) {

  # Check required packages
  required_pkgs <- c("ggplot2", "ggrepel")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed")
    }
  }

  # Extract data
  if (is.data.frame(de_results)) {
    plot_data <- de_results
  } else {
    # Assume it's a limma efit/tfit object
    plot_data <- data.frame(
      gene_id   = rownames(de_results$coefficients),
      logFC     = de_results$coefficients[, 1],
      adj.P.Val = p.adjust(de_results$p.value[, 1], method = "fdr"),
      stringsAsFactors = FALSE
    )
  }

  # Ensure gene_id exists
  if (!"gene_id" %in% colnames(plot_data)) {
    plot_data$gene_id <- rownames(plot_data)
  }

  ## --------- GENE LABEL HANDLING ---------
  # Priority:
  # 1) If user gave gene_label_column and it exists → use that
  # 2) Else if "gene_name" exists → use that
  # 3) Else if "symbol" exists → use that
  # 4) Else try to strip things before first "_" (SYMBOL_ENSG...)
  # 5) Else fall back to gene_id

  if (!is.null(gene_label_column) && gene_label_column %in% colnames(plot_data)) {
    plot_data$gene_label <- plot_data[[gene_label_column]]
  } else if ("gene_name" %in% colnames(plot_data)) {
    plot_data$gene_label <- plot_data$gene_name
  } else if ("symbol" %in% colnames(plot_data)) {
    plot_data$gene_label <- plot_data$symbol
  } else {
    plot_data$gene_label <- ifelse(
      grepl("_", plot_data$gene_id),
      sub("_.*", "", plot_data$gene_id),
      plot_data$gene_id
    )
  }

  # Add significance categories
  plot_data$significance <- "Not Significant"
  plot_data$significance[
    plot_data$adj.P.Val <= fdr_threshold & plot_data$logFC >  lfc_threshold
  ] <- "Upregulated"
  plot_data$significance[
    plot_data$adj.P.Val <= fdr_threshold & plot_data$logFC < -lfc_threshold
  ] <- "Downregulated"

  ## --------- HOW MANY GENES TO LABEL (automatic top-N) ---------
  if (is.null(n_labels_up) || is.null(n_labels_down)) {
    # Split n_labels between up and down if specific numbers not given
    n_labels_up   <- floor(n_labels / 2)
    n_labels_down <- n_labels - n_labels_up
  }

  # Initialize label column
  plot_data$label <- ""

  # Top-N upregulated
  if (n_labels_up > 0) {
    top_up <- plot_data[plot_data$significance == "Upregulated", ]
    if (nrow(top_up) > 0) {
      top_up <- head(top_up[order(top_up$adj.P.Val), ], n_labels_up)
      idx_up <- plot_data$gene_id %in% top_up$gene_id
      plot_data$label[idx_up] <- plot_data$gene_label[idx_up]
    }
  }

  # Top-N downregulated
  if (n_labels_down > 0) {
    top_down <- plot_data[plot_data$significance == "Downregulated", ]
    if (nrow(top_down) > 0) {
      top_down <- head(top_down[order(top_down$adj.P.Val), ], n_labels_down)
      idx_down <- plot_data$gene_id %in% top_down$gene_id
      plot_data$label[idx_down] <- plot_data$gene_label[idx_down]
    }
  }

  ## --------- FORCE-LABEL USER-SPECIFIED GENES ---------
  if (!is.null(highlight_genes)) {
    highlight_genes <- unique(as.character(highlight_genes))

    # Matches by gene_label OR gene_id
    idx_hl <- plot_data$gene_label %in% highlight_genes |
      plot_data$gene_id    %in% highlight_genes

    # For those, ensure we have a label (don’t overwrite existing non-empty labels)
    to_fill <- idx_hl & (is.na(plot_data$label) | plot_data$label == "")
    plot_data$label[to_fill] <- plot_data$gene_label[to_fill]
  }

  # Create plot
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x     = logFC,
      y     = -log10(adj.P.Val),
      color = significance,
      label = label
    )
  ) +
    ggplot2::geom_point(alpha = 0.9, size = point_size) +
    ggplot2::scale_color_manual(
      values = c(
        "Not Significant" = color_ns,
        "Upregulated"     = color_up,
        "Downregulated"   = color_down
      ),
      name = "Significance"
    ) +
    ggplot2::geom_vline(
      xintercept = c(-lfc_threshold, lfc_threshold),
      linetype   = "dashed",
      color      = "grey40"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(fdr_threshold),
      linetype   = "dashed",
      color      = "grey40"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = "Volcano Plot - Differential Gene Expression",
      x     = expression(log[2] ~ "Fold Change"),
      y     = expression(-log[10] ~ "Adjusted P-value")
    ) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )

  # Add labels
  if (any(plot_data$label != "")) {
    p <- p + ggrepel::geom_text_repel(
      data               = plot_data[plot_data$label != "", ],
      max.overlaps       = 20,
      min.segment.length = 0.1,
      box.padding        = 0.5,
      point.padding      = 0.3,
      segment.color      = "grey50",
      size               = label_size,
      show.legend        = FALSE
    )
  }

  # Save if requested
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = p, width = width, height = height)
    message("Volcano plot saved to: ", output_file)
  }

  return(p)
}


#' Create PCA Plot
#'
#' Creates a PCA plot from expression data with sample annotations.
#'
#' @param expr_data Expression data matrix (genes x samples)
#' @param metadata Sample metadata
#' @param color_by Column name in metadata to color samples by
#' @param shape_by Column name in metadata to shape samples by (optional)
#' @param color_mapping Named vector of colors for each level of color_by (optional)
#' @param output_file Output file path (optional)
#' @param width Plot width (default: 8)
#' @param height Plot height (default: 6)
#' @return List containing ggplot object and PCA results
#' @export
create_pca_plot <- function(expr_data, metadata, color_by, shape_by = NULL,
                            color_mapping = NULL, output_file = NULL,
                            width = 8, height = 6) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required")

  message("Starting PCA plot...")

  # ── 1. Transpose & filter genes ──────────────────────────────
  expr_t <- t(expr_data)

  # Keep genes with finite values and non-zero variance
  valid_genes <- apply(expr_t, 2, function(x) all(is.finite(x)) && var(x) > 0)
  expr_filtered <- expr_t[, valid_genes, drop = FALSE]

  if (ncol(expr_filtered) == 0) {
    stop("No valid genes found for PCA (all zero variance or non-finite).")
  }

  # ── 2. Run PCA ──────────────────────────────────────────────
  pca_result <- prcomp(expr_filtered, scale. = TRUE)
  var_explained <- summary(pca_result)$importance[2, 1:2] * 100

  # ── 3. Build plot data ──────────────────────────────────────
  pca_data <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Sample = rownames(pca_result$x)
  )

  # Match metadata to samples
  sample_order <- match(gsub("^X", "", pca_data$Sample), metadata$SampleID)
  metadata_matched <- metadata[sample_order, , drop = FALSE]

  # Add metadata columns for aesthetics
  pca_data[[color_by]] <- metadata_matched[[color_by]]
  if (!is.null(shape_by)) {
    pca_data[[shape_by]] <- metadata_matched[[shape_by]]
  }

  # ── 4. Define color mapping ────────────────────────────────
  if (is.null(color_mapping)) {
    levels_color <- unique(na.omit(metadata_matched[[color_by]]))
    n_levels <- length(levels_color)

    if (n_levels <= 9) {
      color_mapping <- setNames(
        RColorBrewer::brewer.pal(n_levels, "Set1"),
        levels_color
      )
    } else {
      # Use rand_color for >9 categories
      if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
        stop("More than 9 categories detected; please install 'ComplexHeatmap' or provide custom color_mapping.")
      }
      color_mapping <- setNames(
        circlize::rand_color(n_levels),
        levels_color
      )
    }
  }
  p <- ggplot2::ggplot(
    pca_data,
    ggplot2::aes(x = PC1, y = PC2, color = .data[[color_by]])
  ) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::scale_color_manual(values = color_mapping) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = "Principal Component Analysis",
      x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "% variance)"),
      color = color_by
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )

  # Add shape aesthetic if provided
  if (!is.null(shape_by)) {
    p <- p +
      ggplot2::aes(shape = .data[[shape_by]]) +
      ggplot2::labs(shape = shape_by)
  }

  # Optional: add sample labels
  p <- p +
    ggplot2::geom_text(
      ggplot2::aes(label = Sample),
      vjust = -0.5, size = 3, show.legend = FALSE
    )

  # ── 6. Save plot ───────────────────────────────────────────
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = p, width = width, height = height)
    message("PCA plot saved to: ", output_file)
  }

  message("PCA completed successfully")
  return(list(plot = p, pca_result = pca_result, variance_explained = var_explained))
}

#' Create Expression Heatmap
#'
#' Creates a heatmap of gene expression data with sample annotations.
#'
#' @param expr_data Expression data matrix (genes x samples)
#' @param metadata Sample metadata
#' @param annotation_columns Vector of column names in metadata for annotation
#' @param n_genes Number of most variable genes to show (default: 500)
#' @param scale_data Whether to z-score scale the data (default: TRUE)
#' @param color_mapping Named list of color mappings for annotations (optional)
#' @param output_file Output file path (optional)
#' @param width Plot width (default: 10)
#' @param height Plot height (default: 8)
#' @return ComplexHeatmap object or base heatmap
#' @export
create_expression_heatmap <- function(expr_data, metadata, annotation_columns, n_genes = 500,
                                      scale_data = TRUE, color_mapping = NULL,
                                      output_file = NULL, min.width = 5, height = 8,
                                      title = NULL, col.cluster = TRUE, rank.order = TRUE,
                                      long.heatmap = FALSE) {
  message("starting heatmap ")

  gene_vars <- apply(expr_data, 1, var, na.rm = TRUE)
  # !!! Added last three pars to the fn, and here changed the top_genes criterion !!!
  top_genes <- if (rank.order) names(sort(gene_vars, decreasing = TRUE))[1:min(n_genes, nrow(expr_data))] else names(gene_vars)
  expr_subset <- expr_data[top_genes, , drop = FALSE]

  expr_scaled <- if (scale_data) t(scale(t(expr_subset))) else as.matrix(expr_subset)
  message("starting heatmap 1")

  colnames(expr_scaled) <- gsub("^X", "", colnames(expr_scaled))
  sample_order <- match(colnames(expr_scaled), metadata$SampleID)
  metadata_ordered <- metadata[sample_order, , drop = FALSE]
  metadata_ordered <- metadata_ordered[!is.na(metadata_ordered$SampleID), ]
  expr_scaled <- expr_scaled[, colnames(expr_scaled) %in% metadata_ordered$SampleID, drop = FALSE]
  message("starting heatmap 2")
  if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
      requireNamespace("circlize", quietly = TRUE)) {

    # !!! The first comma was missing here => error: ignoring drop !!!
    annotation_data <- metadata_ordered[,annotation_columns, drop = FALSE]
    col_list <- list()
    message("starting heatmap 3")
    print(annotation_columns)
    for (col in annotation_columns) {
      unique_vals <- unique(annotation_data[[col]])

      if (!is.null(color_mapping) && !is.null(color_mapping[[col]])) {
        col_list[[col]] <- color_mapping[[col]]
      } else if (is.factor(unique_vals) || is.character(unique_vals)) {
        # !!! Palette creation within Set1 => further colours can be interpolated !!!
        pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
        # !!! Must check grDevices is installed !!!
        colors <- pal(length(unique_vals))
        col_list[[col]] <- setNames(colors, unique_vals)
      }
    }
    message("starting heatmap 4")
    print(col_list)
    ha <- ComplexHeatmap::HeatmapAnnotation(
      df = annotation_data,
      col = col_list
    )
    message("starting heatmap 5")

    if (!is.null(output_file)) {
        # !!! Creates the folder if it doesn't exist (useful if fn called by itself) !!!
        dir.create(
        dirname(output_file),
        recursive = TRUE,
        showWarnings = FALSE
        )
        # !!! Dynamic width and length!!!
        height <- if (long.heatmap) 4+0.15*lenght(top_genes) else height
        pdf(output_file, width = min.width+0.18*ncol(expr_scaled), height = height)
    }
    message("starting heatmap 6")

    # !!! Fixed colour intensity scale for non-norm. data, show_row_names, row_names_side, row_name_gp, col.cluster, title !!!
    ht <- ComplexHeatmap::Heatmap(
      expr_scaled,
      name = ifelse(scale_data, "Expression\n(Z-score)", "Expression"),
      col = circlize::colorRamp2(
        if (scale_data) c(-2, 0, 2) else {c(min(expr_scaled, na.rm = TRUE), 0, max(expr_scaled, na.rm = TRUE))},
        c("blue", "white", "red")
      ),
      top_annotation = ha,
      show_row_names = long.heatmap,
      row_names_side = "left",
      row_name_gp = grid::gpar(fontsize = 10)
      show_column_names = TRUE,
      cluster_rows = TRUE,
      cluster_columns = col.cluster,
      column_title = if (!is.null(title)) title else paste("Expression Heatmap -", n_genes, "Most Variable Genes"),
      heatmap_legend_param = list(direction = "vertical")
    )

    print(ht)
    if (!is.null(output_file)) {
      dev.off()
      message("Expression heatmap saved to: ", output_file)
    }

    return(ht)

  } else {
    stop("ComplexHeatmap and circlize packages are required for heatmap plotting.")
  }
}
#' Create MA Plot
#'
#' Creates an MA plot (log ratio vs mean average) from differential expression results.
#'
#' @param de_results Differential expression results
#' @param fdr_threshold FDR threshold for significance (default: 0.05)
#' @param output_file Output file path (optional)
#' @param width Plot width (default: 8)
#' @param height Plot height (default: 6)
#' @return ggplot object
#' @export
create_ma_plot <- function(de_results, fdr_threshold = 0.05, output_file = NULL,
                          width = 8, height = 6) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required")
  }

  # Extract data
  if (is.data.frame(de_results)) {
    plot_data <- de_results
    if (!"AveExpr" %in% colnames(plot_data)) {
      stop("AveExpr column not found in results data frame")
    }
  } else {
    # limma efit object
    results_table <- limma::topTable(de_results, number = Inf)
    plot_data <- results_table
  }

  # Add significance
  plot_data$significant <- plot_data$adj.P.Val < fdr_threshold

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = AveExpr, y = logFC, color = significant)) +
    ggplot2::geom_point(alpha = 0.6, size = 1) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "grey50", "TRUE" = "#E31A1C"),
      name = paste0("FDR < ", fdr_threshold),
      labels = c("Not Significant", "Significant")
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = "MA Plot - Mean vs Log Fold Change",
      x = "Average Expression",
      y = expression(log[2]~"Fold Change")
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )

  # Save if requested
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = p, width = width, height = height)
    message("MA plot saved to: ", output_file)
  }

  return(p)
}

#' Create Pie Chart
#'
#' Creates a pie chart showing distribution of up/down regulated genes.
#'
#' @param de_results Differential expression results
#' @param fdr_threshold FDR threshold (default: 0.05)
#' @param lfc_threshold Log fold change threshold (default: 1)
#' @param output_file Output file path (optional)
#' @param width Plot width (default: 6)
#' @param height Plot height (default: 6)
#' @return ggplot object
#' @export
create_pie_chart <- function(de_results, fdr_threshold = 0.05, lfc_threshold = 1,
                            output_file = NULL, width = 6, height = 6) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required")
  }

  # Get significance counts
  if (is.data.frame(de_results)) {
    results_data <- de_results
  } else {
    # Use limma decideTests
    decisions <- limma::decideTests(de_results, p.value = fdr_threshold, lfc = lfc_threshold)
    up_count <- sum(decisions == 1)
    down_count <- sum(decisions == -1)
    ns_count <- sum(decisions == 0)

    results_data <- data.frame(
      category = c("Upregulated", "Downregulated", "Not Significant"),
      count = c(up_count, down_count, ns_count)
    )
  }

  # Calculate percentages
  results_data$percentage <- round(results_data$count / sum(results_data$count) * 100, 1)
  results_data$label <- paste0(results_data$category, "\n", results_data$count, " (",
                              results_data$percentage, "%)")

  # Create pie chart
  p <- ggplot2::ggplot(results_data, ggplot2::aes(x = "", y = count, fill = category)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::scale_fill_manual(
      values = c("Upregulated" = "#E31A1C", "Downregulated" = "#1F78B4",
                "Not Significant" = "grey70"),
      name = "Gene Category"
    ) +
    ggplot2::theme_void() +
    ggplot2::labs(title = "Distribution of Differentially Expressed Genes") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "bottom"
    ) +
    ggplot2::geom_text(ggplot2::aes(label = count),
                      position = ggplot2::position_stack(vjust = 0.5),
                      color = "white", size = 4, fontface = "bold")

  # Save if requested
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = p, width = width, height = height)
    message("Pie chart saved to: ", output_file)
  }

  return(p)
}

#' Create Multi-Panel Summary Plot
#'
#' Creates a comprehensive summary plot combining multiple visualizations.
#'
#' @param expr_data Expression data matrix
#' @param de_results Differential expression results
#' @param metadata Sample metadata
#' @param condition_column Condition column name in metadata
#' @param output_file Output file path
#' @param width Plot width (default: 16)
#' @param height Plot height (default: 12)
#' @return Combined plot object
#' @export
create_summary_plot <- function(expr_data, de_results, metadata, condition_column,
                               output_file, width = 16, height = 12) {

  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package 'gridExtra' is required for summary plots")
  }

  # Create individual plots
  message("Creating PCA plot...")
  pca_plot <- create_pca_plot(expr_data, metadata, color_by = condition_column)$plot +
    ggplot2::theme(legend.position = "none")

  message("Creating volcano plot...")
  volcano_plot <- create_volcano_plot(de_results, n_labels = 5,fdr_threshold = 0.1, lfc_threshold = 1) +
    ggplot2::theme(legend.position = "none")

  message("Creating MA plot...")
  ma_plot <- create_ma_plot(de_results) +
    ggplot2::theme(legend.position = "none")

  message("Creating pie chart...")
  pie_plot <- create_pie_chart(de_results) +
    ggplot2::theme(legend.position = "none")

  # Combine plots
  pdf(output_file, width = width, height = height)

  combined_plot <- gridExtra::grid.arrange(
    pca_plot, volcano_plot,
    ma_plot, pie_plot,
    ncol = 2, nrow = 2,
    top = grid::textGrob("RNA-seq Analysis Summary",
                        gp = grid::gpar(fontsize = 16, fontface = "bold"))
  )

  print(combined_plot)
  dev.off()

  message("Summary plot saved to: ", output_file)

  return(combined_plot)
}#' Create Volcano Plot
#'
#' Creates a volcano plot from limma differential expression results.
#
