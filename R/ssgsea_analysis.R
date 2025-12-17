#' Run ssGSEA analysis
#'
#' Performs single-sample GSEA (ssGSEA) using GSVA on Hallmark gene sets
#' from MSigDB by default, with the option to add custom pathways.
#'
#' @param expr_data Expression matrix (genes x samples), typically log-TPM
#' @param metadata Sample metadata data.frame
#' @param condition_column Column in metadata with group labels (e.g. "condition")
#' @param sample_id_column Column in metadata with sample IDs matching colnames(expr_data)
#' @param group1 Name of group 1 (e.g. "mut")
#' @param group2 Name of group 2 (e.g. "wt")
#' @param species "mouse" or "human"
#' @param extra_pathways Optional named list of extra pathways
#'        (names = pathway names, values = vectors of Ensembl IDs).
#'        These are added on top of the Hallmark sets.
#' @param min_size Minimum gene set size (default: 15)
#' @param max_size Maximum gene set size (default: 500)
#' @param output_dir Output directory
#' @param experiment_name Experiment name (used for file names)
#' @param n_boxplot_pathways Number of pathways to include in boxplot (top by FDR, default: 20)
#'
#' @return List with:
#'   \item{scores}{matrix of ssGSEA scores (pathways x samples)}
#'   \item{stats}{data.frame with per-pathway Wilcoxon test results}
#'   \item{pathways}{list of pathways used}
#' @export
run_ssgsea_analysis <- function(
    expr_data,
    metadata,
    condition_column,
    sample_id_column,
    group1,
    group2,
    species          = "mouse",
    extra_pathways   = NULL,
    min_size         = 15,
    max_size         = 500,
    output_dir,
    experiment_name,
    n_boxplot_pathways = 20
) {
  
  required_pkgs <- c("GSVA", "msigdbr", "ggplot2", "pheatmap", "ggpubr")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed")
    }
  }
  
  message("Running ssGSEA analysis...")
  
  # Ensure sample order consistency
  if (!all(colnames(expr_data) %in% metadata[[sample_id_column]])) {
    stop("Not all expression sample names are present in metadata.")
  }
  
  # Reorder metadata to match expression columns
  metadata <- metadata[match(colnames(expr_data), metadata[[sample_id_column]]), ]
  
  # ----- 1) Build pathway list -----
  if (species == "mouse") {
    species1  <- "Mus musculus"
    category  <- "MH"
  } else {
    species1  <- "human"
    category  <- "H"
  }
  
  message("Loading Hallmark gene sets for ", species1, " (", category, ")")
  msig_df <- msigdbr::msigdbr(
    db_species = ifelse(species == "mouse", "MM", "HS"),
    species    = species1,
    collection = category
  )
  
  hallmark_sets <- split(msig_df$ensembl_gene, msig_df$gs_name)
  # Add user-provided pathways if given
  if (!is.null(extra_pathways)) {
    if (!is.list(extra_pathways) || is.null(names(extra_pathways))) {
      stop("extra_pathways must be a *named* list: list(PATHWAY1 = c('ENSG...'), PATHWAY2 = ...)")
    }
    hallmark_sets <- c(hallmark_sets, extra_pathways)
    message("Added ", length(extra_pathways), " user-defined pathways to Hallmark sets")
  }
  
  # Filter to genes present in expression matrix
  # --- Normalize gene IDs: "TP53_ENSG00000..." -> "ENSG00000..."
  rn <- rownames(expr_data)
  
  # take the last "_" chunk (works even if symbol contains underscores)
  ensg <- sub(".*_", "", rn)
  
  # sanity check
  if (!all(grepl("^ENSG", ensg))) {
    warning("Not all rownames look like SYMBOL_ENSG...; check rownames(expr_data).")
  }
  
  rownames(expr_data) <- ensg
  gene_ids <- ensg
  
  
  
  #gene_ids <- rownames(expr_data)
  pathways_filtered <- lapply(hallmark_sets, function(gs) {
    intersect(gs, gene_ids)
  })
  pathways_filtered <- pathways_filtered[sapply(pathways_filtered, length) >= min_size &
                                           sapply(pathways_filtered, length) <= max_size]
  if (length(pathways_filtered) == 0) {
    stop(
      "0 pathways after filtering. Likely gene ID mismatch.\n",
      "Example expr rownames: ", paste(head(rownames(expr_data)), collapse=", "), "\n",
      "Expected either ENSG... (for msigdbr$ensembl_gene) or symbols (for msigdbr$gene_symbol)."
    )
  }
  message("Using ", length(pathways_filtered),
          " pathways after filtering by size and gene overlap")
  
  # ----- 2) Run ssGSEA via GSVA -----
  message("Running GSVA::gsva (method = 'ssgsea')...")
  
  ssgsea_param <- GSVA::ssgseaParam(
    exprData  = as.matrix(expr_data),
    geneSets  = pathways_filtered,
    minSize   = min_size,
    maxSize   = max_size,
    alpha     = 0.25,   # classic ssGSEA tau; you can expose this as an argument later
    normalize = TRUE
  )
  
  ssgsea_scores <- GSVA::gsva(ssgsea_param, verbose = TRUE)
  ssgsea_scores <-as.matrix(ssgsea_scores)
  # ssgsea_scores: pathways x samples
  
  # ----- 3) Save matrix -----
  scores_file <- file.path(output_dir, paste0(experiment_name, "_ssGSEA_scores.tsv"))
  write.table(ssgsea_scores, file = scores_file, sep = "\t", quote = FALSE)
  message("ssGSEA scores saved to: ", scores_file)
  
  # ----- 4) Heatmap of ssGSEA scores -----
  message("Creating ssGSEA heatmap...")
  
  # Row-scale scores for visualization
  scores_scaled <- t(scale(t(ssgsea_scores)))
  scores_scaled[is.na(scores_scaled)] <- 0
  
  annotation_col <- data.frame(
    Group = metadata[[condition_column]]
  )
  rownames(annotation_col) <- metadata[[sample_id_column]]
  
  heatmap_file <- file.path(output_dir, paste0(experiment_name, "_ssGSEA_heatmap.pdf"))
  grDevices::pdf(heatmap_file, width = 10, height = 8)
  pheatmap::pheatmap(
    scores_scaled,
    annotation_col = annotation_col,
    cluster_rows   = TRUE,
    cluster_cols   = TRUE,
    show_rownames  = TRUE,
    show_colnames  = FALSE,
    fontsize_row   = 6,
    main           = "ssGSEA – pathway activity"
  )
  print(heatmap_file)
  grDevices::dev.off()
  message("ssGSEA heatmap saved to: ", heatmap_file)
  
  # ----- 5) Per-pathway statistics (group1 vs group2) -----
  message("Computing pathway score statistics (", group1, " vs ", group2, ")...")
  
  scores_long <- data.frame(
    sample  = rep(colnames(ssgsea_scores), each = nrow(ssgsea_scores)),
    pathway = rep(rownames(ssgsea_scores), times = ncol(ssgsea_scores)),
    score   = as.vector(ssgsea_scores),
    stringsAsFactors = FALSE
  )
  
  meta_sub <- metadata[, c(sample_id_column, condition_column)]
  colnames(meta_sub)[colnames(meta_sub) == sample_id_column] <- "sample"
  
  scores_long <- merge(scores_long, meta_sub, by = "sample", all.x = TRUE)
  colnames(scores_long)[colnames(scores_long) == condition_column] <- "Group"
  
  pathways <- unique(scores_long$pathway)
  stats_list <- vector("list", length(pathways))
  
  for (i in seq_along(pathways)) {
    pw <- pathways[i]
    df_pw <- scores_long[scores_long$pathway == pw, ]
    g1 <- df_pw$score[df_pw$Group == group1]
    g2 <- df_pw$score[df_pw$Group == group2]
    
    pval <- NA_real_
    if (length(g1) > 1 && length(g2) > 1) {
      pval <- tryCatch(
        stats::wilcox.test(g1, g2)$p.value,
        error = function(e) NA_real_
      )
    }
    
    stats_list[[i]] <- data.frame(
      pathway      = pw,
      median_group1 = median(g1, na.rm = TRUE),
      median_group2 = median(g2, na.rm = TRUE),
      p_value       = pval,
      stringsAsFactors = FALSE
    )
  }
  
  stats_df <- do.call(rbind, stats_list)
  stats_df$padj <- p.adjust(stats_df$p_value, method = "fdr")
  
  stats_file <- file.path(output_dir, paste0(experiment_name, "_ssGSEA_stats.tsv"))
  write.table(stats_df, file = stats_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("ssGSEA statistics saved to: ", stats_file)
  
  # ----- 6) Boxplots for top pathways -----
  message("Creating ssGSEA boxplots for top pathways...")
  
  # pick top by adjusted p-value
  stats_df <- stats_df[order(stats_df$padj), ]
  n_plot <- min(n_boxplot_pathways, nrow(stats_df))
  top_pw <- stats_df$pathway[seq_len(n_plot)]
  
  # Clean pathway names for plotting
  scores_long$pathway_clean <- gsub("HALLMARK_", "", scores_long$pathway)
  scores_long$pathway_clean <- gsub("_", " ", scores_long$pathway_clean)
  
  stats_df$pathway_clean <- gsub("HALLMARK_", "", stats_df$pathway)
  stats_df$pathway_clean <- gsub("_", " ", stats_df$pathway_clean)
  
  scores_plot <- scores_long[scores_long$pathway %in% top_pw, ]
  scores_plot$pathway_clean <- factor(
    gsub("HALLMARK_", "", scores_plot$pathway),
    levels = gsub("HALLMARK_", "", stats_df$pathway[seq_len(n_plot)])
  )
  
  boxplot_file <- file.path(output_dir, paste0(experiment_name, "_ssGSEA_boxplots.pdf"))
  grDevices::pdf(boxplot_file, width = 11, height = 30)
  
  p <- ggplot2::ggplot(
    scores_plot,
    ggplot2::aes(x = Group, y = score, fill = Group)
  ) +
    ggplot2::geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    ggplot2::geom_jitter(
      width = 0.2, size = 0.8, alpha = 0.5, color = "black"
    ) +
    ggpubr::stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      comparisons = list(c(group1, group2))
    ) +
    ggplot2::facet_wrap(~ pathway_clean, scales = "free_y") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 8),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    ggplot2::labs(
      title = "ssGSEA scores – top pathways",
      x     = NULL,
      y     = "ssGSEA score"
    )
  
  print(p)
  grDevices::dev.off()
  message("ssGSEA boxplots saved to: ", boxplot_file)
  
  message("ssGSEA analysis finished.")
  
  return(list(
    scores   = ssgsea_scores,
    stats    = stats_df,
    pathways = pathways_filtered
  ))
}
