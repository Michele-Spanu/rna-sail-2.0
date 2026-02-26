#' Run GSEA Analysis
#'
#' Performs Gene Set Enrichment Analysis using fgsea on differential expression results.
#'
#' @param de_results Differential expression results from limma (efit object or data frame)
#' @param species Species for gene sets ("Mus musculus" or "Homo sapiens")
#' @param category MSigDB category (default: "H" for Hallmark)
#' @param subcategory MSigDB subcategory (default: NULL)
#' @param min_size Minimum gene set size (default: 15)
#' @param max_size Maximum gene set size (default: 500)
#' @param n_perm Number of permutations (default: 100000)
#' @return fgsea results data frame
#' @export

run_gsea_analysis <- function(de_results, species = "MM", category = "H",
                              subcategory = NULL, min_size = 15, max_size = 500,
                              n_perm = 100000, customed_pathways = NULL) {
  # !!! species="MM" not "Mus musculus" + added customed_pathways !!!
  # Check required packages
  required_pkgs <- c("fgsea", "msigdbr", "limma")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed")
    }
  }

  # choose correct species string for msigdbr
  if (species == "MM") {
    species1 <- "Mus musculus"
    category <- "MH"
  } else {
    species1 <- "Homo sapiens"
    category <- "H"
  }

  # Get gene sets
  message("Loading gene sets for ", species, ", category: ", category)
  if (is.null(subcategory)) {
    gene_sets_df <- msigdbr::msigdbr(
      db_species  = species,
      species     = species1,
      collection  = category
    )
  } else {
    gene_sets_df <- msigdbr::msigdbr(
      db_species  = species,
      species     = species1,
      collection  = category,
      subcollection = subcategory
    )
  }

  gene_sets <- split(gene_sets_df$ensembl_gene, gene_sets_df$gs_name)
  # !!! Adding the customed pathways !!!
  if (!is.null(customed_pathways)) gene_sets <- c(gene_sets, customed_pathways)
  
  message("Loaded ", length(gene_sets), " gene sets")

  # Prepare rankings
  if (is.data.frame(de_results)) {
    ranks <- de_results$logFC
    names(ranks) <- sub(".*_", "", rownames(de_results))
  } else {
    efit_results <- limma::topTreat(de_results, coef = 1, n = Inf)
    ranks <- efit_results$logFC
    names(ranks) <- sub(".*_", "", rownames(efit_results))
  }

  ranks <- ranks[!is.na(ranks)]
  message("Created rankings for ", length(ranks), " genes")

  # Run GSEA
  message("Running GSEA with ", n_perm, " permutations...")
  set.seed(123)
  gsea_results <- fgsea::fgsea(
    pathways    = gene_sets,
    stats       = ranks,
    minSize     = min_size,
    maxSize     = max_size,
    nPermSimple = n_perm
  )

  gsea_results <- gsea_results[order(-abs(gsea_results$NES)), ]

  # Attach objects needed for enrichment plots
  attr(gsea_results, "gene_sets")  <- gene_sets
  attr(gsea_results, "gene_ranks") <- ranks

  message("GSEA completed. Found ",
          sum(gsea_results$padj < 0.05, na.rm = TRUE),
          " significant pathways (padj < 0.05)")

  return(gsea_results)
}

#' Plot GSEA enrichment curves
#'
#' Create classical GSEA running enrichment plots for top pathways
#' or for user-specified pathways.
#'
#' @param gsea_results Results from fgsea (output of run_gsea_analysis())
#' @param pathways Optional character vector of pathway names to plot
#' @param n_up Number of top upregulated pathways to plot (default: 3)
#' @param n_down Number of top downregulated pathways to plot (default: 3)
#' @param padj_threshold Adjusted p-value cutoff for significance (default: 0.05)
#' @param output_file Optional PDF file to save all plots (multi-page)
#' @param width,height Plot size in inches when saving to file
#' @return Named list of ggplot objects (one per pathway)
#' @export
#' Plot GSEA enrichment curves
#'
#' Create classical GSEA running enrichment plots for top pathways
#' or for user-specified pathways.
#'
#' @param gsea_results Results from fgsea (output of run_gsea_analysis())
#' @param pathways Optional character vector of pathway names to plot
#' @param n_up Number of top upregulated pathways to plot (default: 3)
#' @param n_down Number of top downregulated pathways to plot (default: 3)
#' @param padj_threshold Adjusted p-value cutoff for significance (default: 0.05)
#' @param output_file Optional PDF file to save all plots (multi-page)
#' @param width,height Plot size in inches when saving to file
#' @param gene_sets Optional list of gene sets (if not supplied, taken from gsea_results attributes)
#' @param gene_ranks Optional named vector of gene ranks (if not supplied, taken from gsea_results attributes)
#' @return Named list of ggplot objects (one per pathway)
#' @export
plot_gsea_enrichment <- function(
    gsea_results,
    pathways       = NULL,
    n_up           = 3,
    n_down         = 3,
    padj_threshold = 0.05,
    output_file    = NULL,
    width          = 7,
    height         = 5,
    gene_sets      = NULL,
    gene_ranks     = NULL
) {

  if (!requireNamespace("fgsea", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Packages 'fgsea' and 'ggplot2' are required but not installed")
  }

  # Prefer explicitly passed gene_sets/gene_ranks;
  # if missing, fall back to attributes
  if (is.null(gene_sets)) {
    gene_sets <- attr(gsea_results, "gene_sets")
  }
  if (is.null(gene_ranks)) {
    gene_ranks <- attr(gsea_results, "gene_ranks")
  }

  if (is.null(gene_sets) || is.null(gene_ranks)) {
    stop("gene_sets and/or gene_ranks not found.\n",
         "Pass them explicitly or re-run run_gsea_analysis() so they are stored as attributes.")
  }

  # If user didn't specify pathways, choose top significant up/down
  if (is.null(pathways)) {
    sig <- gsea_results[!is.na(gsea_results$padj) &
                          gsea_results$padj < padj_threshold, ]

    up   <- sig[sig$NES > 0, , drop = FALSE]
    down <- sig[sig$NES < 0, , drop = FALSE]

    up   <- up[order(-up$NES), , drop = FALSE]
    down <- down[order(down$NES), , drop = FALSE]  # NES is negative, so ascending

    sel_up   <- head(up$pathway,   n_up)
    sel_down <- head(down$pathway, n_down)

    pathways <- c(sel_up, sel_down)
  }

  pathways <- intersect(pathways, gsea_results$pathway)
  if (length(pathways) == 0L) {
    stop("No pathways to plot (after filtering). Check 'pathways' names.")
  }

  plots <- list()

  if (!is.null(output_file)) {
    grDevices::pdf(output_file, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  for (pw in pathways) {
    genes <- gene_sets[[pw]]
    if (is.null(genes)) next

    p <- fgsea::plotEnrichment(genes, gene_ranks)

    # Get NES / padj for title
    row <- gsea_results[gsea_results$pathway == pw, ][1, ]
    pretty_name <- gsub("HALLMARK_", "", pw)
    pretty_name <- gsub("_", " ", pretty_name)

    p <- p +
      ggplot2::labs(
        title = paste0(pretty_name,
                       "\nNES = ", signif(row$NES, 3),
                       ", padj = ", signif(row$padj, 3))
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )

    plots[[pw]] <- p

    if (!is.null(output_file)) {
      print(p)
    }
  }

  if (!is.null(output_file)) {
    message("GSEA enrichment plots saved to: ", output_file)
  }

  return(plots)
}


#' Create GSEA Barplot
#'
#' Creates a barplot visualization of top GSEA results.
#'
#' @param gsea_results Results from fgsea
#' @param n_pathways Number of top pathways to show (default: 20)
#' @param output_file Output PDF file path (optional)
#' @param width Plot width in inches (default: 8)
#' @param height Plot height in inches (default: 6)
#' @return ggplot object
#' @export
plot_gsea_barplot <- function(
    gsea_results,
    n_pathways     = 20,
    output_file    = NULL,
    width          = 8,
    height         = 6,
    padj_threshold = 0.05,
    color_up       = "#E31A1C",
    color_down     = "#1F78B4",
    color_ns       = "grey70"
) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed")
  }

  # --- CONSISTENT selection: Top by absolute NES (same as dotplot) ---
  plot_data <- gsea_results[order(-abs(gsea_results$NES)), , drop = FALSE]
  plot_data <- head(plot_data, min(n_pathways, nrow(plot_data)))
  plot_data <- plot_data[!is.na(plot_data$pathway) & !is.na(plot_data$NES), , drop = FALSE]

  # Clean pathway names
  plot_data$pathway_clean <- gsub("HALLMARK_", "", plot_data$pathway)
  plot_data$pathway_clean <- gsub("_", " ", plot_data$pathway_clean)

  # Direction + significance
  plot_data$category <- "Not Significant"
  plot_data$category[!is.na(plot_data$padj) & plot_data$padj < padj_threshold & plot_data$NES > 0] <- "Upregulated"
  plot_data$category[!is.na(plot_data$padj) & plot_data$padj < padj_threshold & plot_data$NES < 0] <- "Downregulated"

  # Order by NES for plotting (nice left-to-right / top-to-bottom)
  plot_data <- plot_data[order(plot_data$NES), , drop = FALSE]
  plot_data$pathway_clean <- factor(plot_data$pathway_clean, levels = plot_data$pathway_clean)

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = pathway_clean, y = NES, fill = category)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = c(
        "Upregulated"     = color_up,
        "Downregulated"   = color_down,
        "Not Significant" = color_ns
      ),
      name = "Category"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(
      title = "GSEA Results - Top Enriched Pathways",
      x     = "Pathway",
      y     = "Normalized Enrichment Score (NES)"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10),
      legend.position = "bottom"
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = p, width = width, height = height)
    message("GSEA barplot saved to: ", output_file)
  }

  return(p)
}

#' Create GSEA Dotplot
#'
#' Creates a dotplot visualization of GSEA results showing NES, significance, and gene set size.
#'
#' @param gsea_results Results from fgsea
#' @param n_pathways Number of top pathways to show (default: 30)
#' @param output_file Output PDF file path (optional)
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
#' @return ggplot object
#' @export
plot_gsea_dotplot <- function(
    gsea_results,
    n_pathways     = 30,
    output_file    = NULL,
    width          = 10,
    height         = 8,
    padj_threshold = 0.05,
    color_up       = "#E31A1C",
    color_down     = "#1F78B4",
    color_ns       = "grey70"
) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed")
  }

  # --- CONSISTENT selection: Top by absolute NES (same as barplot) ---
  plot_data <- gsea_results[order(-abs(gsea_results$NES)), , drop = FALSE]
  plot_data <- head(plot_data, min(n_pathways, nrow(plot_data)))
  plot_data <- plot_data[!is.na(plot_data$pathway) & !is.na(plot_data$NES), , drop = FALSE]

  # Clean names
  plot_data$pathway_clean <- gsub("HALLMARK_", "", plot_data$pathway)
  plot_data$pathway_clean <- gsub("_", " ", plot_data$pathway_clean)

  # Order by NES
  plot_data <- plot_data[order(plot_data$NES), , drop = FALSE]
  plot_data$pathway_clean <- factor(plot_data$pathway_clean,
                                    levels = plot_data$pathway_clean)

  # Category
  plot_data$category <- "Not Significant"
  plot_data$category[!is.na(plot_data$padj) & plot_data$padj < padj_threshold & plot_data$NES > 0] <- "Upregulated"
  plot_data$category[!is.na(plot_data$padj) & plot_data$padj < padj_threshold & plot_data$NES < 0] <- "Downregulated"

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = NES, y = pathway_clean)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(size = size, color = category),
      alpha = 0.9
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Upregulated"     = color_up,
        "Downregulated"   = color_down,
        "Not Significant" = color_ns
      ),
      name = "Category"
    ) +
    ggplot2::scale_size_continuous(name = "Gene set size", range = c(2, 8)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(
      title = "GSEA Results - Pathway Enrichment",
      x     = "Normalized Enrichment Score (NES)",
      y     = NULL
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      axis.text.y = ggplot2::element_text(size = 9)
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = p, width = width, height = height)
    message("GSEA dotplot saved to: ", output_file)
  }

  return(p)
}


#' Create GSEA Table Plot
#'
#' Creates a table-style visualization of GSEA results using fgsea's plotGseaTable.
#'
#' @param gsea_results Results from fgsea
#' @param gene_sets Gene sets list used for GSEA
#' @param gene_ranks Named vector of gene rankings
#' @param n_pathways Number of pathways to show (default: 30)
#' @param output_file Output PDF file path (optional)
#' @param width Plot width in inches (default: 11)
#' @param height Plot height in inches (default: 15)
#' @return None (creates plot)
#' @export
create_gsea_table_plot <- function(gsea_results, gene_sets, gene_ranks,
                                   n_pathways = 30,
                                   output_file = NULL,
                                   width = 11, height = 15) {

  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required but not installed")
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Package 'grid' is required but not installed")
  }

  # Select top pathways
  top_up   <- gsea_results[gsea_results$NES > 0, ]
  top_up   <- top_up[order(-top_up$NES), ]

  top_down <- gsea_results[gsea_results$NES < 0, ]
  top_down <- top_down[order(top_down$NES), ]

  n_up   <- min(n_pathways / 2, nrow(top_up))
  n_down <- min(n_pathways / 2, nrow(top_down))

  top_pathways_up   <- head(top_up$pathway,   n_up)
  top_pathways_down <- head(top_down$pathway, n_down)
  top_pathways      <- c(top_pathways_up, rev(top_pathways_down))

  # If nothing to plot, bail out early
  if (length(top_pathways) == 0L) {
    warning("No pathways selected for plotGseaTable (check NES values).")
    return(invisible(NULL))
  }

  # Open device if requested
  if (!is.null(output_file)) {
    grDevices::pdf(output_file, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  # Create grob
  p <- fgsea::plotGseaTable(
    pathways  = gene_sets[top_pathways],
    stats     = gene_ranks,
    fgseaRes  = gsea_results,
    gseaParam = 0.5
  )

  # Actually draw it
  grid::grid.newpage()
  grid::grid.draw(p)

  if (!is.null(output_file)) {
    message("GSEA table plot saved to: ", output_file)
  }

  invisible(p)
}

#' Run Camera Gene Set Testing
#'
#' Performs competitive gene set testing using limma's camera method.
#'
#' @param voom_object Voom-transformed expression object
#' @param design_matrix Design matrix
#' @param contrast_vector Contrast vector
#' @param gene_sets List of gene sets (indices)
#' @param species Species for gene sets ("Mus musculus" or "Homo sapiens")
#' @param min_size Minimum gene set size (default: 15)
#' @return Camera results data frame
#' @export
run_camera_analysis <- function(voom_object, design_matrix, contrast_vector, gene_sets = NULL,
                               species = "Mus musculus", min_size = 15) {

  if (!requireNamespace("limma", quietly = TRUE) ||
      !requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Required packages (limma, msigdbr) not installed")
  }

  if (species == "MM") {
    species1 <- "Mus musculus"
    category <- "MH"
  } else {
    species1 <- "Homo sapiens"
    category <- "H"
  }

  # If gene_sets not provided, use Hallmark sets
  if (is.null(gene_sets)) {
    message("Loading Hallmark gene sets for ", species)
    h_gene_sets <- msigdbr::msigdbr(db_species = species, species = species1, collection  = "H")
    gene_sets_list <- split(h_gene_sets$ensembl_gene, h_gene_sets$gs_name)

    # Extract ENSEMBL IDs from rownames
    rownames_voom <- rownames(voom_object$E)
    ensembl_ids <- sub(".*_(ENSMUSG[0-9]+).*", "\\1", rownames_voom)

    # Convert to indices
    gene_sets_indices <- lapply(gene_sets_list, function(gs) {
      which(!is.na(ensembl_ids) & ensembl_ids %in% gs)
    })

    # Filter by minimum size
    gene_sets_indices <- gene_sets_indices[sapply(gene_sets_indices, length) >= min_size]
  } else {
    gene_sets_indices <- gene_sets
  }

  message("Running camera analysis with ", length(gene_sets_indices), " gene sets")

  # Run camera
  camera_results <- limma::camera(
    y = voom_object$E,
    index = gene_sets_indices,
    design = design_matrix,
    contrast = contrast_vector
  )

  camera_results_df <- as.data.frame(camera_results)
  camera_results_df$pathway <- rownames(camera_results_df)

  message("Camera analysis completed. Found ", sum(camera_results_df$PValue < 0.05),
          " significant pathways (p < 0.05)")

  return(camera_results_df)
}

#' Save Pathway Analysis Results
#'
#' Saves pathway analysis results to files with proper formatting.
#'
#' @param gsea_results GSEA results data frame
#' @param camera_results Camera results data frame (optional)
#' @param experiment_name Name of the experiment
#' @param output_dir Output directory
#' @return None (invisible)
#' @export
save_pathway_results <- function(gsea_results, camera_results = NULL, experiment_name, output_dir) {

  create_output_dir(output_dir)

  # Save GSEA results
  gsea_file <- file.path(output_dir, paste0(experiment_name, "_GSEA_results.tsv"))
  gsea_results$leadingEdge <- sapply(gsea_results$leadingEdge, function(x) paste(x, collapse = ";"))
  write.table(gsea_results, file = gsea_file, sep = "\t", quote = FALSE, row.names = FALSE)

  message("GSEA results saved to: ", gsea_file)

  # Save Camera results if provided
  if (!is.null(camera_results)) {
    camera_file <- file.path(output_dir, paste0(experiment_name, "_Camera_results.tsv"))
    write.table(camera_results, file = camera_file, sep = "\t", quote = FALSE, row.names = TRUE)
    message("Camera results saved to: ", camera_file)
  }

  invisible()
}

# !!! Remember to check package existance: readxl and biomaRt !!!
retrieve_pathway <- function(file.path) {
  pathway_list <- list()
  
  pathway_names <- readx::excel_sheets(file.path, species)
  for (pathway in pathway_names) {
      # File must contain a col with the NCBI Entrez ID (named "NCBI") or one with the ENSEMBL IDs (named "ENSEMBL") or one containing the gene Symbol (named "Symbol")  
    excel.df <- readxl::read_excel(file.path, sheet=pathway)
  
    if ("ENSEMBL" %in% colnames(excel.df)) pathway_list[[pathway]] <- unique(excel.df$ENSEMBL) else {
      ensembl <- biomaRt::useEnsembl(biomart = "genes", 
                                     dataset = if (species=="mouse")  "mmusculus_gene_ensembl" else "hsapiens_gene_symbol")
      if ("NCBI" %in% colnames(excel.df)) {
        gene_IDs <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id", "external_gene_name"),
                          filters = "entrezgene_id", values = unique(excel.df$NCBI), 
                          mart = ensembl, useCache = FALSE)
      } else {
        if ("Symbol" %in% colnames(excel.df)) {
          gene_IDs <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
                            filters = "hgnc_symbol", values = unique(excel.df$Symbol),
                            mart = ensembl, useCache = TRUE)
        } else  stop("There is no column with NCBI, ENSEMBL, or Symbol ID in sheet:", pathway)
      }

      pathway_list[[pathway]] <- unique(gene_IDs$ensembl_gene_id)
    }
  }

return(pathway_list)
}
