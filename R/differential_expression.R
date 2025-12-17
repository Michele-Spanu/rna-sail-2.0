#' Run Limma-Voom Differential Expression Analysis
#'
#' Performs differential expression analysis using the limma-voom pipeline.
#'
#' @param expr_data Expression count data (genes x samples)
#' @param sample_groups Character vector of sample group labels
#' @param contrast_matrix Contrast matrix from makeContrasts
#' @param design_matrix Design matrix from model.matrix
#' @param experiment_name Name for the experiment (used in output files)
#' @param output_dir Output directory path
#' @return List containing tfit (treat results) and efit (eBayes results)
#' @export
run_limma_voom <- function(expr_data, sample_groups, contrast_matrix, design_matrix,
                           experiment_name, output_dir,
                           lfc_threshold = 1,
                           fdr_threshold = 0.05) {

  # Check required packages
  required_pkgs <- c("edgeR", "limma")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed")
    }
  }

  # Create output directory
  create_output_dir(output_dir)


  # EdgeR filtering and normalization
  d0 <- edgeR::DGEList(expr_data)
  keep.exprs <- edgeR::filterByExpr(d0, group = sample_groups)
  d0 <- d0[keep.exprs, ]
  d0 <- edgeR::calcNormFactors(d0)

  # Voom transformation
  y <- limma::voom(d0, design_matrix)

  # Linear model fitting
  vfit <- limma::lmFit(y, design_matrix)
  vfit <- limma::contrasts.fit(vfit, contrasts = contrast_matrix)

  # eBayes and treat
  efit <- limma::eBayes(vfit)
  tfit <- limma::treat(vfit, lfc = lfc_threshold)

  # ---- Full DE tables for BOTH efit and tfit ----
  results_efit <- limma::topTable(
    efit,
    coef   = 1,
    number = Inf,
    sort.by = "none"
  )

  results_tfit <- limma::topTreat(
    tfit,
    coef = 1,
    n    = Inf
  )

  # Save both tables
  efit_file  <- file.path(output_dir, paste0(experiment_name, "_DE_efit_results.tsv"))
  tfit_file  <- file.path(output_dir, paste0(experiment_name, "_DE_tfit_results.tsv"))

  write.table(results_efit, file = efit_file, sep = "\t", quote = FALSE, row.names = TRUE)
  write.table(results_tfit, file = tfit_file, sep = "\t", quote = FALSE, row.names = TRUE)

  n_sig_tfit <- sum(results_tfit$adj.P.Val < fdr_threshold, na.rm = TRUE)
  message("Differential expression (treat) results saved to: ", tfit_file)
  message("Differential expression (eBayes) results saved to: ", efit_file)
  message("Found ", n_sig_tfit, " significant genes in tfit (adj.P.Val < ", fdr_threshold, ")")

  return(list(
    tfit          = tfit,
    efit          = efit,
    results_tfit  = results_tfit,
    results_efit  = results_efit,
    results       = results_tfit,
    lfc_threshold = lfc_threshold,
    fdr_threshold = fdr_threshold
  ))
}


#' Run Complete Differential Expression Pipeline
#'
#' Runs the complete differential expression analysis pipeline for two groups,
#' optionally adjusting for covariates and/or using a user-defined design formula.
#'
#' @param counts_data Count expression data (genes x samples)
#' @param metadata Sample metadata data frame
#' @param group1_condition Name of condition 1 (treatment group)
#' @param group2_condition Name of condition 2 (control group)
#' @param condition_column Column name in metadata containing condition information (default: "condition")
#' @param sample_id_column Column name in metadata containing sample IDs (default: "SampleID")
#' @param experiment_name Name for the experiment
#' @param output_dir Output directory
#' @param fdr_threshold FDR threshold for significance (default: 0.1)
#' @param lfc_threshold Log fold change threshold (default: 1)
#' @param covariates Optional character vector of covariate columns in metadata (e.g. c("cell_line","batch"))
#' @param design_formula Optional design formula string (advanced), e.g. "~ cell_line + condition".
#'        Must include condition_column. If it contains ':' or '*', provide contrast_string.
#' @param contrast_string Optional limma contrast string. If NULL and design has no interactions,
#'        a default group1 vs group2 contrast is used.
#' @return List containing DE results and summary statistics
#' @export
run_differential_expression <- function(counts_data, metadata, group1_condition, group2_condition,
                                        condition_column = "condition",
                                        sample_id_column = "SampleID",
                                        experiment_name,
                                        output_dir,
                                        fdr_threshold = 0.1,
                                        lfc_threshold = 1,
                                        covariates = NULL,
                                        design_formula = NULL,
                                        contrast_string = NULL) {

  message("Starting differential expression analysis: ", experiment_name)

  # Packages needed here
  required_pkgs <- c("limma")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed")
    }
  }

  # Match metadata to expression data
  metadata_ordered <- match_metadata_to_expression(counts_data, metadata, sample_id_column)

  # Get sample names for each group (keep your old behaviour)
  group1_samples <- metadata_ordered[[sample_id_column]][metadata_ordered[[condition_column]] == group1_condition]
  group2_samples <- metadata_ordered[[sample_id_column]][metadata_ordered[[condition_column]] == group2_condition]

  group1_samples_clean <- paste0(group1_samples)
  group2_samples_clean <- paste0(group2_samples)

  # Check if samples exist in data
  all_samples <- c(group1_samples_clean, group2_samples_clean)
  missing_samples <- setdiff(all_samples, colnames(counts_data))
  if (length(missing_samples) > 0) {
    warning("Missing samples in expression data: ", paste(missing_samples, collapse = ", "))
  }

  # Subset expression data
  available_samples <- intersect(all_samples, colnames(counts_data))
  if (length(available_samples) < 2) {
    stop("Fewer than 2 samples available after subsetting. Check SampleID matching.")
  }
  expr_subset <- counts_data[, available_samples, drop = FALSE]

  # Subset & reorder metadata to match expression cols
  metadata_sub <- metadata_ordered[metadata_ordered[[sample_id_column]] %in% available_samples, , drop = FALSE]
  metadata_sub <- metadata_sub[match(colnames(expr_subset), metadata_sub[[sample_id_column]]), , drop = FALSE]

  # Create group labels (kept for filterByExpr)
  sample_groups <- ifelse(colnames(expr_subset) %in% group1_samples_clean,
                          group1_condition, group2_condition)
  sample_groups <- factor(sample_groups, levels = c(group2_condition, group1_condition))

  # ---- NEW DESIGN HANDLING ----
  # Validate covariates (if used)
  if (!is.null(covariates)) {
    if (!is.character(covariates)) {
      stop("covariates must be a character vector of metadata column names (e.g. c('cell_line','batch')).")
    }
    missing_cov <- setdiff(covariates, colnames(metadata_sub))
    if (length(missing_cov) > 0) {
      stop("Covariate column(s) not found in metadata: ", paste(missing_cov, collapse = ", "))
    }
  }

  # Ensure condition is a factor with group2 as reference (important!)
  metadata_sub[[condition_column]] <- factor(
    metadata_sub[[condition_column]],
    levels = c(group2_condition, group1_condition)
  )

  # Coerce covariates to factors (recommended)
  if (!is.null(covariates) && length(covariates) > 0) {
    for (cv in covariates) {
      metadata_sub[[cv]] <- factor(metadata_sub[[cv]])
      if (nlevels(metadata_sub[[cv]]) < 2) {
        stop("Covariate '", cv, "' has <2 levels after subsetting to the two groups.")
      }
    }
  }

  # Decide design formula
  if (!is.null(design_formula)) {
    ftxt <- trimws(design_formula)
    if (!startsWith(ftxt, "~")) ftxt <- paste("~", ftxt)

    if (!grepl(paste0("\\b", condition_column, "\\b"), ftxt)) {
      stop("design_formula must include the condition column '", condition_column, "'. ",
           "Example: '~ cell_line + ", condition_column, "'.")
    }

    has_interaction <- grepl(":", ftxt) || grepl("\\*", ftxt)
    if (has_interaction && is.null(contrast_string)) {
      stop("design_formula contains interaction terms (':' or '*'). ",
           "Please provide contrast_string explicitly.")
    }

    design_f <- stats::as.formula(ftxt)

  } else {
    if (is.null(covariates) || length(covariates) == 0) {
      design_f <- stats::as.formula(paste0("~ ", condition_column))
    } else {
      design_f <- stats::as.formula(
        paste0("~ ", paste(covariates, collapse = " + "), " + ", condition_column)
      )
    }
  }

  design <- stats::model.matrix(design_f, data = metadata_sub)

  # Decide contrast
  if (!is.null(contrast_string)) {
    contrast_matrix <- limma::makeContrasts(contrasts = contrast_string, levels = design)
  } else {
    # Auto-contrast: for intercept design, treated-vs-control is coefficient condition<group1>
    coef_raw  <- paste0(condition_column, group1_condition)
    coef_name <- make.names(coef_raw)

    if (!coef_name %in% colnames(design)) {
      stop(
        "Auto contrast failed. Expected coefficient '", coef_name, "' not found in design.\n",
        "Design columns are: ", paste(colnames(design), collapse = ", "), "\n",
        "Provide contrast_string explicitly if needed."
      )
    }

    contrast_matrix <- limma::makeContrasts(contrasts = coef_name, levels = design)
  }
  print(contrast_matrix)
  message("Using design: ", paste(deparse(design_f), collapse = " "))
  message("Design columns: ", paste(colnames(design), collapse = ", "))
  message("Contrast columns: ", paste(colnames(contrast_matrix), collapse = ", "))

  # Run limma-voom
  de_results <- run_limma_voom(
    expr_data       = expr_subset,
    sample_groups   = sample_groups,
    contrast_matrix = contrast_matrix,
    design_matrix   = design,
    experiment_name = experiment_name,
    output_dir      = output_dir,
    fdr_threshold   = fdr_threshold,
    lfc_threshold   = lfc_threshold
  )

  # Summaries for both efit and tfit
  summary_efit <- summarize_de_results(
    de_results$efit,
    fdr_threshold = fdr_threshold,
    lfc_threshold = lfc_threshold,
    fit_label     = "efit"
  )

  summary_tfit <- summarize_de_results(
    de_results$tfit,
    fdr_threshold = fdr_threshold,
    lfc_threshold = lfc_threshold,
    fit_label     = "tfit"
  )

  # Save gene lists for both fits
  save_gene_lists(
    de_results$efit,
    experiment_name = experiment_name,
    output_dir      = output_dir,
    fdr_threshold   = fdr_threshold,
    lfc_threshold   = lfc_threshold,
    fit_label       = "efit"
  )

  save_gene_lists(
    de_results$tfit,
    experiment_name = experiment_name,
    output_dir      = output_dir,
    fdr_threshold   = fdr_threshold,
    lfc_threshold   = lfc_threshold,
    fit_label       = "tfit"
  )

  message("Differential expression analysis completed!")

  return(list(
    de_results   = de_results,
    summary_efit = summary_efit,
    summary_tfit = summary_tfit,
    sample_info  = list(
      group1           = group1_samples,
      group2           = group2_samples,
      group1_condition = group1_condition,
      group2_condition = group2_condition,
      covariates       = covariates,
      design_formula   = if (is.null(design_formula)) NULL else design_formula,
      contrast_string  = contrast_string
    )
  ))
}


#' Summarize Differential Expression Results
#'
#' @keywords internal
summarize_de_results <- function(fit,
                                 fdr_threshold = 0.1,
                                 lfc_threshold = 1,
                                 fit_label = "efit") {

  if (inherits(fit, "MArrayLM2")) {  # treat object
    results <- limma::decideTests(fit, p.value = fdr_threshold)
  } else {
    results <- limma::decideTests(fit, p.value = fdr_threshold, lfc = lfc_threshold)
  }

  n_up    <- sum(results == 1)
  n_down  <- sum(results == -1)
  n_total <- nrow(fit$coefficients)

  summary_stats <- list(
    total_genes         = n_total,
    significant_genes   = n_up + n_down,
    upregulated         = n_up,
    downregulated       = n_down,
    percent_significant = round((n_up + n_down) / n_total * 100, 2),
    fdr_threshold       = fdr_threshold,
    lfc_threshold       = lfc_threshold,
    fit_label           = fit_label
  )

  message("[", fit_label, "] Summary: ", n_up, " up, ", n_down,
          " down out of ", n_total, " genes")

  return(summary_stats)
}


#' Save Gene Lists
#'
#' @keywords internal
save_gene_lists <- function(fit, experiment_name, output_dir,
                            fdr_threshold, lfc_threshold,
                            fit_label = "efit") {

  if (inherits(fit, "MArrayLM2")) {  # treat
    results <- limma::decideTests(fit, p.value = fdr_threshold)
  } else {
    results <- limma::decideTests(fit, p.value = fdr_threshold, lfc = lfc_threshold)
  }

  up_genes   <- rownames(results)[results == 1]
  down_genes <- rownames(results)[results == -1]

  up_file   <- file.path(output_dir, paste0(experiment_name, "_UP_",   fit_label, ".tsv"))
  down_file <- file.path(output_dir, paste0(experiment_name, "_DOWN_", fit_label, ".tsv"))

  write.table(up_genes,   file = up_file,   quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(down_genes, file = down_file, quote = FALSE, row.names = FALSE, col.names = FALSE)

  message("[", fit_label, "] Gene lists saved: ",
          basename(up_file), ", ", basename(down_file))

  invisible()
}
