#' Build limma design matrix and contrast
#' @keywords internal
build_design_and_contrast <- function(metadata,
                                      condition_column,
                                      group1_condition,
                                      group2_condition,
                                      covariates = NULL,
                                      design_formula = NULL,
                                      contrast_string = NULL) {
  
  if (!condition_column %in% colnames(metadata)) {
    stop("condition_column '", condition_column, "' not found in metadata.")
  }
  
  # Validate covariates exist
  if (!is.null(covariates)) {
    if (!all(covariates %in% colnames(metadata))) {
      missing <- covariates[!covariates %in% colnames(metadata)]
      stop("Covariate column(s) not found in metadata: ", paste(missing, collapse = ", "))
    }
  }
  
  # Coerce to factors (recommended for limma; prevents surprises)
  metadata[[condition_column]] <- factor(metadata[[condition_column]],
                                         levels = c(group2_condition, group1_condition))
  
  if (!is.null(covariates)) {
    for (cv in covariates) {
      metadata[[cv]] <- factor(metadata[[cv]])
      if (nlevels(metadata[[cv]]) < 2) {
        stop("Covariate '", cv, "' has <2 levels in the subsetted data.")
      }
    }
  }
  
  # Decide formula
  if (!is.null(design_formula)) {
    # allow user to pass either "~ ..." or "..." (we'll normalize)
    ftxt <- trimws(design_formula)
    if (!startsWith(ftxt, "~")) ftxt <- paste("~", ftxt)
    design_f <- stats::as.formula(ftxt)
    
    # must include the condition term somewhere
    if (!grepl(paste0("\\b", condition_column, "\\b"), ftxt)) {
      stop("design_formula must include the condition column '", condition_column, "'.\n",
           "Example: '~ cell_line + ", condition_column, "'")
    }
    
    # If interaction present, we cannot auto-pick a single coefficient safely
    has_interaction <- grepl(":", ftxt) || grepl("\\*", ftxt)
    if (has_interaction && is.null(contrast_string)) {
      stop("design_formula contains interaction terms (':' or '*'). ",
           "Please provide contrast_string explicitly for limma::makeContrasts().")
    }
    
  } else {
    # default guided formula
    if (is.null(covariates) || length(covariates) == 0) {
      design_f <- stats::as.formula(paste0("~ ", condition_column))
    } else {
      design_f <- stats::as.formula(paste0("~ ", paste(covariates, collapse = " + "),
                                           " + ", condition_column))
    }
  }
  
  design <- stats::model.matrix(design_f, data = metadata)
  
  # Contrast
  if (!is.null(contrast_string)) {
    contrast_matrix <- limma::makeContrasts(contrasts = contrast_string, levels = design)
  } else {
    # With intercept design (~ cov + condition), and condition ref = group2,
    # the coefficient name is usually like "conditionmut"
    coef_raw <- paste0(condition_column, group1_condition)
    coef_name <- make.names(coef_raw)
    
    if (!coef_name %in% colnames(design)) {
      stop("Cannot auto-build contrast. Expected coefficient '", coef_name,
           "' not found in design matrix.\n",
           "Available coefficients: ", paste(colnames(design), collapse = ", "), "\n",
           "Provide contrast_string explicitly (e.g. '", coef_name, "').")
    }
    
    contrast_matrix <- limma::makeContrasts(contrasts = coef_name, levels = design)
  }
  
  list(design = design, contrast_matrix = contrast_matrix, design_formula = design_f)
}
