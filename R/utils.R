#' RNASAIL: Comprehensive RNA-seq Analysis Pipeline
#'
#' This package provides functions for differential expression analysis,
#' pathway analysis, WGCNA co-expression networks, and visualization of RNA-seq data.
#'
#' @docType package
#' @name rna-sail

# Package startup message
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("RNASAIL loaded successfully!")
  packageStartupMessage("Use ?RNASAIL to get started.")

  # Set up conflicts preferences
  if (requireNamespace("conflicted", quietly = TRUE)) {
    conflicted::conflicts_prefer(stats::var)
    conflicted::conflicts_prefer(dplyr::select)
  }
}

#' Install Missing Packages
#'
#' Helper function to install missing CRAN packages.
#'
#' @param packages Character vector of package names to check and install
#' @return None (invisible)
#' @export
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
  invisible()
}

#' Install Missing Bioconductor Packages
#'
#' Helper function to install missing Bioconductor packages.
#'
#' @param packages Character vector of Bioconductor package names
#' @return None (invisible)
#' @export
bioc_install_if_missing <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
    library(pkg, character.only = TRUE)
  }
  invisible()
}

#' Check Column Validity for PCA
#'
#' Helper function to check if a column has finite values and non-zero variance.
#'
#' @param x Numeric vector to check
#' @return Logical indicating if column is valid
#' @keywords internal
is_good_column <- function(x) {
  all(is.finite(x)) && var(x) > 0
}

#' Create Output Directory
#'
#' Creates output directory if it doesn't exist.
#'
#' @param dir_path Path to directory to create
#' @return Character path to directory
#' @export
create_output_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Created directory: ", dir_path)
  }
  return(dir_path)
}

#' Get Default Required Packages
#'
#' Returns the default list of required packages for the pipeline.
#'
#' @return Character vector of package names
#' @export
get_required_packages <- function() {
  c(
    "rtracklayer", "limma", "edgeR", "fgsea", "msigdbr",
    "dplyr", "ggplot2", "ggrepel", "RColorBrewer", "colorspace",
    "ComplexHeatmap", "circlize", "data.table", "stringr",
    "conflicted", "scales", "tidyverse", "matrixStats"
  )
}

#' Get Default Bioconductor Packages
#'
#' Returns the default list of required Bioconductor packages.
#'
#' @return Character vector of Bioconductor package names
#' @export
get_required_bioc_packages <- function() {
  c(
    "limma", "edgeR", "msigdbr", "DOSE", "clusterProfiler",
    "WGCNA", "viper", "dorothea", "rtracklayer", "org.Mm.eg.db",
    "org.Hs.eg.db", "biomaRt"
  )
}
