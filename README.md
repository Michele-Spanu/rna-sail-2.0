# RNASAIL

RNASAIL is an R package providing a full RNA‑seq analysis workflow, including differential expression, pathway analysis, transcription factor activity, immune infiltration, WGCNA, and LINCS connectivity scoring.

## Installation

You can install RNASAIL via devtools:

## Install directly from GitHub

```r
install.packages("devtools")
devtools::install_github("OncologyMedunigraz/rna-sail")
library(RNASAIL)
```



# Verify Installation

```r
library(RNASAIL)

rna-sail::get_required_packages()

## Required Input Files

### 1. Counts matrix
Tab-separated file: genes x samples, raw counts.

### 2. Metadata
Must include: SampleID and a condition column.

### 3. Annotation GTF (optional)
Used for filtering protein‑coding genes.

## Quick Start

```r
library(rna-sail)

counts <- "counts.tsv"
tpm <- "tpm.tsv"
metadata <- "metadata.tsv"
gtf <- "Mus_musculus.GRCm39.111.gtf" # or gencode.v46.annotation.gtf for human

results <- run_complete_pipeline(
  counts_file = counts,
  tpm_file = tpm,
  metadata_file = metadata,
  group1_condition = "treated",
  group2_condition = "control",
  annotation_gtf = gtf,
  output_dir = "Results/",
  experiment_name = "treated_vs_control",
  fdr_threshold=0.05,
  lfc_threshold=1
)
```

## Key Functions
- run_differential_expression
- run_gsea_analysis
- run_tf_analysis
- run_immune_deconvolution
- run_wgcna_analysis
- run_complete_pipeline

## Vignette
```r
vignette("RNASAIL")
```

## Issues
Submit issues at:
https://github.com/OncologyMedunigraz/rna-sail/issues

