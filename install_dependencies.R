#!/usr/bin/env Rscript

cran_packages <- c(
  "data.table", "dplyr", "future", "ggplot2", "ggtext", "harmony",
  "hdf5r", "Matrix", "openxlsx", "patchwork", "scales", "Seurat",
  "stringr", "tidyr", "writexl"
)

bioconductor_packages <- c(
  "AnnotationDbi", "BiocParallel", "clusterProfiler", "GO.db",
  "scDblFinder", "SingleCellExperiment"
)

cran_missing <- setdiff(cran_packages, rownames(installed.packages()))
if (length(cran_missing)) {
  install.packages(cran_missing, repos = "https://cloud.r-project.org")
}

bioc_missing <- setdiff(
  bioconductor_packages,
  rownames(installed.packages())
)
if (length(bioc_missing)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(bioc_missing, ask = FALSE, update = FALSE)
}
