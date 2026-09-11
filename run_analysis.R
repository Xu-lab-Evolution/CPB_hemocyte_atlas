#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg)) {
  setwd(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))))
}

stages <- c(
  "R/01_preprocess_and_cluster.R",
  "R/02_cell_cycle.R",
  "R/03_cluster_validation.R",
  "R/04_replicate_consensus_markers_GO.R",
  "R/05_hemocyte_annotation.R",
  "R/06_main_figures.R",
  "R/07_GO_and_annotation_markers.R"
)

for (stage in stages) {
  message("Running ", stage)
  source(stage, local = new.env(parent = globalenv()), chdir = FALSE)
}

if (tolower(Sys.getenv("EXPORT_SHINYCELL2", unset = "false")) == "true") {
  message("Running R/08_export_shinycell2.R")
  source("R/08_export_shinycell2.R", local = new.env(parent = globalenv()))
}

message("Analysis completed.")
