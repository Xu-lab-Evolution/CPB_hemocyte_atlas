# Cell-cycle gene-set scoring for the selected CPB hemocyte clusters.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
})

source(file.path("R", "config.R"))

ANNOTATION_DIR <- file.path(OUT, "annotation_transfer")
CELL_CYCLE_DIR <- file.path(OUT, "cell_cycle")
dir.create(CELL_CYCLE_DIR, recursive = TRUE, showWarnings = FALSE)

OBJECT_FILE <- file.path(OUT, "combined_selected_resolution.rds")
MAPPING_FILE <- file.path(
  ANNOTATION_DIR, "cell_cycle_old_to_v5_mapping.tsv"
)
require_files(c(OBJECT_FILE, MAPPING_FILE))

set.seed(SEED)
combined <- JoinLayers(readRDS(OBJECT_FILE))
DefaultAssay(combined) <- "RNA"
mapping <- read.delim(MAPPING_FILE, check.names = FALSE)

required_metadata <- c("replicate", "res_0.2")
missing_metadata <- setdiff(required_metadata, colnames(combined[[]]))
if (length(missing_metadata)) {
  stop("Missing metadata: ", paste(missing_metadata, collapse = ", "),
       call. = FALSE)
}

cycle_gene_table <- mapping |>
  filter(mapping_status == "mapped_and_detected") |>
  distinct(phase, gene)
cycle_sets <- split(cycle_gene_table$gene, cycle_gene_table$phase)

for (phase in c("S", "G2/M")) {
  if (length(cycle_sets[[phase]]) < 2L) {
    stop("Too few mapped genes for the ", phase, " score.", call. = FALSE)
  }
}

combined <- AddModuleScore(
  combined, features = list(cycle_sets[["S"]]), name = "cc_s_tmp",
  assay = "RNA", seed = SEED
)
combined <- AddModuleScore(
  combined, features = list(cycle_sets[["G2/M"]]), name = "cc_g2m_tmp",
  assay = "RNA", seed = SEED
)
combined$cell_cycle_s_score <- combined$cc_s_tmp1
combined$cell_cycle_g2m_score <- combined$cc_g2m_tmp1
combined$cc_s_tmp1 <- NULL
combined$cc_g2m_tmp1 <- NULL

s_cutoff <- unname(quantile(combined$cell_cycle_s_score, 0.75, na.rm = TRUE))
g2m_cutoff <- unname(quantile(
  combined$cell_cycle_g2m_score, 0.75, na.rm = TRUE
))
combined$cell_cycle_s_high <- combined$cell_cycle_s_score > s_cutoff
combined$cell_cycle_g2m_high <- combined$cell_cycle_g2m_score > g2m_cutoff
combined$cell_cycle_signal <- case_when(
  combined$cell_cycle_s_high & combined$cell_cycle_g2m_high ~ "S and G2/M high",
  combined$cell_cycle_s_high ~ "S high only",
  combined$cell_cycle_g2m_high ~ "G2/M high only",
  TRUE ~ "Below upper-quartile reference"
)

metadata <- combined[[]] |>
  tibble::rownames_to_column("cell") |>
  mutate(fine_cluster = as.character(res_0.2))

score_summary <- metadata |>
  group_by(fine_cluster, replicate) |>
  summarise(
    n_cells = n(),
    mean_s_score = mean(cell_cycle_s_score),
    median_s_score = median(cell_cycle_s_score),
    s_high_percent = 100 * mean(cell_cycle_s_high),
    mean_g2m_score = mean(cell_cycle_g2m_score),
    median_g2m_score = median(cell_cycle_g2m_score),
    g2m_high_percent = 100 * mean(cell_cycle_g2m_high),
    .groups = "drop"
  )

selected_cycle_genes <- tibble::tribble(
  ~phase, ~gene, ~display_name,
  "G2/M", "LDECv5g09615", "Topoisomerase II",
  "G2/M", "LDECv5g03544", "Cdc25/string-like",
  "G2/M", "LDECv5g10369", "UBE2C",
  "G2/M", "LDECv5g04564", "Anillin",
  "G2/M", "LDECv5g04603", "CDK1",
  "G2/M", "LDECv5g13420", "BUB1",
  "S", "LDECv5g10569", "PCNA",
  "S", "LDECv5g12228", "DNA polymerase alpha",
  "S", "LDECv5g07247", "FEN1",
  "S", "LDECv5g13349", "RNR-M2",
  "S", "LDECv5g09040", "MCM4",
  "S", "LDECv5g02665", "CDC6"
) |>
  filter(gene %in% rownames(combined)) |>
  mutate(label = paste0(phase, ": ", display_name)) |>
  left_join(
    mapping |>
      filter(mapping_status == "mapped_and_detected") |>
      distinct(phase, gene, dm_orth, old_gene, crosswalk_annotation),
    by = c("phase", "gene")
  )

write_tsv(score_summary, file.path(CELL_CYCLE_DIR, "scores_by_cluster.tsv"))
write_tsv(
  selected_cycle_genes,
  file.path(CELL_CYCLE_DIR, "selected_cell_cycle_genes.tsv")
)
write_tsv(
  data.frame(
    score = c("S", "G2/M"),
    upper_quartile_cutoff = c(s_cutoff, g2m_cutoff)
  ),
  file.path(CELL_CYCLE_DIR, "score_thresholds.tsv")
)

saveRDS(combined, file.path(OUT, "combined_cell_cycle.rds"))
capture.output(sessionInfo(), file = file.path(CELL_CYCLE_DIR, "sessionInfo.txt"))
