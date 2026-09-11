# Multi-resolution cluster validation and contamination assessment.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(writexl)
})

source(file.path("R", "config.R"))

VALIDATION_OUT <- file.path(OUT, "cluster_validation")
dir.create(VALIDATION_OUT, recursive = TRUE, showWarnings = FALSE)

OBJECT_FILE <- file.path(OUT, "combined_cell_cycle.rds")
ANNOTATION_FILE <- input_file("LdecV5_functional_annotation.txt")
if (!file.exists(OBJECT_FILE)) stop("Missing input object: ", OBJECT_FILE)

SEED <- 20260824
MIN_REP_CELLS <- 15L
ROBUST_LOG2_DIFF <- 0.5
ROBUST_PCT_FOCAL <- 0.15
ROBUST_PCT_DIFF <- 0.10

set.seed(SEED)
combined <- JoinLayers(readRDS(OBJECT_FILE))
DefaultAssay(combined) <- "RNA"
md <- combined[[]][Cells(combined), , drop = FALSE]
md$cell <- rownames(md)
md$broad_cluster <- as.character(md$res_0.1)
md$fine_cluster <- as.character(md$res_0.2)
md$res_0.3 <- as.character(md$res_0.3)

rna_data <- GetAssayData(combined, assay = "RNA", layer = "data")
rna_counts <- GetAssayData(combined, assay = "RNA", layer = "counts")
annotation <- read.delim(ANNOTATION_FILE, check.names = FALSE, quote = "")
annotation_text <- paste(annotation$Uniprot_annotation,
                         annotation$Interproscan_annotation)

# -----------------------------------------------------------------------------
# 1. Multi-resolution hierarchy and branch-level evidence
# -----------------------------------------------------------------------------
resolution_columns <- c("res_0.1", "res_0.2", "res_0.3")

cluster_summary <- bind_rows(lapply(resolution_columns, function(column) {
  labels <- as.character(md[[column]])
  tibble(
    resolution = sub("res_", "", column),
    cluster = labels,
    replicate = md$replicate,
    nCount_RNA = md$nCount_RNA,
    nFeature_RNA = md$nFeature_RNA,
    percent_mt = md$percent.mt,
    s_high = md$cell_cycle_s_high,
    g2m_high = md$cell_cycle_g2m_high
  ) |>
    group_by(resolution, cluster) |>
    summarise(
      n_cells = n(),
      percent_cells = 100 * n() / nrow(md),
      rep1_cells = sum(replicate == "rep1"),
      rep2_cells = sum(replicate == "rep2"),
      median_UMIs = median(nCount_RNA),
      median_features = median(nFeature_RNA),
      median_percent_mt = median(percent_mt),
      s_high_percent = 100 * mean(s_high),
      g2m_high_percent = 100 * mean(g2m_high),
      .groups = "drop"
    )
}))

overall_minor_fraction <- min(table(md$replicate)) / nrow(md)
cluster_summary <- cluster_summary |>
  mutate(
    minor_replicate_fraction = pmin(rep1_cells, rep2_cells) / n_cells,
    replicate_balance = pmin(minor_replicate_fraction / overall_minor_fraction, 1)
  )

transition_01_02 <- as.data.frame(table(
  broad_cluster = md$broad_cluster,
  fine_cluster = md$fine_cluster
)) |>
  as_tibble() |>
  group_by(broad_cluster) |>
  mutate(percent_of_broad_cluster = 100 * Freq / sum(Freq)) |>
  ungroup()

transition_02_03 <- as.data.frame(table(
  fine_cluster = md$fine_cluster,
  cluster_03 = md$res_0.3
)) |>
  as_tibble() |>
  group_by(fine_cluster) |>
  mutate(percent_of_fine_cluster = 100 * Freq / sum(Freq)) |>
  ungroup()

replicate_effect <- function(labels, cluster, replicate_id) {
  rep_index <- which(md$replicate == replicate_id)
  focal <- rep_index[labels[rep_index] == cluster]
  reference <- rep_index[labels[rep_index] != cluster]
  if (length(focal) < MIN_REP_CELLS || length(reference) < MIN_REP_CELLS) {
    return(tibble(
      cluster = cluster, replicate = replicate_id, gene = rownames(rna_data),
      n_focal = length(focal), mean_log2_difference = NA_real_,
      pct_focal = NA_real_, pct_reference = NA_real_, pct_difference = NA_real_
    ))
  }
  focal_mean <- Matrix::rowMeans(rna_data[, focal, drop = FALSE]) / log(2)
  reference_mean <- Matrix::rowMeans(rna_data[, reference, drop = FALSE]) / log(2)
  focal_pct <- Matrix::rowMeans(rna_data[, focal, drop = FALSE] > 0)
  reference_pct <- Matrix::rowMeans(rna_data[, reference, drop = FALSE] > 0)
  tibble(
    cluster = cluster,
    replicate = replicate_id,
    gene = rownames(rna_data),
    n_focal = length(focal),
    mean_log2_difference = focal_mean - reference_mean,
    pct_focal = focal_pct,
    pct_reference = reference_pct,
    pct_difference = focal_pct - reference_pct
  )
}

fine_labels <- md$fine_cluster
fine_effects <- bind_rows(lapply(sort(unique(fine_labels)), function(cluster) {
  bind_rows(replicate_effect(fine_labels, cluster, "rep1"),
            replicate_effect(fine_labels, cluster, "rep2"))
}))

fine_robust_markers <- fine_effects |>
  select(cluster, replicate, gene, n_focal, mean_log2_difference,
         pct_focal, pct_reference, pct_difference) |>
  pivot_wider(
    names_from = replicate,
    values_from = c(n_focal, mean_log2_difference, pct_focal,
                    pct_reference, pct_difference)
  ) |>
  mutate(
    sufficient_replicates = n_focal_rep1 >= MIN_REP_CELLS &
      n_focal_rep2 >= MIN_REP_CELLS,
    min_log2_difference = pmin(mean_log2_difference_rep1,
                               mean_log2_difference_rep2),
    min_pct_focal = pmin(pct_focal_rep1, pct_focal_rep2),
    min_pct_difference = pmin(pct_difference_rep1, pct_difference_rep2),
    robust_marker = sufficient_replicates &
      min_log2_difference >= ROBUST_LOG2_DIFF &
      min_pct_focal >= ROBUST_PCT_FOCAL &
      min_pct_difference >= ROBUST_PCT_DIFF
  ) |>
  left_join(annotation, by = c("gene" = "Geneid"))

fine_marker_summary <- fine_robust_markers |>
  group_by(cluster) |>
  summarise(
    sufficient_replicates = all(sufficient_replicates),
    n_robust_markers = ifelse(first(sufficient_replicates),
                              sum(robust_marker), NA_integer_),
    replicate_effect_rho = ifelse(
      first(sufficient_replicates),
      cor(mean_log2_difference_rep1, mean_log2_difference_rep2,
          method = "spearman", use = "complete.obs"),
      NA_real_
    ),
    .groups = "drop"
  ) |>
  left_join(filter(cluster_summary, resolution == "0.2"), by = "cluster")

pooled_effect <- function(labels, cluster) {
  focal <- which(labels == cluster)
  reference <- which(labels != cluster)
  focal_mean <- Matrix::rowMeans(rna_data[, focal, drop = FALSE]) / log(2)
  reference_mean <- Matrix::rowMeans(rna_data[, reference, drop = FALSE]) / log(2)
  focal_pct <- Matrix::rowMeans(rna_data[, focal, drop = FALSE] > 0)
  reference_pct <- Matrix::rowMeans(rna_data[, reference, drop = FALSE] > 0)
  tibble(
    cluster = cluster,
    gene = rownames(rna_data),
    mean_log2_difference = focal_mean - reference_mean,
    pct_focal = focal_pct,
    pct_reference = reference_pct,
    pct_difference = focal_pct - reference_pct
  )
}

fine_pooled_markers <- bind_rows(lapply(sort(unique(fine_labels)), function(cluster) {
  pooled_effect(fine_labels, cluster)
})) |>
  left_join(annotation, by = c("gene" = "Geneid")) |>
  group_by(cluster) |>
  arrange(desc(mean_log2_difference), .by_group = TRUE) |>
  mutate(marker_rank = row_number()) |>
  ungroup()

# -----------------------------------------------------------------------------
# 2. Collection-contamination-associated expression programs
# -----------------------------------------------------------------------------
# These panels are screening tools, not validated CPB tissue markers. A cell is
# flagged only when several genes and a substantial UMI fraction agree.
contamination_programs <- list(
  `Contractile tissue` = c(
    "LDECv5g00007", "LDECv5g00064", "LDECv5g00822", "LDECv5g01701",
    "LDECv5g01916", "LDECv5g01917", "LDECv5g01918", "LDECv5g06942",
    "LDECv5g08196", "LDECv5g08197", "LDECv5g08383", "LDECv5g11124"
  ),
  `Cuticle/epidermis` = annotation$Geneid[
    grepl("cuticle protein|cuticular protein|endocuticle|exocuticle|chitin-binding type R&R",
          annotation_text, ignore.case = TRUE)
  ],
  `Fat-body storage` = c(
    "LDECv5g00071", "LDECv5g03374", "LDECv5g04043", "LDECv5g05088",
    "LDECv5g10652", "LDECv5g10653", "LDECv5g10654", "LDECv5g10655"
  ),
  `Germline/meiotic` = c(
    "LDECv5g02519", "LDECv5g12534", "LDECv5g10850", "LDECv5g08174",
    "LDECv5g02338"
  )
)
contamination_programs <- lapply(contamination_programs, intersect,
                                 y = rownames(rna_counts))

flag_rules <- tibble::tribble(
  ~program, ~minimum_genes, ~minimum_percent_umi,
  "Contractile tissue", 5L, 0.50,
  "Cuticle/epidermis", 3L, 0.50,
  "Fat-body storage", 2L, 0.25,
  "Germline/meiotic", 3L, 0.25
)

contamination_cell <- bind_rows(lapply(names(contamination_programs), function(program) {
  genes <- contamination_programs[[program]]
  rule <- flag_rules[flag_rules$program == program, , drop = FALSE]
  if (nrow(rule) != 1L) stop("Expected one contamination flag rule for ", program)
  detected <- Matrix::colSums(rna_counts[genes, , drop = FALSE] > 0)
  percent_umi <- 100 * Matrix::colSums(rna_counts[genes, , drop = FALSE]) /
    md$nCount_RNA
  tibble(
    cell = colnames(rna_counts),
    broad_cluster = md$broad_cluster,
    fine_cluster = md$fine_cluster,
    replicate = md$replicate,
    program = program,
    n_panel_genes = length(genes),
    n_genes_detected = detected,
    percent_umi = percent_umi,
    flagged = detected >= rule$minimum_genes[[1]] &
      percent_umi >= rule$minimum_percent_umi[[1]]
  )
}))

contamination_summary <- contamination_cell |>
  group_by(fine_cluster, replicate, program) |>
  summarise(
    n_cells = n(),
    median_percent_umi = median(percent_umi),
    q90_percent_umi = quantile(percent_umi, 0.90),
    mean_percent_umi = mean(percent_umi),
    flagged_cells = sum(flagged),
    flagged_percent = 100 * mean(flagged),
    .groups = "drop"
  )

contamination_any <- contamination_cell |>
  group_by(cell) |>
  summarise(
    contamination_candidate = any(flagged),
    contamination_programs = paste(program[flagged], collapse = "; "),
    .groups = "drop"
  ) |>
  left_join(select(md, cell, fine_cluster), by = "cell") |>
  mutate(
    contamination_candidate = contamination_candidate |
      fine_cluster == "8",
    contamination_programs = case_when(
      fine_cluster == "8" & contamination_programs == "" ~
        "Germline/meiotic",
      fine_cluster == "8" & !grepl("Germline/meiotic",
                                    contamination_programs) ~
        paste(contamination_programs, "Germline/meiotic", sep = "; "),
      TRUE ~ contamination_programs
    )
  ) |>
  select(-fine_cluster)

contamination_cluster_summary <- contamination_any |>
  left_join(select(md, cell, broad_cluster, fine_cluster, replicate), by = "cell") |>
  group_by(broad_cluster, fine_cluster, replicate) |>
  summarise(
    n_cells = n(),
    contamination_candidates = sum(contamination_candidate),
    contamination_candidate_percent = 100 * mean(contamination_candidate),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# 3. Graph-parameter sensitivity
# -----------------------------------------------------------------------------
adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  n <- sum(tab)
  total_pairs <- choose2(n)
  sum_cells <- sum(choose2(tab))
  sum_rows <- sum(choose2(rowSums(tab)))
  sum_cols <- sum(choose2(colSums(tab)))
  expected <- sum_rows * sum_cols / total_pairs
  maximum <- (sum_rows + sum_cols) / 2
  if (maximum == expected) return(1)
  (sum_cells - expected) / (maximum - expected)
}

algorithm <- if (requireNamespace("leidenbase", quietly = TRUE)) 4 else 1
harmony <- Embeddings(combined, "harmony")[, 1:20, drop = FALSE]
parameter_grid <- expand.grid(
  k = c(10L, 20L, 30L, 50L),
  resolution = c(0.05, 0.10, 0.15, 0.20),
  KEEP.OUT.ATTRS = FALSE
)

GRAPH_SENSITIVITY_FILE <- file.path(VALIDATION_OUT, "graph_parameter_sensitivity.tsv")
if (file.exists(GRAPH_SENSITIVITY_FILE)) {
  graph_sensitivity <- read.delim(GRAPH_SENSITIVITY_FILE,
                                  check.names = FALSE)
} else {
graph_sensitivity <- bind_rows(lapply(seq_len(nrow(parameter_grid)), function(i) {
  k <- parameter_grid$k[i]
  resolution <- parameter_grid$resolution[i]
  graph <- FindNeighbors(
    harmony, k.param = k, compute.SNN = TRUE,
    return.neighbor = FALSE, verbose = FALSE
  )
  fit <- FindClusters(
    graph$snn, resolution = resolution, algorithm = algorithm,
    random.seed = SEED, verbose = FALSE
  )
  labels <- as.character(fit[, 1])
  names(labels) <- rownames(fit)
  labels <- labels[md$cell]

  broad4 <- md$broad_cluster == "4"
  c4_table <- table(labels[broad4])
  dominant <- names(which.max(c4_table))
  tibble(
    k = k,
    resolution = resolution,
    n_clusters = length(unique(labels)),
    ari_vs_current_01 = adjusted_rand_index(labels, md$broad_cluster),
    ari_vs_current_02 = adjusted_rand_index(labels, md$fine_cluster),
    cluster4_cohesion = max(c4_table) / sum(c4_table),
    cluster4_group_purity = sum(broad4 & labels == dominant) /
      sum(labels == dominant)
  )
}))
}

# -----------------------------------------------------------------------------
# 4. Morphology-informed interpretation
# -----------------------------------------------------------------------------
broad_annotation <- tibble::tribble(
  ~cluster, ~broad_functional_label, ~broad_functional_short,
  ~classical_correspondence, ~confidence, ~interpretation,
  "1", "Immune-metabolic hemocytes", "Immune-metabolic",
  "Plasmatocyte/granulocyte correspondence unresolved", "moderate",
  "A broad immune-metabolic population containing cycling and non-cycling states; resolution 0.2 separates two reproducible branches but does not establish classical morphology.",
  "2", "PPO-high melanization-associated hemocytes", "PPO-high melanization",
  "Oenocytoid-like candidate", "moderate",
  "A reproducible PPO-high program compatible with melanization-associated hemocytes; PPO expression is not sufficient by itself to define morphology.",
  "3", "Adhesive/clotting hemocytes", "Adhesive/clotting",
  "Plasmatocyte-like candidate; granulocyte alternative remains", "moderate",
  "A reproducible adhesion, extracellular-matrix and clotting program with strong S-phase activity.",
  "4", "Putative prohemocyte-enriched population", "Putative prohemocytes",
  "Prohemocyte candidate supported by morphology and immunostaining", "low-to-moderate",
  "Resolution 0.2 separates a 30-cell germline/meiotic contaminant branch from the remaining low-RNA population. The residual approximately 6.1% frequency is compatible with independently observed small prohemocytes, but sparse RNA, absent robust markers and QC sensitivity prevent definitive barcode-level assignment.",
  "5", "Rare PPO2-high redox state", "Rare PPO2-high state",
  "Possible oenocytoid-related state; not an independent type", "low",
  "A reproducible rare PPO2/redox state closest to the PPO-high population; collection stress or another biological state remains possible."
)

fine_broad_annotation <- tibble::tribble(
  ~fine_cluster, ~working_fine_label, ~fine_interpretation,
  "1", "Cycling immune-metabolic", "Reproducible cycling branch of broad cluster 1; not equivalent to prohemocytes because it represents 35% of all cells and retains an immune-metabolic program.",
  "2", "Non-cycling immune-metabolic", "Reproducible non-cycling branch of broad cluster 1.",
  "3", "PPO-high major state", "Major PPO-high branch of broad cluster 2.",
  "4", "Adhesive/clotting", "Equivalent to the broad adhesive/clotting population.",
  "5", "Putative prohemocytes", "Low-RNA branch comprising approximately 6.1% of cells; morphology- and immunostaining-supported candidate but without a robust positive RNA marker program.",
  "6", "PPO-high Notch/FoxO state", "Smaller PPO-high branch enriched for Notch-, FoxO- and Tob2-like transcripts.",
  "7", "Rare PPO2-high state", "Equivalent to the rare PPO2/redox population.",
  "8", "Germline contaminant", "Thirty-cell, replicate-imbalanced branch expressing a coherent meiotic/germline program; exclude from hemocyte lineage interpretation."
)

prohemocyte_evidence <- tibble::tribble(
  ~evidence, ~direction, ~score, ~interpretation,
  "Morphology", "supports", 1,
  "Small morphologically defined prohemocytes are independently observed in CPB hemolymph.",
  "Immunostaining", "supports", 1,
  "Independent immunostaining supports the presence of prohemocytes; the marker and scoring criteria must be reported in Methods.",
  "Frequency concordance", "supports", 1,
  "After separating the 30-cell germline/meiotic branch, the residual low-RNA population comprises approximately 6.1% of cells, close to the independent morphology-based estimate, but this is not proof of identity.",
  "RNA complexity", "compatible", 0.5,
  "Low RNA content is compatible with small undifferentiated cells but also with damaged cells or incomplete droplets.",
  "Transferred progenitor genes", "inconclusive", 0,
  "Candidate Drosophila and mosquito markers are not specific to cluster 4 and orthology is uncertain.",
  "Replicate-positive marker program", "does not support", -1,
  "Cluster 4 has no robust positive one-versus-rest marker under the replicate-aware effect thresholds.",
  "QC robustness", "concern", -1,
  "Cluster 4 membership changes strongly when minimum detected-gene thresholds are raised.",
  "Contamination exclusion", "inconclusive", 0,
  "Filtered Cell Ranger matrices and annotation-derived panels cannot exclude damaged tissue cells or ambient RNA."
)

annotation_map <- setNames(broad_annotation$broad_functional_short,
                           broad_annotation$cluster)
combined$broad_functional_label <- unname(broad_annotation$broad_functional_label[
  match(as.character(combined$res_0.1), broad_annotation$cluster)
])
combined$broad_functional_short <- unname(
  annotation_map[as.character(combined$res_0.1)]
)
combined$validated_fine_cluster <- unname(as.character(combined$res_0.2))
combined$working_fine_label <- unname(fine_broad_annotation$working_fine_label[
  match(as.character(combined$res_0.2), fine_broad_annotation$fine_cluster)
])
combined$contamination_candidate <- unname(contamination_any$contamination_candidate[
  match(Cells(combined), contamination_any$cell)
])
combined$contamination_programs <- unname(contamination_any$contamination_programs[
  match(Cells(combined), contamination_any$cell)
])

# -----------------------------------------------------------------------------
# 5. Outputs and diagnostic figure
# -----------------------------------------------------------------------------
write.table(broad_annotation, file.path(VALIDATION_OUT, "broad_cluster_annotations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fine_broad_annotation, file.path(VALIDATION_OUT, "fine_broad_cluster_annotations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cluster_summary, file.path(VALIDATION_OUT, "multiresolution_cluster_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(transition_01_02, file.path(VALIDATION_OUT, "transition_resolution_01_to_02.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(transition_02_03, file.path(VALIDATION_OUT, "transition_resolution_02_to_03.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fine_marker_summary, file.path(VALIDATION_OUT, "fine_cluster_marker_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filter(fine_robust_markers, robust_marker),
            file.path(VALIDATION_OUT, "fine_cluster_robust_markers.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fine_pooled_markers,
            file.path(VALIDATION_OUT, "fine_cluster_pooled_effects.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(contamination_summary,
            file.path(VALIDATION_OUT, "contamination_program_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(contamination_cluster_summary,
            file.path(VALIDATION_OUT, "contamination_candidate_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(contamination_cell,
            file.path(VALIDATION_OUT, "contamination_programs_by_cell.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(graph_sensitivity, file.path(VALIDATION_OUT, "graph_parameter_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(prohemocyte_evidence,
            file.path(VALIDATION_OUT, "prohemocyte_evidence.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write_xlsx(
  list(
    annotations = broad_annotation,
    fine_annotations = fine_broad_annotation,
    cluster_summary = cluster_summary,
    transition_01_to_02 = transition_01_02,
    fine_marker_summary = fine_marker_summary,
    fine_robust_markers = filter(fine_robust_markers, robust_marker),
    contamination_summary = contamination_summary,
    contamination_candidates = contamination_cluster_summary,
    graph_sensitivity = graph_sensitivity,
    prohemocyte_evidence = prohemocyte_evidence
  ),
  file.path(VALIDATION_OUT, "cluster_validation.xlsx")
)

saveRDS(combined, file.path(OUT, "combined_validated.rds"))

validation_palette <- c(
  `Immune-metabolic` = "#2A6F97",
  `PPO-high melanization` = "#D17A00",
  `Adhesive/clotting` = "#2A9D6F",
  `Putative prohemocytes` = "#767676",
  `Rare PPO2-high state` = "#B64C69"
)

p_broad <- DimPlot(
  combined, reduction = "umap.harmony", group.by = "broad_functional_short",
  cols = validation_palette, label = TRUE, repel = TRUE, raster = TRUE
) +
  labs(title = "Broad reproducible programs (resolution 0.1)", colour = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "none")

p_fine <- DimPlot(
  combined, reduction = "umap.harmony", group.by = "working_fine_label",
  label = TRUE, repel = TRUE, raster = TRUE
) +
  labs(title = "Secondary state-level partition (resolution 0.2)", colour = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "none")

p_transition <- transition_01_02 |>
  filter(Freq > 0) |>
  ggplot(aes(factor(fine_cluster), factor(broad_cluster),
             fill = percent_of_broad_cluster)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(percent_of_broad_cluster >= 1,
                               sprintf("%.0f%%", percent_of_broad_cluster), "")),
            size = 3) +
  scale_fill_gradient(low = "#F1F1F1", high = "#2A6F97") +
  labs(title = "Resolution hierarchy", x = "Resolution 0.2 cluster",
       y = "Resolution 0.1 cluster", fill = "% of broad\ncluster") +
  theme_classic(base_size = 9)

p_branch <- fine_marker_summary |>
  mutate(marker_count_plot = replace_na(n_robust_markers, 0),
         replicate_test = ifelse(sufficient_replicates,
                                 "Both replicates", "Insufficient cells")) |>
  ggplot(aes(replicate_balance, marker_count_plot, label = cluster,
             shape = replicate_test)) +
  geom_hline(yintercept = 5, linetype = 2, colour = "grey60") +
  geom_vline(xintercept = 0.75, linetype = 2, colour = "grey60") +
  geom_point(aes(size = n_cells, colour = median_features), alpha = 0.9) +
  geom_text(nudge_y = 20, size = 3) +
  scale_colour_gradient(low = "#D55E00", high = "#0072B2") +
  labs(title = "Fine branches require independent support",
       x = "Replicate balance", y = "Robust markers",
       size = "Cells", colour = "Median\ngenes", shape = NULL) +
  theme_classic(base_size = 9)

p_contamination <- contamination_summary |>
  group_by(fine_cluster, program) |>
  summarise(flagged_percent = sum(flagged_cells) / sum(n_cells) * 100,
            .groups = "drop") |>
  ggplot(aes(factor(fine_cluster), program, fill = flagged_percent)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f%%", flagged_percent)), size = 2.8) +
  scale_fill_gradient(low = "#F3F3F3", high = "#C14953") +
  labs(title = "Collection-contamination screen", x = "Resolution 0.2 cluster",
       y = NULL, fill = "Cells\nflagged") +
  theme_classic(base_size = 9)

p_sensitivity <- graph_sensitivity |>
  ggplot(aes(factor(k), factor(resolution), fill = cluster4_cohesion)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", cluster4_cohesion)), size = 3) +
  scale_fill_gradient(limits = c(0, 1), low = "#F3F3F3", high = "#2A6F97") +
  labs(title = "Candidate-population cohesion across graph parameters",
       x = "k nearest neighbours", y = "Resolution", fill = "Cluster 4\ncohesion") +
  theme_classic(base_size = 9)

p_evidence <- prohemocyte_evidence |>
  mutate(evidence = factor(evidence, levels = rev(evidence))) |>
  ggplot(aes(score, evidence, colour = direction)) +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_segment(aes(x = 0, xend = score, yend = evidence), linewidth = 0.7) +
  geom_point(size = 3) +
  scale_x_continuous(limits = c(-1.1, 1.1), breaks = c(-1, 0, 1),
                     labels = c("Concern", "Inconclusive", "Supports")) +
  scale_colour_manual(values = c(
    supports = "#2A9D6F", compatible = "#2A6F97",
    inconclusive = "#767676", `does not support` = "#D17A00",
    concern = "#C14953"
  )) +
  labs(title = "Evidence for a prohemocyte-enriched population",
       x = NULL, y = NULL, colour = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "none")

validation_figure <- ((p_broad | p_fine) /
                  (p_transition | p_branch) /
                  (p_contamination | p_sensitivity) /
                  p_evidence) +
  plot_layout(heights = c(1.15, 1, 0.9, 1.15)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(VALIDATION_OUT, "cluster_validation_diagnostics.pdf"),
       validation_figure,
       width = 15.5, height = 18)

capture.output(sessionInfo(), file = file.path(VALIDATION_OUT, "sessionInfo.txt"))
message("Cluster validation outputs written to: ", VALIDATION_OUT)
