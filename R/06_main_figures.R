# Generate manuscript Figures 1-4.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stringr)
})

source(file.path("R", "config.R"))

CONSENSUS_OUT <- file.path(OUT, "replicate_consensus_GO")
ANNOTATION_OUT <- file.path(OUT, "final_annotation")
VALIDATION_OUT <- file.path(OUT, "cluster_validation")
FIG_OUT <- file.path(OUT, "figures")
dir.create(FIG_OUT, recursive = TRUE, showWarnings = FALSE)

combined <- JoinLayers(readRDS(file.path(
  OUT, "combined_annotated.rds"
)))
DefaultAssay(combined) <- "RNA"

annotation <- read.delim(
  file.path(ANNOTATION_OUT, "cluster_annotations.tsv"),
  check.names = FALSE
) |>
  mutate(fine_cluster = as.character(fine_cluster))

cluster_levels <- annotation$fine_cluster
short_levels <- annotation$short_label
cluster_short <- setNames(annotation$short_label, annotation$fine_cluster)
cluster_full <- setNames(annotation$cluster_annotation, annotation$fine_cluster)
cluster_axis <- setNames(
  paste0(annotation$fine_cluster, "  ", annotation$short_label),
  annotation$fine_cluster
)

annotation_palette <- c(
  "Proliferating plasmatocytes" = "#0072B2",
  "Granulocytes" = "#009E73",
  "Stress-responsive oenocytoids" = "#D55E00",
  "Adhesive/clotting plasmatocytes" = "#CC79A7",
  "Prohemocytes (provisional)" = "#E69F00",
  "Differentiating oenocytoids" = "#56B4E9",
  "Activated oenocytoids" = "#B79F00",
  "Germline contaminant" = "#7F7F7F"
)
cluster_palette <- setNames(annotation_palette[short_levels], cluster_levels)

md <- combined[[]][Cells(combined), , drop = FALSE]
md$cell <- rownames(md)
md$fine_cluster <- factor(as.character(md$res_0.2), levels = cluster_levels)
md$short_label <- factor(md$short_label, levels = short_levels)
md$cluster_display <- factor(
  cluster_axis[as.character(md$fine_cluster)],
  levels = unname(cluster_axis[cluster_levels])
)

umap <- as.data.frame(Embeddings(combined, reduction = "umap.harmony"))
colnames(umap)[1:2] <- c("UMAP_1", "UMAP_2")
umap$cell <- rownames(umap)
umap <- left_join(umap, md, by = "cell")
umap_centres <- umap |>
  group_by(fine_cluster) |>
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

theme_manuscript <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(colour = "grey30"),
      legend.key.height = unit(0.42, "cm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

annotated_umap <- function(show_legend = TRUE, title = NULL) {
  ggplot(umap, aes(UMAP_1, UMAP_2, colour = fine_cluster)) +
    geom_point(size = 0.16, alpha = 0.78, stroke = 0) +
    geom_label(
      data = umap_centres, aes(label = fine_cluster),
      colour = "black", fill = alpha("white", 0.85),
      size = 3, label.size = 0.15, label.padding = unit(0.12, "lines")
    ) +
    scale_colour_manual(
      values = cluster_palette,
      breaks = cluster_levels,
      labels = unname(cluster_axis[cluster_levels]),
      drop = FALSE
    ) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2", colour = NULL) +
    theme_manuscript(9) +
    theme(legend.position = if (show_legend) "right" else "none")
}

# Figure 1: integrated structure and replicate representation -----------------
p1_umap <- annotated_umap(
  title = "Established hemocyte classes and transcriptional states"
)

replicate_palette <- c(rep1 = "#0072B2", rep2 = "#D55E00")
p1_replicate <- ggplot(umap, aes(UMAP_1, UMAP_2, colour = replicate)) +
  geom_point(size = 0.16, alpha = 0.65, stroke = 0) +
  scale_colour_manual(values = replicate_palette) +
  labs(title = "Cells from both biological replicates", x = "UMAP 1",
       y = "UMAP 2", colour = "Sample") +
  theme_manuscript(9) +
  theme(legend.position = "top")

replicate_composition <- md |>
  count(replicate, fine_cluster, .drop = FALSE) |>
  group_by(replicate) |>
  mutate(percent_of_sample = 100 * n / sum(n)) |>
  ungroup()

p1_composition <- ggplot(
  replicate_composition,
  aes(fine_cluster, percent_of_sample, colour = replicate, group = replicate)
) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = replicate_palette) +
  scale_x_discrete(labels = cluster_levels) +
  labs(title = "Cluster abundance by replicate", x = "Fine cluster",
       y = "Cells in sample (%)", colour = "Sample") +
  theme_manuscript(9) +
  theme(legend.position = "top")

figure_1 <- p1_umap | (p1_replicate / p1_composition)
figure_1 <- figure_1 +
  plot_layout(widths = c(1.45, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(FIG_OUT, "figure_1.pdf"), figure_1,
       width = 14.5, height = 7.0, limitsize = FALSE)

# Figure 2: independent recovery, consensus markers and GO evidence -----------
recovery <- read.delim(
  file.path(CONSENSUS_OUT, "best_independent_cluster_recovery.tsv"),
  check.names = FALSE
) |>
  mutate(fine_cluster = as.character(fine_cluster))

recovery_long <- recovery |>
  select(fine_cluster, starts_with("precision_"), starts_with("recall_")) |>
  pivot_longer(-fine_cluster, names_to = "measure", values_to = "value") |>
  separate(measure, into = c("metric", "replicate"), sep = "_", extra = "merge") |>
  mutate(
    metric = recode(metric, precision = "Purity", recall = "Coverage"),
    replicate = recode(replicate, rep1 = "rep1", rep2 = "rep2"),
    fine_cluster = factor(fine_cluster, levels = cluster_levels)
  )

p2_recovery <- ggplot(
  recovery_long,
  aes(fine_cluster, value, colour = metric, shape = replicate,
      group = interaction(metric, replicate))
) +
  geom_hline(yintercept = 0.60, linetype = 2, colour = "grey55") +
  geom_line(alpha = 0.45, position = position_dodge(width = 0.26)) +
  geom_point(size = 2.3, position = position_dodge(width = 0.26)) +
  scale_colour_manual(values = c(Purity = "#0072B2", Coverage = "#D55E00")) +
  scale_y_continuous(limits = c(0, 1), labels = label_percent()) +
  labs(title = "Recovery in independently clustered replicates",
       subtitle = "Dashed line: pre-specified purity and coverage threshold",
       x = "Integrated fine cluster", y = "Best independent recovery",
       colour = NULL, shape = NULL) +
  theme_manuscript(9) +
  theme(legend.position = "top")

p2_markers <- annotation |>
  mutate(
    n_robust_markers = replace_na(n_robust_markers, 0),
    label_x = if_else(n_robust_markers == 0, 2, n_robust_markers),
    fine_cluster = factor(fine_cluster, levels = rev(cluster_levels)),
    recovered = if_else(independently_recovered, "Recovered", "Not recovered")
  ) |>
  ggplot(aes(n_robust_markers, fine_cluster, fill = recovered)) +
  geom_col(width = 0.7) +
  geom_text(aes(x = label_x, label = n_robust_markers), hjust = 0, size = 2.8) +
  scale_y_discrete(labels = function(x) cluster_axis[as.character(x)]) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
  scale_fill_manual(values = c(Recovered = "#2B8CBE", `Not recovered` = "#BDBDBD")) +
  labs(title = "Markers available for enrichment analysis",
       x = "Replicate-consensus positive markers", y = NULL, fill = NULL) +
  theme_manuscript(9) +
  theme(legend.position = "top", axis.text.y = element_text(size = 7.2))

go_top <- read.delim(
  file.path(CONSENSUS_OUT, "GO_enrichment_top_nonredundant.tsv"),
  check.names = FALSE
) |>
  mutate(consensus_cluster = as.character(consensus_cluster)) |>
  group_by(consensus_cluster) |>
  slice_min(p.adjust, n = 4, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    term_label = str_trunc(paste0(Description, " [", ontology, "]"), 52),
    consensus_cluster = factor(consensus_cluster, levels = cluster_levels)
  )

p2_go <- ggplot(
  go_top,
  aes(consensus_cluster, reorder(term_label, -log10(p.adjust)),
      size = Count, colour = -log10(p.adjust))
) +
  geom_point(alpha = 0.9) +
  scale_colour_viridis_c(option = "C", end = 0.9) +
  labs(title = "GO enrichment of replicate-consensus markers",
       subtitle = "No GO result is shown for cluster 5 (zero qualifying markers)",
       x = "Fine cluster", y = NULL, size = "Marker genes",
       colour = expression(-log[10](adjusted~P))) +
  theme_manuscript(8.5) +
  theme(axis.text.y = element_text(size = 7.2))

p2_umap <- annotated_umap(
  show_legend = TRUE,
  title = "Standardized CPB hemocyte annotation"
) + theme(legend.text = element_text(size = 6.8))

classification_markers <- read.delim(
  file.path(ANNOTATION_OUT, "selected_classification_markers.tsv"),
  check.names = FALSE
) |>
  filter(gene %in% rownames(combined)) |>
  mutate(label = paste0(criterion, ": ", display_name))
classification_marker_labels <- setNames(
  classification_markers$label, classification_markers$gene
)
combined$cluster_display <- md[Cells(combined), "cluster_display"]

p2_class_markers <- DotPlot(
  combined, features = classification_markers$gene,
  group.by = "cluster_display", assay = "RNA", dot.scale = 6,
  cols = c("#3B4CC0", "#F6C141")
) +
  coord_flip() +
  scale_x_discrete(labels = classification_marker_labels) +
  labs(
    title = "Orthologous and functional markers used for classical correspondence",
    subtitle = "Prohemocyte candidates are displayed as compatibility evidence, not specific markers",
    x = NULL, y = NULL, size = "Detected (%)"
  ) +
  guides(colour = guide_colourbar(title = "Scaled average\nexpression")) +
  theme_manuscript(8.5) +
  theme(axis.text.y = element_text(size = 6.8),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 6.8))

figure_2 <- (p2_recovery | p2_markers) / (p2_umap | p2_go) /
  p2_class_markers +
  plot_layout(heights = c(0.78, 1.05, 1.12)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(FIG_OUT, "figure_2.pdf"), figure_2,
       width = 16, height = 16.5, limitsize = FALSE)

# Figure 3: cell-cycle activity at fine-cluster resolution --------------------
hemocyte_md <- md |>
  filter(as.character(fine_cluster) != "8")

s_cutoff <- quantile(hemocyte_md$cell_cycle_s_score, 0.75, na.rm = TRUE)
g2m_cutoff <- quantile(hemocyte_md$cell_cycle_g2m_score, 0.75, na.rm = TRUE)

score_long <- hemocyte_md |>
  select(cell, fine_cluster, cluster_display, cell_cycle_s_score,
         cell_cycle_g2m_score) |>
  pivot_longer(starts_with("cell_cycle_"), names_to = "phase",
               values_to = "score") |>
  mutate(phase = recode(
    phase,
    cell_cycle_s_score = "S",
    cell_cycle_g2m_score = "G2/M"
  ))

cycle_summary <- score_long |>
  group_by(fine_cluster, cluster_display, phase) |>
  summarise(
    n_cells = n(),
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    pct_above_global_q75 = 100 * mean(score > if_else(
      first(phase) == "S", s_cutoff, g2m_cutoff), na.rm = TRUE),
    .groups = "drop"
  )

cutoff_data <- tibble(
  phase = c("S", "G2/M"), cutoff = c(s_cutoff, g2m_cutoff)
)

p3_scores <- ggplot(
  score_long,
  aes(fine_cluster, score, fill = fine_cluster)
) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.18) +
  geom_boxplot(width = 0.10, outlier.shape = NA, fill = "white", linewidth = 0.22) +
  geom_hline(data = cutoff_data, aes(yintercept = cutoff),
             colour = "grey35", linetype = 2, linewidth = 0.35) +
  facet_wrap(~phase, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = cluster_palette, drop = FALSE) +
  labs(title = "Cell-cycle gene-set scores in hemocyte clusters",
       subtitle = "Cluster 8 contaminant excluded; dashed line is the hemocyte-wide Q75",
       x = "Fine cluster", y = "Control-adjusted module score") +
  theme_manuscript(8.5) +
  theme(legend.position = "none")

p3_prevalence <- ggplot(
  cycle_summary,
  aes(fine_cluster, pct_above_global_q75, colour = phase, group = phase)
) +
  geom_hline(yintercept = 25, colour = "grey65", linetype = 2) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c(S = "#009E73", `G2/M` = "#C14953")) +
  labs(title = "Cells above the hemocyte-wide upper quartile",
       x = "Fine cluster", y = "Cells (%)", colour = "Gene set") +
  theme_manuscript(8.5) +
  theme(legend.position = "top")

cycle_features <- read.delim(
  file.path(OUT, "cell_cycle", "selected_cell_cycle_genes.tsv"),
  check.names = FALSE
) |>
  distinct(gene, phase, display_name) |>
  filter(gene %in% rownames(combined)) |>
  mutate(label = paste0(phase, ": ", display_name))
cycle_labels <- setNames(cycle_features$label, cycle_features$gene)

hemocyte_object <- subset(combined, cells = hemocyte_md$cell)
hemocyte_object$cluster_display <- hemocyte_md[
  Cells(hemocyte_object), "cluster_display"
]

p3_genes <- DotPlot(
  hemocyte_object, features = cycle_features$gene,
  group.by = "cluster_display", assay = "RNA", dot.scale = 6,
  cols = c("#3B4CC0", "#F6C141")
) +
  coord_flip() +
  scale_x_discrete(labels = cycle_labels) +
  labs(title = "Transferred cell-cycle orthologs", x = NULL, y = NULL,
       size = "Detected (%)") +
  guides(colour = guide_colourbar(title = "Scaled average\nexpression")) +
  theme_manuscript(8.5) +
  theme(axis.text.y = element_text(size = 7.1),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 7))

p3_s_umap <- FeaturePlot(
  hemocyte_object, features = "cell_cycle_s_score",
  reduction = "umap.harmony",
  min.cutoff = "q05", max.cutoff = "q95", pt.size = 0.18,
  raster = TRUE, cols = c("#F0F0F0", "#7A1F5C")
) +
  labs(title = "S gene-set score", x = "UMAP 1", y = "UMAP 2",
       colour = "Score") +
  theme_manuscript(8.5)

p3_g2m_umap <- FeaturePlot(
  hemocyte_object, features = "cell_cycle_g2m_score",
  reduction = "umap.harmony",
  min.cutoff = "q05", max.cutoff = "q95", pt.size = 0.18,
  raster = TRUE, cols = c("#F0F0F0", "#9D1D35")
) +
  labs(title = "G2/M gene-set score", x = "UMAP 1", y = "UMAP 2",
       colour = "Score") +
  theme_manuscript(8.5)

figure_3 <- ((p3_scores | p3_prevalence) / p3_genes /
                    (p3_s_umap | p3_g2m_umap)) +
  plot_layout(heights = c(1.0, 1.25, 1.0)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(FIG_OUT, "figure_3.pdf"), figure_3,
       width = 14.2, height = 14.8, limitsize = FALSE)

# Figure 4: clustering criteria, low-RNA candidate and contamination ----------
transition <- read.delim(
  file.path(VALIDATION_OUT,
            "transition_resolution_01_to_02.tsv"),
  check.names = FALSE
) |>
  filter(percent_of_broad_cluster >= 0.5) |>
  mutate(
    fine_cluster = factor(fine_cluster, levels = cluster_levels),
    broad_cluster = factor(broad_cluster, levels = 5:1)
  )

broad_axis <- c(
  `1` = "1  Proliferating plasmatocyte / granulocyte",
  `2` = "2  Oenocytoid branch",
  `3` = "3  Adhesive/clotting plasmatocytes",
  `4` = "4  Prohemocyte + contaminant",
  `5` = "5  Activated oenocytoids"
)

p4_hierarchy <- ggplot(
  transition,
  aes(fine_cluster, broad_cluster, fill = percent_of_broad_cluster)
) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = if_else(
    percent_of_broad_cluster >= 1,
    sprintf("%.0f%%", percent_of_broad_cluster), ""
  )), size = 2.7) +
  scale_y_discrete(labels = broad_axis) +
  scale_fill_viridis_c(option = "C", direction = -1) +
  labs(title = "Broad-to-fine clustering hierarchy",
       subtitle = "Percentages are within each resolution-0.1 cluster",
       x = "Resolution-0.2 fine cluster", y = "Resolution-0.1 cluster",
       fill = "% of broad\ncluster") +
  theme_manuscript(8.5)

support <- annotation |>
  left_join(
    recovery |>
      transmute(
        fine_cluster,
        minimum_recovery_f1 = pmin(f1_rep1, f1_rep2)
      ),
    by = "fine_cluster"
  ) |>
  mutate(
    n_robust_markers = replace_na(n_robust_markers, 0),
    recovered = if_else(independently_recovered, "Recovered", "Not recovered")
  )

p4_support <- ggplot(
  support,
  aes(minimum_recovery_f1, n_robust_markers, colour = fine_cluster,
      size = total_cells)
) +
  geom_vline(xintercept = 0.60, linetype = 2, colour = "grey60") +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey60") +
  geom_point(alpha = 0.9) +
  geom_text(aes(label = fine_cluster), colour = "black", size = 2.8,
            nudge_y = 6) +
  scale_colour_manual(values = cluster_palette, guide = "none") +
  scale_size_continuous(range = c(2.5, 8), labels = label_number()) +
  scale_x_continuous(limits = c(0, 1), labels = label_percent()) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.14))) +
  labs(title = "Independent recovery and positive-marker support",
       subtitle = "Cluster 5 is recovered but has no consensus positive markers",
       x = "Minimum recovery F1 across replicates",
       y = "Replicate-consensus markers", size = "Cells") +
  theme_manuscript(8.5)

complexity_long <- md |>
  transmute(fine_cluster, nFeature_RNA, nCount_RNA) |>
  pivot_longer(c(nFeature_RNA, nCount_RNA), names_to = "measure",
               values_to = "value") |>
  mutate(measure = recode(
    measure,
    nFeature_RNA = "Detected genes",
    nCount_RNA = "RNA counts"
  ))

p4_complexity <- ggplot(
  complexity_long,
  aes(fine_cluster, value, fill = fine_cluster)
) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.18) +
  geom_boxplot(width = 0.10, outlier.shape = NA, fill = "white", linewidth = 0.2) +
  facet_wrap(~measure, scales = "free_y", nrow = 1) +
  scale_y_log10(labels = label_number()) +
  scale_fill_manual(values = cluster_palette, guide = "none") +
  labs(title = "RNA complexity differs strongly among fine clusters",
       subtitle = "Low RNA supports neither viability nor lineage identity by itself",
       x = "Fine cluster", y = "Per-cell value (log10 scale)") +
  theme_manuscript(8.5)

p4_abundance <- ggplot(
  replicate_composition,
  aes(fine_cluster, percent_of_sample, fill = replicate)
) +
  geom_col(position = position_dodge(width = 0.78), width = 0.68) +
  scale_fill_manual(values = replicate_palette) +
  labs(title = "Sample-normalized cluster abundance",
       subtitle = "Cluster 8 is strongly replicate-skewed",
       x = "Fine cluster", y = "Cells in sample (%)", fill = "Sample") +
  theme_manuscript(8.5) +
  theme(legend.position = "top")

contamination <- read.delim(
  file.path(VALIDATION_OUT,
            "contamination_program_summary.tsv"),
  check.names = FALSE
) |>
  group_by(fine_cluster, program) |>
  summarise(
    flagged_percent = 100 * sum(flagged_cells) / sum(n_cells),
    .groups = "drop"
  ) |>
  complete(fine_cluster = as.integer(cluster_levels), program,
           fill = list(flagged_percent = 0)) |>
  mutate(fine_cluster = factor(fine_cluster, levels = cluster_levels))

p4_contamination <- ggplot(
  contamination,
  aes(fine_cluster, program, fill = flagged_percent)
) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = if_else(
    flagged_percent >= 0.05, sprintf("%.1f%%", flagged_percent), "0"
  )), size = 2.5) +
  scale_fill_gradient(low = "#F3F3F3", high = "#C14953") +
  labs(title = "Targeted collection-contamination screen",
       subtitle = "Program flags are screening evidence, not barcode-level proof",
       x = "Fine cluster", y = NULL, fill = "Cells\nflagged") +
  theme_manuscript(8.5)

sensitivity <- read.delim(
  file.path(VALIDATION_OUT,
            "graph_parameter_sensitivity.tsv"),
  check.names = FALSE
)

p4_sensitivity <- ggplot(
  sensitivity,
  aes(factor(k), factor(resolution), fill = cluster4_cohesion)
) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.2f", cluster4_cohesion)), size = 2.7) +
  scale_fill_gradient(limits = c(0, 1), low = "#F3F3F3", high = "#2A6F97") +
  labs(title = "Low-RNA branch cohesion across graph parameters",
       subtitle = "Sensitivity analysis of the resolution-0.1 candidate branch",
       x = "k nearest neighbours", y = "Resolution", fill = "Branch\ncohesion") +
  theme_manuscript(8.5)

prohemocyte_evidence <- read.delim(
  file.path(VALIDATION_OUT,
            "prohemocyte_evidence.tsv"),
  check.names = FALSE
) |>
  mutate(evidence = factor(evidence, levels = rev(evidence)))

p4_evidence <- ggplot(
  prohemocyte_evidence,
  aes(score, evidence, colour = direction)
) +
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
       subtitle = paste(
         "Morphology and immunostaining support CPB prohemocytes;",
         "the RNA-low cluster remains provisional"
       ),
       x = NULL, y = NULL, colour = NULL) +
  theme_manuscript(8.5) +
  theme(legend.position = "none")

figure_4 <- ((p4_hierarchy | p4_support) /
                    (p4_complexity | p4_abundance) /
                    (p4_contamination | p4_sensitivity) /
                    p4_evidence) +
  plot_layout(heights = c(1.05, 1.0, 0.95, 1.05)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(FIG_OUT, "figure_4.pdf"), figure_4,
       width = 15.5, height = 18.0, limitsize = FALSE)

# Export the data summarized in the figures.
write.table(
  replicate_composition,
  file.path(FIG_OUT, "figure1_replicate_composition.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  cycle_summary,
  file.path(FIG_OUT, "figure3_cell_cycle_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  support,
  file.path(FIG_OUT, "figure4_cluster_support_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  contamination,
  file.path(FIG_OUT, "figure4_contamination_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

capture.output(sessionInfo(), file = file.path(FIG_OUT, "sessionInfo_figures.txt"))
message("Figures 1-4 written to: ", FIG_OUT)
