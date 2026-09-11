# GO enrichment and annotation-marker expression for CPB hemocyte clusters.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggtext)
  library(patchwork)
  library(stringr)
})

source(file.path("R", "config.R"))

CONSENSUS_OUT <- file.path(OUT, "replicate_consensus_GO")
ANNOTATION_OUT <- file.path(OUT, "final_annotation")
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
cluster_axis <- setNames(
  paste0(annotation$fine_cluster, "  ", annotation$short_label),
  annotation$fine_cluster
)
cell_cluster <- as.character(combined$res_0.2)
rna <- GetAssayData(combined, assay = "RNA", layer = "data")

numeric_id <- function(gene) {
  sub("^LDECv5g", "", gene)
}

format_gene_label <- function(gene_name, gene, drosophila_evidence_class,
                              prefix = NULL) {
  id <- numeric_id(gene)
  if (identical(drosophila_evidence_class, "FlyBase-mapped ortholog")) {
    id <- paste0(
      "<span style='color:#0072B2;font-weight:700'>", id, "</span>"
    )
  } else if (identical(
    drosophila_evidence_class, "Drosophila protein annotation only"
  )) {
    id <- paste0("<span style='font-weight:700'>", id, "</span>")
  }
  label <- paste0(gene_name, " (", id, ")")
  if (!is.null(prefix)) paste0(prefix, " | ", label) else label
}

manual_names <- c(
  LDECv5g13052 = "Uncharacterized protein",
  LDECv5g11412 = "FHA/Ki67-like protein",
  LDECv5g01134 = "La protein homolog",
  LDECv5g01117 = "Myc",
  LDECv5g13007 = "Nucleolar protein 58",
  LDECv5g13594 = "PPO paralog C",
  LDECv5g15334 = "Uncharacterized protein",
  LDECv5g07766 = "Peroxidasin",
  LDECv5g05220 = "Dual oxidase",
  LDECv5g09650 = "BRI3-like protein",
  LDECv5g01707 = "Galactose mutarotase",
  LDECv5g05970 = "Small heat-shock protein (L(2)efl)",
  LDECv5g00509 = "Notch",
  LDECv5g00320 = "FoxO",
  LDECv5g13137 = "Ficolin-like protein",
  LDECv5g03007 = "Hemocytin/Hml-like",
  LDECv5g10301 = "Hemocyte transglutaminase",
  LDECv5g10644 = "Papilin",
  LDECv5g15152 = "Slowpoke potassium channel",
  LDECv5g13227 = "Laminin alpha",
  LDECv5g01578 = "Tob2",
  LDECv5g09082 = "Mth2-like receptor",
  LDECv5g13247 = "Fibrinogen C-domain protein",
  LDECv5g14050 = "DM9-domain protein",
  LDECv5g11330 = "Aminoacylase-1",
  LDECv5g04564 = "Anillin",
  LDECv5g10569 = "PCNA",
  LDECv5g09615 = "Topoisomerase II",
  LDECv5g13695 = "Rac2-like",
  LDECv5g13501 = "Integrin beta-PS",
  LDECv5g09037 = "Integrin alpha-PS3",
  LDECv5g15036 = "PPO paralog A",
  LDECv5g15038 = "PPO paralog B",
  LDECv5g08619 = "Pebbled-like",
  LDECv5g08078 = "Nimrod/Eater-like",
  LDECv5g15611 = "Hemocyte transglutaminase",
  LDECv5g06984 = "SPARC-like",
  LDECv5g00806 = "Ance-like"
)

display_name <- function(gene, uniprot = NA_character_, interpro = NA_character_) {
  manual <- unname(manual_names[gene])
  if (!is.na(manual)) return(manual)
  if (!is.na(uniprot) && nzchar(uniprot)) {
    return(str_trim(str_split_fixed(uniprot, ";", 2)[1]))
  }
  if (!is.na(interpro) && nzchar(interpro)) {
    return(str_trim(str_split_fixed(interpro, "\\|", 2)[1]))
  }
  "Uncharacterized protein"
}

scale_by_gene <- function(values) {
  if (all(is.na(values)) || sd(values, na.rm = TRUE) == 0) {
    return(rep(0, length(values)))
  }
  pmax(-2, pmin(2, as.numeric(scale(values))))
}

theme_marker <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(colour = "grey30"),
      strip.background = element_rect(fill = "#F2F2F2", colour = NA),
      strip.text.y.left = element_text(face = "bold", angle = 0),
      panel.spacing.y = unit(0.10, "lines")
    )
}

# Collate Drosophila support while retaining the provenance of each match.
drosophila_marker_map <- read.delim(
  file.path(OUT, "annotation_transfer", "drosophila_marker_old_to_v5_mapping.tsv"),
  check.names = FALSE
) |>
  filter(
    mapping_status == "mapped_and_detected", !is.na(gene),
    !is.na(`Flybase ID`), nzchar(`Flybase ID`)
  ) |>
  transmute(
    gene,
    flybase_id = `Flybase ID`,
    drosophila_evidence = "Mapped Drosophila hemocyte-marker ortholog"
  )

drosophila_cell_cycle_map <- read.delim(
  file.path(OUT, "cell_cycle", "selected_cell_cycle_genes.tsv"),
  check.names = FALSE
) |>
  filter(!is.na(gene), !is.na(dm_orth), nzchar(dm_orth)) |>
  transmute(
    gene,
    flybase_id = dm_orth,
    drosophila_evidence = "Mapped Drosophila cell-cycle ortholog"
  )

all_consensus_markers <- read.delim(
  file.path(CONSENSUS_OUT, "replicate_consensus_markers.tsv"),
  check.names = FALSE
)

drosophila_protein_map <- all_consensus_markers |>
  filter(
    !is.na(Uniprot_annotation),
    str_detect(Uniprot_annotation, regex("Drosophila", ignore_case = TRUE))
  ) |>
  distinct(gene) |>
  transmute(
    gene,
    flybase_id = NA_character_,
    drosophila_evidence = "Drosophila protein annotation"
  )

drosophila_support <- bind_rows(
  drosophila_marker_map,
  drosophila_cell_cycle_map,
  drosophila_protein_map
) |>
  group_by(gene) |>
  summarise(
    drosophila_supported = TRUE,
    flybase_id = {
      ids <- sort(unique(flybase_id[!is.na(flybase_id) & nzchar(flybase_id)]))
      if (length(ids)) paste(ids, collapse = "; ") else NA_character_
    },
    drosophila_evidence = paste(
      sort(unique(drosophila_evidence)), collapse = "; "
    ),
    .groups = "drop"
  ) |>
  mutate(
    drosophila_evidence_class = if_else(
      !is.na(flybase_id) & nzchar(flybase_id),
      "FlyBase-mapped ortholog",
      "Drosophila protein annotation only"
    )
  )

# Panel A: five most significant Biological Process terms per cluster.
go_group_levels <- c(
  "Cluster 1 | Proliferating plasmatocytes",
  "Cluster 2 | Granulocytes",
  "Cluster 3 | Stress-responsive oenocytoids",
  "Cluster 4 | Adhesive/clotting plasmatocytes",
  "Cluster 6 | Differentiating oenocytoids",
  "Cluster 7 | Activated oenocytoids"
)

go_group_labels <- c(
  `1` = go_group_levels[1], `2` = go_group_levels[2],
  `3` = go_group_levels[3], `4` = go_group_levels[4],
  `6` = go_group_levels[5], `7` = go_group_levels[6]
)

go_plot_data <- read.delim(
  file.path(CONSENSUS_OUT, "GO_enrichment_all.tsv"), check.names = FALSE
) |>
  filter(
    ontology == "BP",
    as.character(consensus_cluster) %in% names(go_group_labels)
  ) |>
  mutate(consensus_cluster = as.character(consensus_cluster)) |>
  arrange(consensus_cluster, p.adjust, desc(fold_enrichment)) |>
  distinct(consensus_cluster, Description, .keep_all = TRUE) |>
  group_by(consensus_cluster) |>
  slice_min(p.adjust, n = 5, with_ties = FALSE) |>
  arrange(p.adjust, .by_group = TRUE) |>
  mutate(term_rank = row_number()) |>
  ungroup() |>
  mutate(
    go_group_label = factor(
      unname(go_group_labels[consensus_cluster]), levels = go_group_levels
    ),
    expression_cluster = factor(consensus_cluster, levels = cluster_levels),
    term_key = paste(consensus_cluster, term_rank, Description, sep = "__"),
    term_label = str_wrap(Description, width = 48),
    significance = -log10(p.adjust)
  )

go_term_order <- go_plot_data |>
  arrange(match(consensus_cluster, names(go_group_labels)), term_rank) |>
  pull(term_key)
go_plot_data$term_key <- factor(
  go_plot_data$term_key, levels = rev(go_term_order)
)

p_go <- ggplot(
  go_plot_data,
  aes(fold_enrichment, term_key, size = Count, colour = significance)
) +
  geom_point(alpha = 0.92) +
  facet_wrap(~go_group_label, ncol = 2, scales = "free_y", drop = FALSE) +
  scale_y_discrete(
    labels = setNames(go_plot_data$term_label, go_plot_data$term_key)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  scale_size_continuous(range = c(2.2, 7.2), breaks = c(4, 8, 12, 16)) +
  scale_colour_viridis_c(option = "C", end = 0.92) +
  labs(
    title = "Biological Process enrichment of replicate-consensus markers",
    subtitle = str_wrap(paste(
      "Five terms with the lowest adjusted P value per supported cluster;",
      "all have adjusted P < 0.05 except cluster 6 rank 5 (0.051).",
      "Cluster 5 has no qualifying markers and cluster 8 is excluded."
    ), width = 125),
    x = "Fold enrichment", y = NULL,
    size = "Marker genes", colour = expression(-log[10](adjusted~P))
  ) +
  theme_marker(8.5) +
  theme(
    strip.text = element_text(size = 7.4, face = "bold"),
    axis.text.y = element_text(size = 7.1, colour = "grey15", lineheight = 0.9),
    panel.grid.major.x = element_line(colour = "grey91", linewidth = 0.25),
    panel.spacing = unit(0.45, "lines"),
    legend.position = "right"
  )

# Panel B: genes used in the main text to support classical annotations.
support_markers <- read.delim(
  file.path(ANNOTATION_OUT, "selected_classification_markers.tsv"),
  check.names = FALSE
) |>
  left_join(drosophila_support, by = "gene") |>
  mutate(drosophila_supported = coalesce(drosophila_supported, FALSE)) |>
  mutate(
    support_cluster = case_when(
      criterion == "Proliferation" ~ "1",
      criterion == "Granulocyte" ~ "2",
      criterion == "Granulocyte/Oenocytoid" ~ "2",
      str_starts(criterion, "Oenocytoid") ~ "3/6/7",
      criterion == "Plasmatocyte" ~ "4",
      str_starts(criterion, "Prohemocyte") ~ "5",
      TRUE ~ "Other"
    ),
    support_order = match(support_cluster, c("1", "2", "3/6/7", "4", "5")),
    support_group_label = recode(
      support_cluster,
      `1` = "Cluster 1 | Proliferating plasmatocytes",
      `2` = "Cluster 2 | Granulocytes",
      `3/6/7` = "Clusters 3/6/7 | Oenocytoid branch",
      `4` = "Cluster 4 | Adhesive/clotting plasmatocytes",
      `5` = "Cluster 5 | Prohemocytes (provisional)"
    ),
    gene_name = vapply(gene, display_name, character(1)),
    gene_label = mapply(
      function(name, id, supported) {
        format_gene_label(name, id, supported)
      },
      gene_name, gene, drosophila_evidence_class,
      USE.NAMES = FALSE
    )
  ) |>
  arrange(support_order, row_number())

dotplot_data <- bind_rows(lapply(seq_len(nrow(support_markers)), function(i) {
  bind_rows(lapply(cluster_levels, function(cluster) {
    values <- rna[support_markers$gene[i], cell_cluster == cluster]
    tibble(
      support_cluster = support_markers$support_cluster[i],
      support_order = support_markers$support_order[i],
      support_group_label = support_markers$support_group_label[i],
      marker_order = i,
      gene = support_markers$gene[i],
      gene_name = support_markers$gene_name[i],
      gene_label = support_markers$gene_label[i],
      drosophila_supported = support_markers$drosophila_supported[i],
      flybase_id = support_markers$flybase_id[i],
      drosophila_evidence = support_markers$drosophila_evidence[i],
      drosophila_evidence_class =
        support_markers$drosophila_evidence_class[i],
      expression_cluster = cluster,
      average_log_normalized_expression = mean(values),
      percent_detected = 100 * mean(values > 0)
    )
  }))
})) |>
  group_by(gene) |>
  mutate(scaled_average_expression = scale_by_gene(
    average_log_normalized_expression
  )) |>
  ungroup() |>
  arrange(marker_order)

dot_label_order <- dotplot_data |>
  distinct(marker_order, gene_label) |>
  arrange(marker_order)
dotplot_data$gene_label <- factor(
  dotplot_data$gene_label, levels = rev(dot_label_order$gene_label)
)
dotplot_data$expression_cluster <- factor(
  dotplot_data$expression_cluster, levels = cluster_levels
)
dotplot_data$support_group_label <- factor(
  dotplot_data$support_group_label,
  levels = c(
    "Cluster 1 | Proliferating plasmatocytes",
    "Cluster 2 | Granulocytes",
    "Clusters 3/6/7 | Oenocytoid branch",
    "Cluster 4 | Adhesive/clotting plasmatocytes",
    "Cluster 5 | Prohemocytes (provisional)"
  )
)

p_dotplot <- ggplot(
  dotplot_data,
  aes(expression_cluster, gene_label,
      size = percent_detected, colour = scaled_average_expression)
) +
  geom_point(alpha = 0.95) +
  facet_grid(
    support_group_label ~ ., scales = "free_y", space = "free_y",
    switch = "y", drop = FALSE
  ) +
  scale_size_continuous(
    range = c(0.35, 6.2), limits = c(0, 100), breaks = c(25, 50, 75, 100)
  ) +
  scale_colour_gradient2(
    low = "#3B6FB6", mid = "#F2F2F2", high = "#C43C39",
    midpoint = 0, limits = c(-2, 2), oob = scales::squish
  ) +
  labs(
    title = "Markers supporting the standardized cluster annotations",
    subtitle = paste(
      "Genes are grouped by the cluster annotation they support;",
      "clusters 3/6/7 form the shared oenocytoid branch"
    ),
    x = "Expression cluster", y = NULL,
    size = "Detected (%)", colour = "Scaled average\nexpression"
  ) +
  theme_marker(8.5) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(size = 7, face = "bold", angle = 0),
    axis.text.y = ggtext::element_markdown(size = 7, colour = "grey15"),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
    panel.spacing.y = unit(0.16, "lines"),
    legend.position = "right"
  )

figure_5_v12 <- (p_go / p_dotplot) +
  plot_layout(heights = c(1, 1.12)) +
  plot_annotation(
    tag_levels = "A",
    caption = paste(
      "Blue bold IDs: explicit FlyBase-mapped Drosophila orthologs; bold black",
      "IDs: Drosophila protein-annotation support only."
    ),
    theme = theme(
      plot.caption = element_text(size = 8, colour = "#0072B2", hjust = 0)
    )
  )

ggsave(
  file.path(FIG_OUT, "figure_5.pdf"), figure_5_v12,
  width = 14.5, height = 16.5, limitsize = FALSE
)

write.table(
  go_plot_data |>
    mutate(term_label = str_replace_all(term_label, "\\n", " ")),
  file.path(FIG_OUT, "figure5_top5_BP_GO_terms.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  dotplot_data,
  file.path(FIG_OUT, "figure5_annotation_marker_dotplot_data.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  drosophila_support |>
    filter(gene %in% support_markers$gene) |>
    arrange(gene),
  file.path(FIG_OUT, "figure5_drosophila_support.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
capture.output(sessionInfo(), file = file.path(FIG_OUT, "sessionInfo_marker_figure.txt"))
message("Figure 5 written to: ", FIG_OUT)
