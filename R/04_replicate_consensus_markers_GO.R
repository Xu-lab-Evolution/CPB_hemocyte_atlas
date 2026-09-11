# Replicate-consensus clustering, marker analysis, and GO enrichment.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(GO.db)
  library(AnnotationDbi)
  library(writexl)
})

source(file.path("R", "config.R"))

CONSENSUS_OUT <- file.path(OUT, "replicate_consensus_GO")
FIG_OUT <- file.path(CONSENSUS_OUT, "figures")
dir.create(FIG_OUT, recursive = TRUE, showWarnings = FALSE)

OBJECT_FILE <- file.path(OUT, "combined_validated.rds")
FUNCTION_FILE <- input_file("LdecV5_functional_annotation.txt")
ORTHOLOG_FILE <- file.path(OUT, "annotation_transfer",
                           "drosophila_marker_expression_by_cluster.tsv")
required_files <- c(OBJECT_FILE, FUNCTION_FILE, ORTHOLOG_FILE)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

SEED <- 20260824
N_HVG <- 3000L
N_PCS <- 50L
DIMS <- 1:20
K_NN <- 20L
INDEPENDENT_RESOLUTION <- 0.2
RESOLUTION_GRID <- seq(0.1, 0.8, by = 0.1)
MIN_REP_CELLS <- 15L
MIN_MATCH_RHO <- 0.50
MIN_FINE_PURITY <- 0.50
MIN_RECOVERY_PRECISION <- 0.60
MIN_RECOVERY_RECALL <- 0.60
ROBUST_LOG2_DIFF <- 0.5
ROBUST_PCT_FOCAL <- 0.15
ROBUST_PCT_DIFF <- 0.10

set.seed(SEED)
combined <- JoinLayers(readRDS(OBJECT_FILE))
DefaultAssay(combined) <- "RNA"
counts <- GetAssayData(combined, assay = "RNA", layer = "counts")
metadata <- combined[[]][Cells(combined), , drop = FALSE]
metadata$fine_cluster <- as.character(metadata$res_0.2)

# -----------------------------------------------------------------------------
# 1. Cluster each replicate independently
# -----------------------------------------------------------------------------
cluster_replicate <- function(replicate_id) {
  cells <- rownames(metadata)[metadata$replicate == replicate_id]
  obj <- CreateSeuratObject(
    counts = counts[, cells, drop = FALSE],
    meta.data = metadata[cells, , drop = FALSE],
    project = paste0("CPB_", replicate_id)
  )
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = N_HVG, verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, npcs = N_PCS, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "pca", dims = DIMS,
                       k.param = K_NN, verbose = FALSE)
  algorithm <- if (requireNamespace("leidenbase", quietly = TRUE)) 4 else 1
  obj <- FindClusters(obj, resolution = INDEPENDENT_RESOLUTION,
                      algorithm = algorithm, random.seed = SEED,
                      verbose = FALSE)
  obj$independent_cluster <- as.character(Idents(obj))
  obj
}

REPLICATE_CACHE <- file.path(CONSENSUS_OUT, "independent_replicate_objects.rds")
if (file.exists(REPLICATE_CACHE)) {
  rep_objects <- readRDS(REPLICATE_CACHE)
} else {
  rep_objects <- list(
    rep1 = cluster_replicate("rep1"),
    rep2 = cluster_replicate("rep2")
  )
  saveRDS(rep_objects, REPLICATE_CACHE)
}

# Resolution is sample-size dependent. Test a grid and measure whether each
# integrated branch is recovered with adequate purity and coverage.
GRID_CACHE <- file.path(CONSENSUS_OUT, "independent_cluster_grid.rds")
if (file.exists(GRID_CACHE)) {
  cluster_grid <- readRDS(GRID_CACHE)
} else {
  algorithm <- if (requireNamespace("leidenbase", quietly = TRUE)) 4 else 1
  cluster_grid <- lapply(rep_objects, function(obj) {
    setNames(lapply(RESOLUTION_GRID, function(resolution) {
      if (isTRUE(all.equal(resolution, INDEPENDENT_RESOLUTION))) {
        return(as.character(obj$independent_cluster))
      }
      clustered <- FindClusters(
        obj, resolution = resolution, algorithm = algorithm,
        random.seed = SEED, verbose = FALSE
      )
      as.character(Idents(clustered))
    }), sprintf("%.1f", RESOLUTION_GRID))
  })
  saveRDS(cluster_grid, GRID_CACHE)
}

cluster_recovery <- bind_rows(lapply(names(rep_objects), function(rep_id) {
  truth <- as.character(rep_objects[[rep_id]]$fine_cluster)
  bind_rows(lapply(names(cluster_grid[[rep_id]]), function(res_key) {
    labels <- cluster_grid[[rep_id]][[res_key]]
    bind_rows(lapply(sort(unique(truth)), function(fine_id) {
      scores <- bind_rows(lapply(sort(unique(labels)), function(cluster_id) {
        overlap <- sum(truth == fine_id & labels == cluster_id)
        precision <- overlap / sum(labels == cluster_id)
        recall <- overlap / sum(truth == fine_id)
        tibble(
          independent_cluster = cluster_id,
          overlap_cells = overlap,
          precision = precision,
          recall = recall,
          f1 = ifelse(precision + recall > 0,
                      2 * precision * recall / (precision + recall), 0)
        )
      }))
      scores |>
        slice_max(f1, n = 1, with_ties = FALSE) |>
        mutate(
          replicate = rep_id,
          resolution = as.numeric(res_key),
          fine_cluster = fine_id,
          fine_cluster_cells = sum(truth == fine_id)
        )
    }))
  }))
}))

best_recovery <- cluster_recovery |>
  group_by(replicate, fine_cluster) |>
  arrange(desc(f1), resolution, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup()

recovery_pairs <- best_recovery |>
  dplyr::select(replicate, fine_cluster, resolution, independent_cluster,
                overlap_cells, precision, recall, f1, fine_cluster_cells) |>
  pivot_wider(
    names_from = replicate,
    values_from = c(resolution, independent_cluster, overlap_cells, precision,
                    recall, f1, fine_cluster_cells)
  ) |>
  mutate(
    passes_recovery = fine_cluster_cells_rep1 >= MIN_REP_CELLS &
      fine_cluster_cells_rep2 >= MIN_REP_CELLS &
      precision_rep1 >= MIN_RECOVERY_PRECISION &
      precision_rep2 >= MIN_RECOVERY_PRECISION &
      recall_rep1 >= MIN_RECOVERY_RECALL &
      recall_rep2 >= MIN_RECOVERY_RECALL
  )

independent_summary <- bind_rows(lapply(names(rep_objects), function(rep_id) {
  md <- rep_objects[[rep_id]][[]]
  as_tibble(md, rownames = "cell") |>
    count(replicate, independent_cluster, fine_cluster, name = "n_cells") |>
    group_by(replicate, independent_cluster) |>
    mutate(
      cluster_total = sum(n_cells),
      fine_fraction = n_cells / cluster_total
    ) |>
    ungroup()
}))

independent_dominant <- independent_summary |>
  group_by(replicate, independent_cluster) |>
  slice_max(n_cells, n = 1, with_ties = FALSE) |>
  ungroup() |>
  transmute(
    replicate,
    independent_cluster,
    dominant_fine_cluster = fine_cluster,
    fine_cluster_purity = fine_fraction,
    n_cells = cluster_total
  )

# Compare cluster-relative pseudobulk profiles. Subtracting the replicate-wide
# mean for each gene prevents housekeeping abundance from driving the match.
pseudobulk_relative <- function(obj, features) {
  data <- GetAssayData(obj, assay = "RNA", layer = "data")[features, , drop = FALSE]
  labels <- obj$independent_cluster
  bulk <- sapply(sort(unique(labels)), function(cluster_id) {
    Matrix::rowMeans(data[, labels == cluster_id, drop = FALSE])
  })
  if (is.null(dim(bulk))) bulk <- matrix(bulk, ncol = 1)
  rownames(bulk) <- features
  colnames(bulk) <- sort(unique(labels))
  sweep(bulk, 1, Matrix::rowMeans(data), FUN = "-")
}

match_features <- intersect(VariableFeatures(rep_objects$rep1),
                            VariableFeatures(rep_objects$rep2))
pb1 <- pseudobulk_relative(rep_objects$rep1, match_features)
pb2 <- pseudobulk_relative(rep_objects$rep2, match_features)
cluster_cor <- cor(pb1, pb2, method = "spearman", use = "pairwise.complete.obs")

cor_long <- as.data.frame(as.table(cluster_cor), stringsAsFactors = FALSE)
names(cor_long) <- c("rep1_cluster", "rep2_cluster", "spearman_rho")
cor_long <- as_tibble(cor_long) |>
  mutate(
    rep1_best = rep2_cluster == colnames(cluster_cor)[max.col(cluster_cor,
                                                               ties.method = "first")][match(rep1_cluster, rownames(cluster_cor))],
    rep2_best = rep1_cluster == rownames(cluster_cor)[max.col(t(cluster_cor),
                                                               ties.method = "first")][match(rep2_cluster, colnames(cluster_cor))],
    reciprocal_best = rep1_best & rep2_best
  )

matches <- cor_long |>
  filter(reciprocal_best) |>
  left_join(
    independent_dominant |>
      filter(replicate == "rep1") |>
      dplyr::select(rep1_cluster = independent_cluster,
             rep1_dominant_fine = dominant_fine_cluster,
             rep1_purity = fine_cluster_purity,
             rep1_cells = n_cells),
    by = "rep1_cluster"
  ) |>
  left_join(
    independent_dominant |>
      filter(replicate == "rep2") |>
      dplyr::select(rep2_cluster = independent_cluster,
             rep2_dominant_fine = dominant_fine_cluster,
             rep2_purity = fine_cluster_purity,
             rep2_cells = n_cells),
    by = "rep2_cluster"
  ) |>
  mutate(
    same_fine_cluster = rep1_dominant_fine == rep2_dominant_fine,
    passes_match_threshold = spearman_rho >= MIN_MATCH_RHO &
      rep1_cells >= MIN_REP_CELLS & rep2_cells >= MIN_REP_CELLS,
    supports_integrated_fine_cluster = passes_match_threshold &
      same_fine_cluster & rep1_purity >= MIN_FINE_PURITY &
      rep2_purity >= MIN_FINE_PURITY,
    consensus_cluster = if_else(supports_integrated_fine_cluster,
                                rep1_dominant_fine, NA_character_)
  ) |>
  arrange(desc(supports_integrated_fine_cluster), consensus_cluster)

# -----------------------------------------------------------------------------
# 2. Replicate-concordant markers within independently recovered populations
# -----------------------------------------------------------------------------
replicate_effect <- function(obj, labels, focal_cluster) {
  data <- GetAssayData(obj, assay = "RNA", layer = "data")
  focal <- which(labels == focal_cluster)
  reference <- which(labels != focal_cluster)
  focal_mean <- Matrix::rowMeans(data[, focal, drop = FALSE]) / log(2)
  reference_mean <- Matrix::rowMeans(data[, reference, drop = FALSE]) / log(2)
  focal_pct <- Matrix::rowMeans(data[, focal, drop = FALSE] > 0)
  reference_pct <- Matrix::rowMeans(data[, reference, drop = FALSE] > 0)
  tibble(
    gene = rownames(data),
    n_focal = length(focal),
    mean_log2_difference = focal_mean - reference_mean,
    pct_focal = focal_pct,
    pct_reference = reference_pct,
    pct_difference = focal_pct - reference_pct
  )
}

functional <- read.delim(FUNCTION_FILE, check.names = FALSE, quote = "",
                         na.strings = c("", "NA"))
ortholog_expression <- read.delim(ORTHOLOG_FILE, check.names = FALSE,
                                  quote = "")
ortholog_catalog <- ortholog_expression |>
  dplyr::select(gene, `Gene name`, `Annotation symbol`, `Flybase ID`, `marker of`, ref) |>
  distinct()

consensus_markers <- bind_rows(lapply(seq_len(nrow(recovery_pairs)), function(i) {
  pair <- recovery_pairs[i, ]
  if (!pair$passes_recovery) return(tibble())
  labels1 <- cluster_grid$rep1[[sprintf("%.1f", pair$resolution_rep1)]]
  labels2 <- cluster_grid$rep2[[sprintf("%.1f", pair$resolution_rep2)]]
  e1 <- replicate_effect(rep_objects$rep1, labels1,
                         pair$independent_cluster_rep1) |>
    rename_with(~paste0(.x, "_rep1"), -gene)
  e2 <- replicate_effect(rep_objects$rep2, labels2,
                         pair$independent_cluster_rep2) |>
    rename_with(~paste0(.x, "_rep2"), -gene)
  inner_join(e1, e2, by = "gene") |>
    mutate(
      consensus_cluster = pair$fine_cluster,
      rep1_resolution = pair$resolution_rep1,
      rep2_resolution = pair$resolution_rep2,
      rep1_cluster = pair$independent_cluster_rep1,
      rep2_cluster = pair$independent_cluster_rep2,
      min_log2_difference = pmin(mean_log2_difference_rep1,
                                 mean_log2_difference_rep2),
      min_pct_focal = pmin(pct_focal_rep1, pct_focal_rep2),
      min_pct_difference = pmin(pct_difference_rep1, pct_difference_rep2),
      robust_marker = min_log2_difference >= ROBUST_LOG2_DIFF &
        min_pct_focal >= ROBUST_PCT_FOCAL &
        min_pct_difference >= ROBUST_PCT_DIFF
    )
})) |>
  left_join(functional, by = c("gene" = "Geneid"))

marker_summary <- consensus_markers |>
  group_by(consensus_cluster, rep1_resolution, rep2_resolution,
           rep1_cluster, rep2_cluster) |>
  summarise(
    n_rep1 = dplyr::first(n_focal_rep1),
    n_rep2 = dplyr::first(n_focal_rep2),
    n_robust_markers = sum(robust_marker),
    replicate_effect_rho = cor(mean_log2_difference_rep1,
                               mean_log2_difference_rep2,
                               method = "spearman", use = "complete.obs"),
    .groups = "drop"
  ) |>
  left_join(recovery_pairs |>
              dplyr::select(fine_cluster, passes_recovery,
                            precision_rep1, precision_rep2,
                            recall_rep1, recall_rep2,
                            f1_rep1, f1_rep2),
            by = c("consensus_cluster" = "fine_cluster")) |>
  mutate(shared_population = passes_recovery &
           replicate_effect_rho >= MIN_MATCH_RHO)

shared_fine_clusters <- marker_summary |>
  filter(shared_population) |>
  pull(consensus_cluster) |>
  unique() |>
  sort()

robust_markers <- consensus_markers |>
  filter(robust_marker, consensus_cluster %in% shared_fine_clusters) |>
  arrange(consensus_cluster, desc(min_log2_difference))

ortholog_support <- robust_markers |>
  inner_join(ortholog_catalog, by = "gene", relationship = "many-to-many") |>
  dplyr::select(consensus_cluster, gene, min_log2_difference, min_pct_focal,
         `Gene name`, `Annotation symbol`, `Flybase ID`, `marker of`, ref) |>
  distinct() |>
  arrange(consensus_cluster, desc(min_log2_difference))

# -----------------------------------------------------------------------------
# 3. GO over-representation analysis using only replicate-consensus markers
# -----------------------------------------------------------------------------
background <- rownames(counts)[Matrix::rowSums(counts > 0) >= 3]
term2gene <- functional |>
  dplyr::select(gene = Geneid, GO_term) |>
  filter(!is.na(GO_term), gene %in% background) |>
  separate_rows(GO_term, sep = ",\\s*") |>
  mutate(GO_term = str_trim(GO_term)) |>
  filter(str_detect(GO_term, "^GO:[0-9]{7}$")) |>
  distinct(GO_term, gene)

go_ids <- intersect(unique(term2gene$GO_term), AnnotationDbi::keys(GO.db::GOTERM))
go_terms <- AnnotationDbi::Term(GO.db::GOTERM[go_ids])
go_ontology <- AnnotationDbi::Ontology(GO.db::GOTERM[go_ids])
term_info <- tibble(
  GO_term = names(go_terms),
  Description = unname(go_terms),
  Ontology = unname(go_ontology)
) |>
  filter(!is.na(Description), Ontology %in% c("BP", "MF"))

generic_go <- c(
  "GO:0008150", "GO:0009987", "GO:0008152", "GO:0050789",
  "GO:0065007", "GO:0050896", "GO:0032501", "GO:0032502",
  "GO:0044237", "GO:0044238", "GO:0003674", "GO:0005488",
  "GO:0003824", "GO:0005215", "GO:0005515", "GO:0003676"
)
term_info <- filter(term_info, !GO_term %in% generic_go)
term2gene <- semi_join(term2gene, term_info, by = "GO_term")

parse_ratio <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(z) {
    as.numeric(z[1]) / as.numeric(z[2])
  }, numeric(1))
}

run_go <- function(cluster_id, ontology) {
  t2g <- term2gene |>
    semi_join(filter(term_info, Ontology == ontology), by = "GO_term")
  genes <- robust_markers |>
    filter(consensus_cluster == cluster_id) |>
    pull(gene) |>
    intersect(unique(t2g$gene))
  universe <- intersect(background, unique(t2g$gene))
  if (length(genes) < 10L) return(tibble())
  fit <- clusterProfiler::enricher(
    gene = genes,
    universe = universe,
    TERM2GENE = dplyr::select(t2g, GO_term, gene),
    TERM2NAME = dplyr::select(term_info, GO_term, Description),
    minGSSize = 5,
    maxGSSize = 500,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH"
  )
  if (is.null(fit) || !nrow(as.data.frame(fit))) return(tibble())
  as_tibble(as.data.frame(fit)) |>
    mutate(
      consensus_cluster = as.character(cluster_id),
      ontology = ontology,
      gene_fraction = parse_ratio(GeneRatio),
      background_fraction = parse_ratio(BgRatio),
      fold_enrichment = gene_fraction / background_fraction
    ) |>
    dplyr::select(consensus_cluster, ontology, ID, Description, GeneRatio, BgRatio,
           fold_enrichment, pvalue, p.adjust, qvalue, geneID, Count) |>
    arrange(p.adjust, desc(fold_enrichment))
}

go_all <- bind_rows(lapply(shared_fine_clusters, function(cluster_id) {
  bind_rows(run_go(cluster_id, "BP"), run_go(cluster_id, "MF"))
}))

reduce_terms <- function(x, max_terms = 10L, max_jaccard = 0.70) {
  x <- x |>
    filter(p.adjust < 0.05, Count >= 3, fold_enrichment >= 1.5) |>
    arrange(p.adjust, desc(fold_enrichment), desc(Count))
  if (!nrow(x)) return(x)
  keep <- integer()
  sets <- list()
  for (i in seq_len(nrow(x))) {
    current <- unique(strsplit(x$geneID[i], "/", fixed = TRUE)[[1]])
    overlap <- if (!length(sets)) numeric() else vapply(sets, function(s) {
      length(intersect(current, s)) / length(union(current, s))
    }, numeric(1))
    if (!length(overlap) || all(overlap <= max_jaccard)) {
      keep <- c(keep, i)
      sets[[length(sets) + 1L]] <- current
    }
    if (length(keep) >= max_terms) break
  }
  x[keep, , drop = FALSE]
}

go_top <- if (nrow(go_all)) {
  go_all |>
    group_by(consensus_cluster, ontology) |>
    group_modify(~reduce_terms(.x)) |>
    ungroup()
} else {
  tibble()
}

# -----------------------------------------------------------------------------
# 4. Evidence-led reannotation
# -----------------------------------------------------------------------------
fine_counts <- as_tibble(metadata, rownames = "cell") |>
  count(fine_cluster, replicate, name = "n_cells") |>
  pivot_wider(names_from = replicate, values_from = n_cells, values_fill = 0) |>
  mutate(total_cells = rep1 + rep2)

annotation_table <- tibble::tribble(
  ~fine_cluster, ~consensus_annotation, ~consensus_short_label, ~consensus_confidence,
  ~classical_correspondence, ~selected_marker_evidence,
  ~selected_GO_evidence, ~ortholog_evidence,
  ~interpretation,
  "1", "Cycling/proliferative hemocytes", "Cycling/proliferative", "high for state",
  "Lineage unresolved",
  "Anillin, abnormal spindle, PCNA, Myc, DNA topoisomerase II and RNA-processing genes",
  "Mitotic cell cycle, chromosome and spindle organization, sister-chromatid segregation, cell division and RNA processing",
  "No lineage-defining Drosophila hemocyte ortholog among the consensus markers",
  "A replicate-shared proliferative state, not a morphology-defined lineage and not evidence for a distinct prehemocyte class.",
  "2", "PPO1-high immune-adhesive hemocytes", "PPO1-high immune-adhesive", "moderate",
  "Melanization-associated population; oenocytoid-like correspondence is possible but not exclusive",
  "PPO1, peroxidasin, dual oxidase, PPO-activating factor, integrin beta-PS, integrin alpha-PS3 and paxillin",
  "Biotic-stimulus and innate-immune responses, integrin signaling, actin regulation, adhesion and integrin binding",
  "A PPO ortholog annotated as a crystal-cell marker plus integrin-associated orthologs reported for other Drosophila hemocyte states",
  "A replicate-shared PPO1-high immune and adhesive program. The mixed ortholog evidence does not justify a one-to-one classical lineage label.",
  "3", "PPO-high Notch/FoxO stress-responsive hemocytes", "PPO-high Notch/FoxO",
  "moderate",
  "Oenocytoid-like melanization branch candidate",
  "High shared PPO paralogs with Notch, FoxO, heat-shock proteins, Mlp84B, ficolin and redox-associated genes",
  "Heat response, transition-metal transport and regulation of cell-population proliferation",
  "Pebbled-like and another crystal-cell marker ortholog occur among the consensus markers",
  "A major replicate-shared PPO-rich stress-responsive state. It is related to cluster 6 but is not equivalent to a definitive morphological oenocytoid assignment.",
  "4", "Adhesive/clotting hemocytes", "Adhesive/clotting", "moderate-to-high",
  "Plasmatocyte-like candidate; granulocyte correspondence cannot be excluded",
  "Hemocytin, hemocyte transglutaminases, papilin, laminin subunits, Nimrod/Eater-like, thrombospondin and integrins",
  "Cell adhesion, regulation of locomotion and motility, receptor signaling and calcium-ion binding",
  "Hemocytin and Nimrod/Eater-like orthologs are reported plasmatocyte markers in Drosophila",
  "The strongest lineage-like transcriptomic population: a replicate-shared adhesion, extracellular-matrix and clotting program, provisionally plasmatocyte-like.",
  "5", "Putative prohemocyte-enriched population", "Putative prohemocytes",
  "low-to-moderate",
  "Prohemocyte candidate supported by CPB morphology and immunostaining",
  "No positive marker satisfies the replicate-consensus thresholds; very low RNA and feature counts define the population",
  "No GO test because no qualifying positive marker set exists",
  "No lineage-defining Drosophila ortholog peaks",
  "Present in both samples and supported by morphology, but RNA-low and lacking a reproducible positive marker program; not a transcriptomically proven lineage.",
  "6", "Notch/FoxO/Tob2-high PPO-associated state", "Notch/FoxO/Tob2 PPO",
  "low-to-moderate",
  "Differentiation-associated state within the PPO-rich branch",
  "Tob2, FoxO, Notch, Mth2-like receptor, ficolin, trehalase and catalase; shared PPO paralogs remain highly expressed",
  "cAMP metabolism, lipid response, tissue development and regulation of multicellular-organismal processes",
  "A crystal-cell marker ortholog is present, but the state-specific genes are not lineage-specific markers",
  "A replicate-shared differentiation-associated PPO state related to cluster 3; treated as a state rather than a separate hemocyte type.",
  "7", "Rare PPO/stress antimicrobial state", "Rare PPO/stress", "low",
  "Rare activated state related to PPO-rich hemocytes",
  "PPO1, small heat-shock proteins, ficolin, Mlp84B, cathepsin, mannose receptor and redox-associated genes",
  "Response and defense response to bacteria, general defense response and stress response",
  "A PPO ortholog annotated as a Drosophila crystal-cell marker",
  "A rare but exactly recovered state with a coherent antimicrobial and stress program; rarity and stress sensitivity preclude a separate-lineage claim.",
  "8", "Germline contaminant", "Germline contaminant", "high",
  "Non-hemocyte contaminant",
  "Meiotic and germline-associated transcripts",
  "Not analyzed as a hemocyte population",
  "No hemocyte-lineage interpretation",
  "Replicate-skewed meiotic/germline population excluded from hemocyte interpretation."
) |>
  left_join(fine_counts, by = "fine_cluster") |>
  left_join(
    marker_summary |>
      dplyr::select(fine_cluster = consensus_cluster, n_robust_markers,
             replicate_effect_rho),
    by = "fine_cluster"
  ) |>
  mutate(
    independently_recovered = fine_cluster %in% shared_fine_clusters,
    annotation_basis = case_when(
      fine_cluster == "8" ~ "Excluded contamination",
      fine_cluster == "5" ~
        "Independently recovered morphology-supported candidate; no positive marker program",
      independently_recovered & coalesce(n_robust_markers, 0L) >= 10L ~
        "Independent grid recovery plus replicate-consensus markers",
      independently_recovered ~
        "Independent grid recovery but insufficient positive markers for GO",
      TRUE ~ "Not independently reproduced"
    )
  )

label_map <- setNames(annotation_table$consensus_annotation,
                      annotation_table$fine_cluster)
short_map <- setNames(annotation_table$consensus_short_label,
                      annotation_table$fine_cluster)
combined$consensus_annotation <- unname(label_map[as.character(combined$res_0.2)])
combined$consensus_short_label <- unname(short_map[as.character(combined$res_0.2)])
combined$replicate_consensus <- ifelse(
  as.character(combined$res_0.2) %in% shared_fine_clusters,
  "Independently recovered", "Not independently recovered"
)

# -----------------------------------------------------------------------------
# 5. Tables and figures
# -----------------------------------------------------------------------------
write.table(independent_summary,
            file.path(CONSENSUS_OUT, "independent_cluster_composition.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cor_long, file.path(CONSENSUS_OUT, "independent_cluster_correlations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(matches, file.path(CONSENSUS_OUT, "reciprocal_cluster_matches.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cluster_recovery,
            file.path(CONSENSUS_OUT, "cluster_recovery_resolution_grid.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(recovery_pairs,
            file.path(CONSENSUS_OUT, "best_independent_cluster_recovery.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(robust_markers,
            file.path(CONSENSUS_OUT, "replicate_consensus_markers.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(marker_summary,
            file.path(CONSENSUS_OUT, "replicate_consensus_marker_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ortholog_support,
            file.path(CONSENSUS_OUT, "drosophila_ortholog_support.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(go_all, file.path(CONSENSUS_OUT, "GO_enrichment_all.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(go_top, file.path(CONSENSUS_OUT, "GO_enrichment_top_nonredundant.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(annotation_table, file.path(CONSENSUS_OUT, "cluster_annotations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

workbook <- list(
  annotations = annotation_table,
  best_recovery = recovery_pairs,
  recovery_grid = cluster_recovery,
  reciprocal_matches = matches,
  marker_summary = marker_summary,
  consensus_markers = robust_markers,
  GO_top = go_top,
  GO_all = go_all,
  fly_ortholog_support = ortholog_support,
  independent_composition = independent_summary,
  cluster_correlations = cor_long
)
writexl::write_xlsx(workbook,
                    file.path(CONSENSUS_OUT, "replicate_consensus_GO.xlsx"))

match_labels <- matches |>
  transmute(rep1_cluster, rep2_cluster,
            match_label = if_else(supports_integrated_fine_cluster,
                                  paste0("shared C", consensus_cluster),
                                  "reciprocal, unmapped"))
heat_data <- cor_long |>
  left_join(match_labels, by = c("rep1_cluster", "rep2_cluster")) |>
  mutate(match_label = replace_na(match_label, ""))

p_heat <- ggplot(heat_data,
                 aes(x = rep2_cluster, y = rep1_cluster, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = if_else(reciprocal_best,
                                sprintf("%.2f*", spearman_rho),
                                sprintf("%.2f", spearman_rho))), size = 3) +
  scale_fill_viridis_c(option = "C", limits = c(-1, 1)) +
  labs(x = "Replicate 2 independent cluster",
       y = "Replicate 1 independent cluster",
       fill = "Spearman rho",
       title = "Independent-cluster transcriptomic similarity",
       subtitle = "Asterisk: reciprocal best match") +
  theme_classic(base_size = 10)

recovery_plot_data <- best_recovery |>
  pivot_longer(c(precision, recall), names_to = "metric",
               values_to = "value") |>
  mutate(
    metric = recode(metric, precision = "Purity", recall = "Coverage"),
    fine_cluster = factor(fine_cluster,
                          levels = sort(unique(fine_cluster)))
  )
p_recovery <- ggplot(recovery_plot_data,
                     aes(x = fine_cluster, y = value, colour = metric,
                         shape = replicate,
                         group = interaction(metric, replicate))) +
  geom_hline(yintercept = MIN_RECOVERY_PRECISION, linetype = 2,
             colour = "#777777") +
  geom_line(alpha = 0.45, position = position_dodge(width = 0.28)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.28)) +
  scale_colour_manual(values = c(Purity = "#0072B2", Coverage = "#D55E00")) +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
  labs(x = "Integrated cluster", y = "Best independent recovery",
       colour = NULL, shape = NULL,
       title = "Recovery across the independent resolution grid",
       subtitle = "Dashed line: required purity and coverage") +
  theme_classic(base_size = 10) +
  theme(legend.position = "top")

p_markers <- annotation_table |>
  mutate(n_robust_markers = replace_na(n_robust_markers, 0L),
         fine_cluster = factor(fine_cluster,
                               levels = rev(sort(unique(fine_cluster))))) |>
  ggplot(aes(x = n_robust_markers, y = fine_cluster,
             fill = independently_recovered)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(`TRUE` = "#2B8CBE", `FALSE` = "#BDBDBD")) +
  labs(x = "Replicate-consensus positive markers", y = "Integrated cluster",
       fill = "Independent recovery",
       title = "Evidence available for GO analysis") +
  theme_classic(base_size = 10)

if (nrow(go_top)) {
  go_plot_data <- go_top |>
    group_by(consensus_cluster) |>
    slice_min(p.adjust, n = 4, with_ties = FALSE) |>
    ungroup() |>
    mutate(
      term_label = str_trunc(paste0(Description, " [", ontology, "]"), 52),
      consensus_cluster = factor(consensus_cluster,
                                 levels = sort(unique(consensus_cluster)))
    )
  p_go <- ggplot(go_plot_data,
                 aes(x = consensus_cluster,
                     y = reorder(term_label, -log10(p.adjust)),
                     size = Count, colour = -log10(p.adjust))) +
    geom_point(alpha = 0.9) +
    scale_colour_viridis_c(option = "C", end = 0.9) +
    labs(x = "Consensus cluster", y = NULL, size = "Marker genes",
         colour = expression(-log[10](adjusted~P)),
         title = "GO enrichment of replicate-consensus markers") +
    theme_classic(base_size = 9) +
    theme(axis.text.y = element_text(size = 8))
} else {
  p_go <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No significant GO terms") +
    theme_void()
}

annotation_palette <- c(
  "Cycling/proliferative" = "#0072B2",
  "PPO1-high immune-adhesive" = "#009E73",
  "PPO-high Notch/FoxO" = "#D55E00",
  "Adhesive/clotting" = "#CC79A7",
  "Putative prohemocytes" = "#E69F00",
  "Notch/FoxO/Tob2 PPO" = "#56B4E9",
  "Rare PPO/stress" = "#F0E442",
  "Germline contaminant" = "#7F7F7F"
)
p_umap <- DimPlot(combined, reduction = "umap.harmony",
                  group.by = "consensus_short_label", cols = annotation_palette,
                  pt.size = 0.12, raster = TRUE) +
  labs(title = "Evidence-led CPB hemocyte annotation", colour = NULL,
       x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 10) +
  theme(legend.text = element_text(size = 7))

ggsave(file.path(FIG_OUT, "independent_cluster_similarity.pdf"), p_heat,
       width = 7.2, height = 5.8)
ggsave(file.path(FIG_OUT, "independent_cluster_recovery_grid.pdf"),
       p_recovery, width = 7.2, height = 5.2)
ggsave(file.path(FIG_OUT, "consensus_marker_counts.pdf"), p_markers,
       width = 7.2, height = 4.8)
ggsave(file.path(FIG_OUT, "consensus_GO_dotplot.pdf"), p_go,
       width = 10, height = 7.5)
ggsave(file.path(FIG_OUT, "reannotated_UMAP.pdf"), p_umap,
       width = 10.5, height = 6.5)
ggsave(file.path(FIG_OUT, "replicate_consensus_summary.pdf"),
       (p_recovery | p_markers) / (p_umap | p_go) +
         plot_annotation(tag_levels = "A"),
       width = 16, height = 12, limitsize = FALSE)

saveRDS(combined, file.path(OUT,
                            "combined_consensus.rds"))
capture.output(sessionInfo(), file = file.path(CONSENSUS_OUT, "sessionInfo.txt"))

message("Shared integrated clusters: ",
        paste(shared_fine_clusters, collapse = ", "))
message("Results: ", CONSENSUS_OUT)
