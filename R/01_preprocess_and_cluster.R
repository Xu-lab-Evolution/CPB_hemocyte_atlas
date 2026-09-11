# Preprocessing, integration, clustering, and resolution assessment.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(writexl)
})

if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
  options(future.globals.maxSize = 4 * 1024^3) 
}

source(file.path("R", "config.R"))

set.seed(SEED)

# -----------------------------------------------------------------------------
# 0. Output paths and analysis settings
# -----------------------------------------------------------------------------
QC  <- file.path(OUT, "QC")
FIG <- file.path(OUT, "figures")
MRK <- file.path(OUT, "markers")
VAL <- file.path(OUT, "resolution_validation")
ANN <- file.path(OUT, "annotation_transfer")
invisible(lapply(c(OUT, QC, FIG, MRK, VAL, ANN), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

# QC thresholds: permissive before scDblFinder; final after doublet removal.
MIN_FEATURES_PRE <- 100
MIN_COUNTS_PRE   <- 200
MAX_MT_PRE       <- 30
MIN_FEATURES     <- 200
MIN_COUNTS       <- 500
MAX_MT           <- 15

# The elbow plot supports using the first 20 principal components.
N_HVG <- 3000
N_PCS <- 50
DIMS  <- 1:20
K_NN  <- 20

# The complete grid is retained for resolution diagnostics.
RES_GRID <- seq(0.1, 1.0, by = 0.1)

# Subsampling stability is estimated from 50 independent 80% cell samples.
N_SUBSAMPLE    <- 50
SUBSAMPLE_FRAC <- 0.80

# Tang et al. use Jaccard 0.75 as a useful stable-cluster cutoff. Rather than
# requiring every cluster to be stable, use the proportion of CELLS contained
# in stable clusters. This avoids one tiny unstable cluster forcing severe
# under-clustering.
STABLE_JACCARD      <- 0.75
HIGH_STABLE_JACCARD <- 0.85
MIN_MEDIAN_JACCARD  <- 0.75
MIN_STABLE_CELL_FRAC <- 0.90
MIN_CLUSTER_N_WARN  <- 20      # warning/report only; not a hard rejection rule

# scClustViz-style local separation. These are explicit effect-size filters, not
# universal constants. They prevent thousands of trivial but significant genes
# in this large dataset from defining cluster separability.
LOCAL_PADJ     <- 0.01
LOCAL_LOG2FC   <- 0.50
LOCAL_PCT_DIFF <- 0.15
MIN_LOCAL_GENES <- 5
MIN_LOCAL_SEP_FRAC <- 0.80

# One-vs-rest marker quality used as an additional resolution diagnostic.
RES_MARKER_LOG2FC   <- 0.50
RES_MARKER_PCT1     <- 0.20
RES_MARKER_PCT_DIFF <- 0.15
MIN_STRONG_MARKERS_PER_CLUSTER <- 5
MIN_MARKER_RICH_FRAC <- 0.80

# Replicate directional support for local DE. We do not run a second set of
# pseudoreplicated significance tests. Instead, a globally meaningful DE gene
# must show the same mean-expression direction in rep1 and rep2.
MIN_PAIR_CELLS <- 15
REP_MEAN_DIFF  <- 0.05
MIN_REP_GENES  <- 1
MIN_REP_SUPPORTED_FRAC <- 0.80

# Balanced cell cap for resolution-validation DE only. Final marker discovery
# below uses the complete selected dataset.
DE_MAX_CELLS_PER_IDENT <- 750

S1_DIR <- input_file("filtered_feature_bc_matrix_s1", directory = TRUE)
S2_DIR <- input_file("filtered_feature_bc_matrix_s2", directory = TRUE)
GTF_FILE <- input_file("CPBv5_single_mod_mito_prot_coding.gtf")
ANNO_FILE <- input_file("LdecV5_functional_annotation.txt")
CROSSWALK_FILE <- input_file("table_LDECv5_vs_LdNA.txt")
DROSO_MARKER_FILE <- input_file("marker_genes_droso_cpb.txt")
CELL_CYCLE_FILE <- input_file("cpb_droso_cellcycle.tsv")

# Record key parameters for reproducibility.
params <- data.frame(
  parameter = c(
    "DIMS", "K_NN", "RES_GRID", "N_SUBSAMPLE", "SUBSAMPLE_FRAC",
    "STABLE_JACCARD", "HIGH_STABLE_JACCARD", "MIN_MEDIAN_JACCARD",
    "MIN_STABLE_CELL_FRAC", "MIN_CLUSTER_N_WARN",
    "LOCAL_PADJ", "LOCAL_LOG2FC", "LOCAL_PCT_DIFF", "MIN_LOCAL_GENES",
    "MIN_LOCAL_SEP_FRAC", "RES_MARKER_LOG2FC", "RES_MARKER_PCT1",
    "RES_MARKER_PCT_DIFF", "MIN_STRONG_MARKERS_PER_CLUSTER",
    "MIN_MARKER_RICH_FRAC", "MIN_PAIR_CELLS", "REP_MEAN_DIFF",
    "MIN_REP_GENES", "MIN_REP_SUPPORTED_FRAC", "DE_MAX_CELLS_PER_IDENT"
  ),
  value = c(
    paste(range(DIMS), collapse = ":"), K_NN, paste(RES_GRID, collapse = ","),
    N_SUBSAMPLE, SUBSAMPLE_FRAC, STABLE_JACCARD, HIGH_STABLE_JACCARD,
    MIN_MEDIAN_JACCARD, MIN_STABLE_CELL_FRAC, MIN_CLUSTER_N_WARN,
    LOCAL_PADJ, LOCAL_LOG2FC, LOCAL_PCT_DIFF, MIN_LOCAL_GENES,
    MIN_LOCAL_SEP_FRAC, RES_MARKER_LOG2FC, RES_MARKER_PCT1,
    RES_MARKER_PCT_DIFF, MIN_STRONG_MARKERS_PER_CLUSTER, MIN_MARKER_RICH_FRAC,
    MIN_PAIR_CELLS, REP_MEAN_DIFF, MIN_REP_GENES, MIN_REP_SUPPORTED_FRAC,
    DE_MAX_CELLS_PER_IDENT
  )
)
write.table(params, file.path(OUT, "analysis_parameters.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# -----------------------------------------------------------------------------
# 1. Read data and calculate mitochondrial fraction
# -----------------------------------------------------------------------------
gtf <- read.delim(GTF_FILE, header = FALSE, comment.char = "#", quote = "")
mt_attr <- gtf |>
  dplyr::filter(V1 == "MZ189364", V3 == "CDS") |>
  dplyr::pull(V9)
mt_genes <- unique(stats::na.omit(
  stringr::str_match(mt_attr, 'gene_id "([^"]+)"')[, 2]
))

make_object <- function(path, sample) {
  x <- CreateSeuratObject(Read10X(path), project = sample)
  x$replicate <- sample
  x[["percent.mt"]] <- PercentageFeatureSet(
    x, features = intersect(mt_genes, rownames(x)), assay = "RNA"
  )
  x
}

s1 <- make_object(S1_DIR, "rep1")
s2 <- make_object(S2_DIR, "rep2")

# -----------------------------------------------------------------------------
# 2. QC and doublet detection
# -----------------------------------------------------------------------------
get_meta <- function(x, stage) {
  dplyr::mutate(x[[]], cell = rownames(x[[]]), stage = stage)
}

qc_thresholds <- tibble::tribble(
  ~metric,        ~value,           ~threshold,
  "nFeature_RNA", MIN_FEATURES_PRE, "permissive",
  "nFeature_RNA", MIN_FEATURES,     "final",
  "nCount_RNA",   MIN_COUNTS_PRE,   "permissive",
  "nCount_RNA",   MIN_COUNTS,       "final",
  "percent.mt",   MAX_MT_PRE,       "permissive",
  "percent.mt",   MAX_MT,           "final"
)

plot_qc <- function(df, title, thresholds = TRUE) {
  z <- df |>
    dplyr::select(replicate, nFeature_RNA, nCount_RNA, percent.mt) |>
    tidyr::pivot_longer(-replicate, names_to = "metric", values_to = "value")
  p <- ggplot(z, aes(replicate, value, fill = replicate)) +
    geom_violin(trim = FALSE, alpha = 0.75) +
    facet_wrap(~metric, scales = "free_y", nrow = 1) +
    labs(title = title, x = NULL, y = NULL, fill = "Sample") +
    theme_classic() + theme(legend.position = "top")
  if (thresholds) {
    p <- p + geom_hline(data = qc_thresholds,
                        aes(yintercept = value, linetype = threshold),
                        inherit.aes = FALSE)
  }
  p
}

qc_input <- dplyr::bind_rows(get_meta(s1, "input"), get_meta(s2, "input"))
p_qc_input <- plot_qc(qc_input, "QC before additional filtering")
p_qc_scatter <- ggplot(qc_input,
                        aes(nCount_RNA, nFeature_RNA, colour = replicate)) +
  geom_point(size = 0.35, alpha = 0.35) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Input complexity", x = "UMIs (log10)",
       y = "Detected genes (log10)", colour = "Sample") +
  theme_classic() + theme(legend.position = "top")

prefilter <- function(x) {
  subset(x, subset = nFeature_RNA >= MIN_FEATURES_PRE &
                    nCount_RNA >= MIN_COUNTS_PRE & percent.mt < MAX_MT_PRE)
}
s1 <- prefilter(s1)
s2 <- prefilter(s2)
qc_prefilter <- dplyr::bind_rows(get_meta(s1, "prefilter"),
                                 get_meta(s2, "prefilter"))

run_doublets <- function(x, seed) {
  set.seed(seed)
  sce <- as.SingleCellExperiment(x, assay = "RNA")
  sce <- scDblFinder(sce, verbose = TRUE,
                     BPPARAM = BiocParallel::SerialParam(RNGseed = seed))
  x$scDblFinder.class <- colData(sce)$scDblFinder.class
  x$scDblFinder.score <- colData(sce)$scDblFinder.score
  x
}
s1 <- run_doublets(s1, SEED)
s2 <- run_doublets(s2, SEED + 1)

qc_doublet <- dplyr::bind_rows(get_meta(s1, "doublet_called"),
                               get_meta(s2, "doublet_called"))
doublet_summary <- qc_doublet |>
  dplyr::count(replicate, scDblFinder.class, name = "n") |>
  dplyr::group_by(replicate) |>
  dplyr::mutate(percent = 100 * n / sum(n)) |>
  dplyr::ungroup()
write.table(doublet_summary, file.path(QC, "doublet_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

p_doublet <- ggplot(doublet_summary,
                    aes(scDblFinder.class, percent, fill = replicate)) +
  geom_col(position = "dodge") +
  labs(title = "scDblFinder calls", x = NULL, y = "% cells", fill = "Sample") +
  theme_classic() + theme(legend.position = "top")

combined <- merge(s1, y = s2, add.cell.ids = c("rep1", "rep2"),
                  project = "CPBv5_hemocytes")
combined <- subset(combined, subset = scDblFinder.class == "singlet" &
                                      nFeature_RNA >= MIN_FEATURES &
                                      nCount_RNA >= MIN_COUNTS &
                                      percent.mt < MAX_MT)
qc_final <- get_meta(combined, "final")
p_qc_final <- plot_qc(qc_final,
                       "QC after doublet removal and final filtering", FALSE)

count_stage <- function(df, stage) {
  df |> dplyr::count(replicate, name = "n_cells") |>
    dplyr::mutate(stage = stage)
}
retention <- dplyr::bind_rows(
  count_stage(qc_input, "01_input"),
  count_stage(qc_prefilter, "02_after_permissive_QC"),
  count_stage(dplyr::filter(qc_doublet, scDblFinder.class == "singlet"),
              "03_singlets"),
  count_stage(qc_final, "04_final_QC")
)
input_n <- retention |>
  dplyr::filter(stage == "01_input") |>
  dplyr::select(replicate, n_input = n_cells)
retention <- retention |>
  dplyr::left_join(input_n, by = "replicate") |>
  dplyr::mutate(
    percent_input = 100 * n_cells / n_input,
    stage = factor(stage, levels = c("01_input", "02_after_permissive_QC",
                                     "03_singlets", "04_final_QC"))
  ) |>
  dplyr::arrange(replicate, stage)
write.table(retention, file.path(QC, "cell_retention.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

p_retention <- ggplot(retention,
                      aes(stage, percent_input, group = replicate,
                          colour = replicate)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  labs(title = "Cell retention", x = NULL, y = "% of input cells",
       colour = "Sample") +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 25, hjust = 1))

qc_summary <- dplyr::bind_rows(qc_input, qc_final) |>
  dplyr::group_by(stage, replicate) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    median_features = median(nFeature_RNA),
    median_UMIs = median(nCount_RNA),
    median_percent_mt = median(percent.mt),
    .groups = "drop"
  )
write.table(qc_summary, file.path(QC, "QC_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(QC, "QC_report.pdf"), width = 11, height = 6.5, onefile = TRUE)
print(p_qc_input)
print(p_qc_scatter)
print(p_doublet)
print(p_qc_final)
print(p_retention)
dev.off()

# -----------------------------------------------------------------------------
# 3. Gene filtering, normalization, PCA and Harmony
# -----------------------------------------------------------------------------
combined <- JoinLayers(combined)
counts <- GetAssayData(combined, assay = "RNA", layer = "counts")
keep <- Matrix::rowSums(counts > 0) >= 3
combined <- subset(combined, features = rownames(counts)[keep])

combined[["RNA"]] <- split(combined[["RNA"]], f = combined$replicate)
combined <- NormalizeData(combined, verbose = FALSE)
combined <- FindVariableFeatures(combined, nfeatures = N_HVG, verbose = FALSE)
combined <- ScaleData(combined, features = VariableFeatures(combined),
                      verbose = FALSE)
combined <- RunPCA(combined, npcs = N_PCS, verbose = FALSE)

ggsave(file.path(FIG, "PCA_elbow.pdf"),
       ElbowPlot(combined, ndims = N_PCS), width = 6, height = 4.5)

combined <- RunUMAP(combined, reduction = "pca", dims = DIMS,
                    reduction.name = "umap.unintegrated",
                    seed.use = SEED, verbose = FALSE)
ggsave(file.path(FIG, "UMAP_unintegrated_samples.pdf"),
       DimPlot(combined, reduction = "umap.unintegrated",
               group.by = "replicate"), width = 6, height = 5)

combined <- IntegrateLayers(combined, method = HarmonyIntegration,
                            orig.reduction = "pca",
                            new.reduction = "harmony", verbose = FALSE)
combined <- RunUMAP(combined, reduction = "harmony", dims = DIMS,
                    reduction.name = "umap.harmony",
                    seed.use = SEED, verbose = FALSE)
ggsave(file.path(FIG, "UMAP_Harmony_samples.pdf"),
       DimPlot(combined, reduction = "umap.harmony", group.by = "replicate"),
       width = 6, height = 5)

# -----------------------------------------------------------------------------
# 4. Full-data clustering across the resolution grid
# -----------------------------------------------------------------------------
algorithm <- if (requireNamespace("leidenbase", quietly = TRUE)) 4 else 1
res_key <- function(r) sprintf("%.1f", r)
res_col <- function(r) paste0("res_", res_key(r))

harmony <- Embeddings(combined, "harmony")[, DIMS, drop = FALSE]
full_graph <- FindNeighbors(harmony, k.param = K_NN, compute.SNN = TRUE,
                            return.neighbor = FALSE, verbose = FALSE)

cluster_graph <- function(snn, resolution, seed) {
  z <- FindClusters(snn, resolution = resolution, algorithm = algorithm,
                    random.seed = seed, verbose = FALSE)
  labs <- as.character(z[, 1])
  names(labs) <- rownames(z)
  labs
}

full_labels <- list()
for (r in RES_GRID) {
  key <- res_key(r)
  full_labels[[key]] <- cluster_graph(full_graph$snn, r, SEED)
  combined[[res_col(r)]] <- full_labels[[key]][Cells(combined)]
}

# Diagnostic resolution grid, replacing clustree (avoids ggplot2 compatibility
# problems and shows the actual cluster geometry directly).
p_grid <- lapply(RES_GRID, function(r) {
  DimPlot(combined, reduction = "umap.harmony", group.by = res_col(r),
          label = TRUE, repel = TRUE, raster = TRUE) +
    ggtitle(paste0("resolution = ", r)) + NoLegend()
})
ggsave(file.path(FIG, "UMAP_resolution_grid.pdf"),
       patchwork::wrap_plots(p_grid, ncol = 3), width = 15, height = 16)

# -----------------------------------------------------------------------------
# 5. Subsampling/Jaccard cluster stability
# -----------------------------------------------------------------------------
# For each reference cluster, match the best-overlapping cluster after 80%
# subsampling. The Jaccard is calculated on cells present in the subsample, so
# the omitted 20% do not artificially reduce the maximum possible score.
max_jaccard_by_reference_cluster <- function(reference, test) {
  reference <- as.character(reference)
  test <- as.character(test)
  tab <- table(reference, test)
  n_ref <- rowSums(tab)
  n_test <- colSums(tab)

  dplyr::bind_rows(lapply(seq_len(nrow(tab)), function(i) {
    inter <- tab[i, ]
    jac <- inter / (n_ref[i] + n_test - inter)
    tibble::tibble(cluster = rownames(tab)[i], jaccard = max(jac))
  }))
}

message("Running ", N_SUBSAMPLE, " subsampling iterations across ",
        length(RES_GRID), " resolutions...")
stability_runs <- vector("list", N_SUBSAMPLE * length(RES_GRID))
k <- 1L

for (b in seq_len(N_SUBSAMPLE)) {
  set.seed(SEED + b)
  idx <- sort(sample.int(nrow(harmony),
                         size = floor(SUBSAMPLE_FRAC * nrow(harmony)),
                         replace = FALSE))
  emb_sub <- harmony[idx, , drop = FALSE]
  g_sub <- FindNeighbors(emb_sub, k.param = K_NN, compute.SNN = TRUE,
                         return.neighbor = FALSE, verbose = FALSE)

  for (r in RES_GRID) {
    key <- res_key(r)
    sub_lab <- cluster_graph(g_sub$snn, r, SEED + b)
    ref_lab <- full_labels[[key]][names(sub_lab)]
    stability_runs[[k]] <- max_jaccard_by_reference_cluster(ref_lab, sub_lab) |>
      dplyr::mutate(resolution = r, iteration = b)
    k <- k + 1L
  }
  if (b %% 10 == 0) message("  completed ", b, "/", N_SUBSAMPLE)
}

stability_long <- dplyr::bind_rows(stability_runs)
write.table(stability_long, file.path(VAL, "subsampling_jaccard_all.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cluster_stability <- stability_long |>
  dplyr::group_by(resolution, cluster) |>
  dplyr::summarise(
    median_jaccard = median(jaccard),
    q10_jaccard = stats::quantile(jaccard, 0.10),
    q90_jaccard = stats::quantile(jaccard, 0.90),
    .groups = "drop"
  )

# Add full-data cluster sizes so that stability can be summarized by the
# fraction of CELLS in stable clusters, as recommended by Tang et al.
cluster_sizes <- dplyr::bind_rows(lapply(RES_GRID, function(r) {
  labs <- full_labels[[res_key(r)]][Cells(combined)]
  tibble::tibble(cluster = as.character(names(table(labs))),
                 n_cells = as.integer(table(labs)),
                 resolution = r)
}))
cluster_stability <- cluster_stability |>
  dplyr::left_join(cluster_sizes, by = c("resolution", "cluster")) |>
  dplyr::mutate(
    stable_075 = median_jaccard >= STABLE_JACCARD,
    stable_085 = median_jaccard >= HIGH_STABLE_JACCARD
  )
write.table(cluster_stability, file.path(VAL, "cluster_stability.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

balanced_silhouette <- function(labels, emb, n_per_cluster = 100) {
  if (length(unique(labels)) < 2) return(NA_real_)
  set.seed(SEED)
  idx <- unlist(lapply(split(seq_along(labels), labels), function(i) {
    sample(i, min(length(i), n_per_cluster))
  }), use.names = FALSE)
  lab <- factor(labels[idx])
  sil <- cluster::silhouette(as.integer(lab),
                             stats::dist(emb[idx, , drop = FALSE]))
  mean(tapply(sil[, "sil_width"], lab, mean))
}

resolution_stability <- dplyr::bind_rows(lapply(RES_GRID, function(r) {
  key <- res_key(r)
  labs <- full_labels[[key]][Cells(combined)]
  n_tab <- table(labs)
  stab <- dplyr::filter(cluster_stability, resolution == r)
  n_total <- sum(stab$n_cells)

  tibble::tibble(
    resolution = r,
    n_clusters = length(n_tab),
    min_cluster_n = min(n_tab),
    n_clusters_lt_warn = sum(n_tab < MIN_CLUSTER_N_WARN),
    median_cluster_jaccard = median(stab$median_jaccard),
    min_cluster_jaccard = min(stab$median_jaccard),
    stable_cell_fraction_075 = sum(stab$n_cells[stab$stable_075]) / n_total,
    stable_cell_fraction_085 = sum(stab$n_cells[stab$stable_085]) / n_total,
    stable_cluster_fraction_075 = mean(stab$stable_075),
    balanced_silhouette = balanced_silhouette(labs, harmony)
  )
})) |>
  dplyr::mutate(
    stability_pass = median_cluster_jaccard >= MIN_MEDIAN_JACCARD &
      stable_cell_fraction_075 >= MIN_STABLE_CELL_FRAC
  )

write.table(resolution_stability,
            file.path(VAL, "resolution_stability_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Candidate resolutions: those supported by cell-weighted stability plus the
# immediately finer boundary, which tells us where biological separability
# starts to deteriorate.
stable_res <- resolution_stability |>
  dplyr::filter(stability_pass) |>
  dplyr::pull(resolution)

if (length(stable_res)) {
  next_finer <- RES_GRID[RES_GRID > max(stable_res)]
  next_finer <- if (length(next_finer)) min(next_finer) else numeric()
  BIO_RES_GRID <- sort(unique(c(stable_res, next_finer)))
} else {
  BIO_RES_GRID <- resolution_stability |>
    dplyr::arrange(dplyr::desc(stable_cell_fraction_075),
                   dplyr::desc(median_cluster_jaccard),
                   dplyr::desc(balanced_silhouette)) |>
    dplyr::slice_head(n = 4) |>
    dplyr::pull(resolution) |>
    sort()
}
writeLines(paste(BIO_RES_GRID, collapse = ","),
           file.path(VAL, "biological_validation_resolutions.txt"))
message("Biological validation will be run for resolutions: ",
        paste(BIO_RES_GRID, collapse = ", "))

# -----------------------------------------------------------------------------
# 6. scClustViz-style transcriptional validation of candidate resolutions
# -----------------------------------------------------------------------------
combined <- JoinLayers(combined)
DefaultAssay(combined) <- "RNA"

# Normalized RNA matrix used only for a directional replicate check. This avoids
# repeated per-replicate significance tests and the associated future/memory
# overhead. With two replicates this is a reproducibility check, not inference.
rna_data <- GetAssayData(combined, assay = "RNA", layer = "data")
replicate_vec <- stats::setNames(as.character(combined$replicate), Cells(combined))

replicate_direction_support <- function(genes, cluster, neighbor, labels) {
  genes <- intersect(genes, rownames(rna_data))
  if (!length(genes)) return(character())

  diffs <- list()
  for (rp in sort(unique(replicate_vec))) {
    c1 <- names(labels)[labels == cluster & replicate_vec[names(labels)] == rp]
    c2 <- names(labels)[labels == neighbor & replicate_vec[names(labels)] == rp]
    if (length(c1) < MIN_PAIR_CELLS || length(c2) < MIN_PAIR_CELLS) {
      return(character())
    }
    diffs[[rp]] <- Matrix::rowMeans(rna_data[genes, c1, drop = FALSE]) -
      Matrix::rowMeans(rna_data[genes, c2, drop = FALSE])
  }

  dm <- do.call(cbind, diffs)
  genes[apply(dm, 1, function(x) all(is.finite(x) & x >= REP_MEAN_DIFF))]
}

# Pairwise DE is used to define each cluster's nearest TRANSCRIPTOMIC neighbour:
# the other cluster against which it has the fewest meaningful positive DE
# genes. This follows scClustViz more closely than centroid distance in Harmony.
pairwise_local_de <- function(obj, labels, resolution) {
  clusters <- sort(unique(as.character(labels)))
  Idents(obj) <- factor(labels[Cells(obj)])
  pair_summary <- list()
  gene_results <- list()
  q <- 1L
  gq <- 1L

  if (length(clusters) < 2) return(list(summary = tibble::tibble(), genes = tibble::tibble()))

  for (i in seq_len(length(clusters) - 1L)) {
    for (j in (i + 1L):length(clusters)) {
      a <- clusters[i]
      b <- clusters[j]
      de <- FindMarkers(
        obj, ident.1 = a, ident.2 = b, assay = "RNA", slot = "data",
        test.use = "wilcox", min.pct = 0.05, logfc.threshold = 0,
        max.cells.per.ident = DE_MAX_CELLS_PER_IDENT, random.seed = SEED,
        densify = FALSE, verbose = FALSE
      )
      if (!nrow(de)) next
      de$gene <- rownames(de)

      de_a <- de |>
        dplyr::filter(p_val_adj < LOCAL_PADJ,
                      avg_log2FC >= LOCAL_LOG2FC,
                      (pct.1 - pct.2) >= LOCAL_PCT_DIFF) |>
        dplyr::mutate(cluster = a, neighbor = b,
                      cluster_log2FC = avg_log2FC,
                      cluster_pct_diff = pct.1 - pct.2)
      de_b <- de |>
        dplyr::filter(p_val_adj < LOCAL_PADJ,
                      avg_log2FC <= -LOCAL_LOG2FC,
                      (pct.2 - pct.1) >= LOCAL_PCT_DIFF) |>
        dplyr::mutate(cluster = b, neighbor = a,
                      cluster_log2FC = -avg_log2FC,
                      cluster_pct_diff = pct.2 - pct.1)

      sup_a <- replicate_direction_support(de_a$gene, a, b, labels)
      sup_b <- replicate_direction_support(de_b$gene, b, a, labels)

      if (nrow(de_a)) {
        gene_results[[gq]] <- de_a |>
          dplyr::mutate(resolution = resolution,
                        replicate_supported = gene %in% sup_a)
        gq <- gq + 1L
      }
      if (nrow(de_b)) {
        gene_results[[gq]] <- de_b |>
          dplyr::mutate(resolution = resolution,
                        replicate_supported = gene %in% sup_b)
        gq <- gq + 1L
      }

      pair_summary[[q]] <- tibble::tibble(
        resolution = resolution, cluster = a, neighbor = b,
        n_meaningful_DE = nrow(de_a), n_rep_supported = length(sup_a)
      )
      q <- q + 1L
      pair_summary[[q]] <- tibble::tibble(
        resolution = resolution, cluster = b, neighbor = a,
        n_meaningful_DE = nrow(de_b), n_rep_supported = length(sup_b)
      )
      q <- q + 1L
    }
  }

  list(summary = dplyr::bind_rows(pair_summary),
       genes = dplyr::bind_rows(gene_results))
}

# One-vs-rest marker quality for a candidate resolution. This is the second
# scClustViz-like criterion: most clusters should have multiple useful markers.
resolution_marker_quality <- function(obj, labels, resolution) {
  Idents(obj) <- factor(labels[Cells(obj)])
  mk <- FindAllMarkers(
    obj, assay = "RNA", slot = "data", only.pos = TRUE, test.use = "wilcox",
    min.pct = 0.10, logfc.threshold = 0.10,
    max.cells.per.ident = DE_MAX_CELLS_PER_IDENT,
    random.seed = SEED, verbose = FALSE
  )
  if (!nrow(mk)) {
    return(list(markers = tibble::tibble(), summary = tibble::tibble()))
  }
  mk <- mk |>
    dplyr::mutate(cluster = as.character(cluster),
                  pct_diff = pct.1 - pct.2,
                  strong_marker = p_val_adj < LOCAL_PADJ &
                    avg_log2FC >= RES_MARKER_LOG2FC &
                    pct.1 >= RES_MARKER_PCT1 &
                    pct_diff >= RES_MARKER_PCT_DIFF,
                  resolution = resolution)

  clusters <- sort(unique(as.character(labels)))
  per_cluster <- tibble::tibble(cluster = clusters) |>
    dplyr::left_join(
      mk |>
        dplyr::group_by(cluster) |>
        dplyr::summarise(n_strong_markers = sum(strong_marker), .groups = "drop"),
      by = "cluster"
    ) |>
    dplyr::mutate(n_strong_markers = tidyr::replace_na(n_strong_markers, 0L),
                  resolution = resolution)

  list(markers = mk, summary = per_cluster)
}

local_all <- list()
local_gene_all <- list()
res_marker_all <- list()
res_marker_counts_all <- list()
li <- gi <- mi <- mci <- 1L

for (r in BIO_RES_GRID) {
  message("Biological resolution validation: ", r)
  labels <- full_labels[[res_key(r)]][Cells(combined)]

  loc <- pairwise_local_de(combined, labels, r)
  if (nrow(loc$summary)) {
    # The nearest transcriptomic neighbour is the comparison with the FEWEST
    # meaningful positive DE genes for each focal cluster.
    nearest <- loc$summary |>
      dplyr::group_by(cluster) |>
      dplyr::arrange(n_meaningful_DE, dplyr::desc(n_rep_supported), neighbor) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup()
    local_all[[li]] <- nearest
    li <- li + 1L
  }
  if (nrow(loc$genes)) {
    local_gene_all[[gi]] <- loc$genes
    gi <- gi + 1L
  }

  mq <- resolution_marker_quality(combined, labels, r)
  if (nrow(mq$markers)) {
    res_marker_all[[mi]] <- mq$markers
    mi <- mi + 1L
  }
  if (nrow(mq$summary)) {
    res_marker_counts_all[[mci]] <- mq$summary
    mci <- mci + 1L
  }
}

nearest_local <- dplyr::bind_rows(local_all)
nearest_local_genes <- dplyr::bind_rows(local_gene_all)
resolution_marker_candidates <- dplyr::bind_rows(res_marker_all)
resolution_marker_counts <- dplyr::bind_rows(res_marker_counts_all)

write.table(nearest_local, file.path(VAL, "nearest_transcriptomic_neighbor.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(nearest_local_genes,
            file.path(VAL, "nearest_neighbor_meaningful_DE_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(resolution_marker_counts,
            file.path(VAL, "strong_marker_counts_by_resolution_cluster.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

local_summary <- nearest_local |>
  dplyr::group_by(resolution) |>
  dplyr::summarise(
    min_local_DE = min(n_meaningful_DE),
    median_local_DE = median(n_meaningful_DE),
    frac_clusters_locally_separable = mean(n_meaningful_DE >= MIN_LOCAL_GENES),
    frac_clusters_rep_supported = mean(n_rep_supported >= MIN_REP_GENES),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    local_sep_pass = frac_clusters_locally_separable >= MIN_LOCAL_SEP_FRAC,
    replicate_pass = frac_clusters_rep_supported >= MIN_REP_SUPPORTED_FRAC
  )

marker_quality_summary <- resolution_marker_counts |>
  dplyr::group_by(resolution) |>
  dplyr::summarise(
    min_strong_markers = min(n_strong_markers),
    median_strong_markers = median(n_strong_markers),
    frac_clusters_marker_rich = mean(n_strong_markers >= MIN_STRONG_MARKERS_PER_CLUSTER),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    marker_quality_pass = frac_clusters_marker_rich >= MIN_MARKER_RICH_FRAC
  )

# -----------------------------------------------------------------------------
# 7. Resolution recommendation
# -----------------------------------------------------------------------------
resolution_summary <- resolution_stability |>
  dplyr::left_join(local_summary, by = "resolution") |>
  dplyr::left_join(marker_quality_summary, by = "resolution") |>
  dplyr::mutate(
    full_pass = stability_pass &
      dplyr::coalesce(local_sep_pass, FALSE) &
      dplyr::coalesce(marker_quality_pass, FALSE) &
      dplyr::coalesce(replicate_pass, FALSE),
    # Transparent descriptive score for fallback/ranking only; it is not a
    # statistical objective function.
    diagnostic_score = stable_cell_fraction_075 *
      pmin(median_cluster_jaccard / HIGH_STABLE_JACCARD, 1) *
      dplyr::coalesce(frac_clusters_locally_separable, 0) *
      dplyr::coalesce(frac_clusters_marker_rich, 0) *
      dplyr::coalesce(frac_clusters_rep_supported, 0)
  )

passing <- dplyr::filter(resolution_summary, full_pass)
if (nrow(passing)) {
  SELECTED_RES <- max(passing$resolution)
  selection_status <- "passes_stability_localDE_marker_quality_replicate_support"
} else {
  # Prefer the highest resolution satisfying stability and the largest number of
  # biological checks. This is explicitly flagged rather than presented as an
  # automatic optimum.
  pick <- resolution_summary |>
    dplyr::filter(resolution %in% BIO_RES_GRID) |>
    dplyr::mutate(
      n_checks_passed = as.integer(stability_pass) +
        as.integer(dplyr::coalesce(local_sep_pass, FALSE)) +
        as.integer(dplyr::coalesce(marker_quality_pass, FALSE)) +
        as.integer(dplyr::coalesce(replicate_pass, FALSE))
    ) |>
    dplyr::arrange(dplyr::desc(n_checks_passed),
                   dplyr::desc(diagnostic_score),
                   dplyr::desc(resolution))
  SELECTED_RES <- pick$resolution[1]
  selection_status <- "fallback_best_supported_candidate_manual_review_required"
  warning("No resolution passed every criterion; inspect resolution_selection_summary.tsv before accepting the fallback.")
}

resolution_summary <- resolution_summary |>
  dplyr::mutate(recommended = resolution == SELECTED_RES,
                selection_status = ifelse(recommended, selection_status, ""))
write.table(resolution_summary,
            file.path(VAL, "resolution_selection_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(
  c(paste0("Recommended resolution: ", SELECTED_RES),
    paste0("Selection status: ", selection_status),
    "Primary stability metric: fraction of cells in clusters with median Jaccard >= 0.75.",
    paste0("Stability requirement: >=", MIN_STABLE_CELL_FRAC * 100,
           "% of cells stable and median cluster Jaccard >=", MIN_MEDIAN_JACCARD, "."),
    "Biological validation: scClustViz-style nearest transcriptomic-neighbour DE plus one-vs-rest marker quality.",
    "Replicate validation: same mean-expression direction in both biological replicates.",
    "Minimum cluster size and worst-cluster Jaccard are warnings, not automatic rejection criteria.",
    "Inspect resolution_selection_summary.tsv, cluster_stability.tsv, and UMAP_resolution_grid.pdf before annotation."),
  file.path(VAL, "resolution_recommendation.txt")
)

# Diagnostic plots.
p_stable_cells <- ggplot(resolution_summary,
                         aes(resolution, stable_cell_fraction_075)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = MIN_STABLE_CELL_FRAC, linetype = 3) +
  geom_vline(xintercept = SELECTED_RES, linetype = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Cells in stable clusters", x = "Resolution",
       y = "Fraction of cells in Jaccard >= 0.75 clusters") + theme_classic()

p_jac <- ggplot(resolution_summary,
                aes(resolution, median_cluster_jaccard)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = MIN_MEDIAN_JACCARD, linetype = 3) +
  geom_vline(xintercept = SELECTED_RES, linetype = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Median cluster stability", x = "Resolution",
       y = "Median Jaccard") + theme_classic()

p_local <- ggplot(resolution_summary,
                  aes(resolution, frac_clusters_locally_separable)) +
  geom_line(na.rm = TRUE) + geom_point(na.rm = TRUE) +
  geom_hline(yintercept = MIN_LOCAL_SEP_FRAC, linetype = 3) +
  geom_vline(xintercept = SELECTED_RES, linetype = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Local transcriptional separation", x = "Resolution",
       y = "Fraction of clusters with >=5 local DE genes") + theme_classic()

p_mark <- ggplot(resolution_summary,
                 aes(resolution, frac_clusters_marker_rich)) +
  geom_line(na.rm = TRUE) + geom_point(na.rm = TRUE) +
  geom_hline(yintercept = MIN_MARKER_RICH_FRAC, linetype = 3) +
  geom_vline(xintercept = SELECTED_RES, linetype = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Cluster marker quality", x = "Resolution",
       y = "Fraction clusters with >=5 strong markers") + theme_classic()

p_rep <- ggplot(resolution_summary,
                aes(resolution, frac_clusters_rep_supported)) +
  geom_line(na.rm = TRUE) + geom_point(na.rm = TRUE) +
  geom_hline(yintercept = MIN_REP_SUPPORTED_FRAC, linetype = 3) +
  geom_vline(xintercept = SELECTED_RES, linetype = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Replicate support for local splits", x = "Resolution",
       y = "Fraction of clusters supported in both replicates") + theme_classic()

p_sil <- ggplot(resolution_summary, aes(resolution, balanced_silhouette)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = SELECTED_RES, linetype = 2) +
  labs(title = "Silhouette (diagnostic only)", x = "Resolution",
       y = "Cluster-balanced mean silhouette") + theme_classic()

ggsave(file.path(FIG, "resolution_selection.pdf"),
       (p_stable_cells | p_jac) / (p_local | p_mark) / (p_rep | p_sil),
       width = 12, height = 14)

# -----------------------------------------------------------------------------
# 8. Apply recommended clustering and inspect sample representation
# -----------------------------------------------------------------------------
selected_col <- res_col(SELECTED_RES)
combined$seurat_clusters <- as.character(combined[[selected_col]][, 1])
Idents(combined) <- "seurat_clusters"

p_clusters <- DimPlot(combined, reduction = "umap.harmony",
                      label = TRUE, repel = TRUE)
p_samples <- DimPlot(combined, reduction = "umap.harmony",
                     group.by = "replicate")
p_split <- DimPlot(combined, reduction = "umap.harmony",
                   split.by = "replicate", group.by = "seurat_clusters",
                   label = TRUE, repel = TRUE)

ggsave(file.path(FIG, "UMAP_selected_clusters.pdf"), p_clusters,
       width = 6, height = 5)
ggsave(file.path(FIG, "UMAP_Harmony_samples.pdf"), p_samples,
       width = 6, height = 5)
ggsave(file.path(FIG, "UMAP_clusters_by_sample.pdf"), p_split,
       width = 11, height = 5)
ggsave(file.path(FIG, "preliminary_clustering_overview.pdf"),
       (p_clusters + labs(title = "Selected clusters",
                          x = "UMAP 1", y = "UMAP 2") |
          p_samples + labs(title = "Biological replicates",
                           x = "UMAP 1", y = "UMAP 2")) +
         patchwork::plot_annotation(tag_levels = "A"),
       width = 10.5, height = 5.2)

cluster_counts <- combined[[]] |>
  dplyr::count(replicate, seurat_clusters, name = "n_cells") |>
  dplyr::group_by(replicate) |>
  dplyr::mutate(percent = 100 * n_cells / sum(n_cells)) |>
  dplyr::ungroup()
write.table(cluster_counts, file.path(OUT, "cells_per_cluster_per_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------------------------------------------------------
# 9. Final marker discovery
# -----------------------------------------------------------------------------
# Genes are ranked independently within each cluster. Replicate consistency is
# reported separately and is not used to replace weak marker evidence.
Idents(combined) <- "seurat_clusters"

markers_raw <- FindAllMarkers(
  combined, assay = "RNA", slot = "data", only.pos = TRUE,
  test.use = "wilcox", min.pct = 0.05, logfc.threshold = 0.10,
  random.seed = SEED, verbose = FALSE
) |>
  dplyr::mutate(
    cluster = as.character(cluster),
    pct_diff = pct.1 - pct.2,
    marker_class = dplyr::case_when(
      p_val_adj < 0.01 & avg_log2FC >= 0.50 & pct.1 >= 0.20 &
        pct_diff >= 0.15 ~ "strong",
      p_val_adj < 0.01 & avg_log2FC >= 0.25 & pct.1 >= 0.15 &
        pct_diff >= 0.10 ~ "standard",
      TRUE ~ "weak"
    ),
    class_rank = dplyr::recode(marker_class,
                               strong = 1L, standard = 2L, weak = 3L),
    specificity_score = avg_log2FC * pmax(pct_diff, 0)
  )

# Calculate replicate direction only for marker-quality candidates, not for every
# weak DE gene. This is much faster and avoids unnecessary large matrix subsets.
marker_candidates <- markers_raw |>
  dplyr::filter(marker_class %in% c("strong", "standard"))

marker_rep_support <- function(marker_tbl) {
  out <- list()
  q <- 1L
  for (cl in sort(unique(marker_tbl$cluster))) {
    genes <- unique(marker_tbl$gene[marker_tbl$cluster == cl])
    genes <- intersect(genes, rownames(rna_data))
    if (!length(genes)) next
    for (rp in sort(unique(replicate_vec))) {
      in_cells <- Cells(combined)[combined$seurat_clusters == cl &
                                    combined$replicate == rp]
      out_cells <- Cells(combined)[combined$seurat_clusters != cl &
                                     combined$replicate == rp]
      if (!length(in_cells) || !length(out_cells)) next
      md <- Matrix::rowMeans(rna_data[genes, in_cells, drop = FALSE]) -
        Matrix::rowMeans(rna_data[genes, out_cells, drop = FALSE])
      out[[q]] <- tibble::tibble(gene = genes, cluster = cl,
                                 replicate = rp, mean_expr_diff = as.numeric(md))
      q <- q + 1L
    }
  }
  dplyr::bind_rows(out)
}

rep_dir_long <- marker_rep_support(marker_candidates)
if (nrow(rep_dir_long)) {
  rep_dir_wide <- rep_dir_long |>
    tidyr::pivot_wider(names_from = replicate, values_from = mean_expr_diff,
                       names_prefix = "mean_diff_")
  markers <- markers_raw |>
    dplyr::left_join(rep_dir_wide, by = c("gene", "cluster"))
} else {
  markers <- markers_raw
}

rep_cols <- grep("^mean_diff_", colnames(markers), value = TRUE)
if (length(rep_cols) >= 2) {
  markers$replicate_consistent <- apply(markers[, rep_cols, drop = FALSE], 1,
                                        function(x) all(is.finite(x) & x > 0))
} else {
  markers$replicate_consistent <- NA
}

anno <- read.delim(ANNO_FILE, header = TRUE, quote = "", check.names = FALSE)
colnames(anno)[1] <- "gene"
anno <- dplyr::distinct(anno, gene, .keep_all = TRUE)
markers_anno <- markers |>
  dplyr::left_join(anno, by = "gene") |>
  dplyr::arrange(cluster, class_rank, dplyr::desc(replicate_consistent),
                 dplyr::desc(specificity_score), dplyr::desc(avg_log2FC))

# Marker table for interpretation: exclude weak/background DE, but retain both
# strong and standard markers. No cross-cluster uniqueness constraint.
marker_genes <- markers_anno |>
  dplyr::filter(marker_class %in% c("strong", "standard"))

write.table(markers_anno, file.path(MRK, "all_markers_annotated.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(marker_genes, file.path(MRK, "marker_genes_filtered.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write_xlsx(list(all_markers = markers_anno,
                filtered_markers = marker_genes),
           file.path(MRK, "marker_genes_tables.xlsx"))

marker_strength <- tibble::tibble(cluster = sort(unique(combined$seurat_clusters))) |>
  dplyr::left_join(
    markers |>
      dplyr::group_by(cluster) |>
      dplyr::summarise(
        n_strong = sum(marker_class == "strong"),
        n_standard = sum(marker_class == "standard"),
        n_rep_consistent = sum(replicate_consistent %in% TRUE, na.rm = TRUE),
        .groups = "drop"
      ), by = "cluster"
  ) |>
  dplyr::mutate(dplyr::across(where(is.numeric), ~tidyr::replace_na(.x, 0)))
write.table(marker_strength, file.path(MRK, "marker_strength_by_cluster.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------------------------------------------------------
# 10. Top markers for heatmap/dotplot and GO-enrichment exports
# -----------------------------------------------------------------------------
select_top_markers <- function(tbl, n) {
  tbl |>
    dplyr::arrange(cluster, class_rank,
                   dplyr::desc(replicate_consistent),
                   dplyr::desc(specificity_score),
                   dplyr::desc(avg_log2FC)) |>
    dplyr::group_by(cluster) |>
    dplyr::slice_head(n = n) |>
    dplyr::mutate(rank_in_cluster = dplyr::row_number()) |>
    dplyr::ungroup()
}

top10 <- select_top_markers(marker_genes, 10)
top5  <- select_top_markers(marker_genes, 5)

write.table(top10, file.path(MRK, "Top10_markers_for_heatmap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(top5, file.path(MRK, "Top5_markers_for_dotplot.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write_xlsx(list(Top10_heatmap = top10, Top5_dotplot = top5),
           file.path(MRK, "top_marker_tables.xlsx"))

# A gene can genuinely be a marker of more than one related cluster. Keep that
# information in the tables; only de-duplicate the feature vector used by Seurat
# plotting functions, which cannot usefully plot the same row repeatedly.
heatmap_features <- top10$gene[!duplicated(top10$gene)]
dotplot_features <- top5$gene[!duplicated(top5$gene)]

if (length(heatmap_features)) {
  combined <- ScaleData(combined, features = heatmap_features, verbose = FALSE)
  # A balanced cell subset makes the heatmap readable and prevents the largest
  # cluster from dominating the image. Marker statistics still use all cells.
  set.seed(SEED)
  heatmap_cells <- unlist(lapply(split(Cells(combined), Idents(combined)), function(x) {
    sample(x, min(length(x), 200))
  }), use.names = FALSE)

  p_heat <- DoHeatmap(combined, features = heatmap_features,
                      cells = heatmap_cells, group.by = "seurat_clusters",
                      raster = TRUE) + NoLegend()
  ggsave(file.path(FIG, "heatmap_Top10_markers.pdf"), p_heat,
         width = 13, height = max(8, 0.10 * length(heatmap_features) + 3),
         limitsize = FALSE)
}

if (length(dotplot_features)) {
  p_dot <- DotPlot(combined, features = dotplot_features,
                   group.by = "seurat_clusters") +
    RotatedAxis() + labs(x = "Marker gene", y = "Cluster")
  ggsave(file.path(FIG, "dotplot_Top5_markers.pdf"), p_dot,
         width = max(12, 0.28 * length(dotplot_features) + 5),
         height = 7, limitsize = FALSE)
}

# Export permissive and replicate-consistent marker lists. Expressed genes form
# the enrichment universe.
go_markers_all <- marker_genes |>
  dplyr::select(cluster, gene, p_val_adj, avg_log2FC, pct.1, pct.2, pct_diff,
                marker_class, replicate_consistent, dplyr::everything())

go_markers_rep <- go_markers_all |>
  dplyr::filter(replicate_consistent %in% TRUE)

write.table(go_markers_all, file.path(MRK, "GO_marker_genes_by_cluster.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(go_markers_rep,
            file.path(MRK, "GO_marker_genes_rep_consistent_by_cluster.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

GO_DIR <- file.path(MRK, "GO_gene_lists")
dir.create(GO_DIR, recursive = TRUE, showWarnings = FALSE)
for (cl in sort(unique(go_markers_all$cluster))) {
  genes_cl <- unique(go_markers_all$gene[go_markers_all$cluster == cl])
  writeLines(genes_cl, file.path(GO_DIR, paste0("cluster_", cl, "_markers.txt")))

  genes_rep <- unique(go_markers_rep$gene[go_markers_rep$cluster == cl])
  writeLines(genes_rep,
             file.path(GO_DIR, paste0("cluster_", cl, "_markers_rep_consistent.txt")))
}

writeLines(rownames(combined), file.path(MRK, "GO_background_expressed_genes.txt"))

# One workbook with one sheet per cluster for convenient manual inspection/GO.
go_sheets <- split(go_markers_all, go_markers_all$cluster)
names(go_sheets) <- paste0("cluster_", names(go_sheets))
if (length(go_sheets)) {
  write_xlsx(go_sheets, file.path(MRK, "GO_markers_by_cluster.xlsx"))
}

# -----------------------------------------------------------------------------
# 11. Transfer old Drosophila marker and cell-cycle annotations to CPBv5 IDs
# -----------------------------------------------------------------------------
# The Drosophila marker and cell-cycle tables were curated against the old
# LdNA identifiers. Transfer them through table_LDECv5_vs_LdNA.txt before
# interpreting the resulting Seurat object.
extract_ldna_ids <- function(x) {
  x <- ifelse(is.na(x), "", x)
  stringr::str_extract_all(x, "LdNA_[0-9]+")
}

crosswalk_raw <- read.delim(CROSSWALK_FILE, header = TRUE, quote = "",
                            check.names = FALSE)
crosswalk <- dplyr::bind_rows(lapply(seq_len(nrow(crosswalk_raw)), function(i) {
  old <- unique(unlist(extract_ldna_ids(crosswalk_raw$LdNA_IDs[i])))
  if (!length(old)) return(NULL)
  tibble::tibble(
    old_gene = old,
    gene = crosswalk_raw$Geneid[i],
    crosswalk_annotation = crosswalk_raw$blast_annotation[i]
  )
})) |>
  dplyr::distinct(old_gene, gene, .keep_all = TRUE)

expand_old_id_table <- function(tbl, columns) {
  dplyr::bind_rows(lapply(seq_len(nrow(tbl)), function(i) {
    old <- unique(unlist(extract_ldna_ids(unlist(tbl[i, columns], use.names = FALSE))))
    if (!length(old)) old <- NA_character_
    dplyr::bind_cols(tbl[i, , drop = FALSE],
                     tibble::tibble(old_gene = old))
  }))
}

summarise_gene_expression <- function(genes) {
  genes <- intersect(unique(genes), rownames(rna_data))
  if (!length(genes)) return(tibble::tibble())
  rna_counts <- GetAssayData(combined, assay = "RNA", layer = "counts")
  clusters <- sort(unique(as.character(combined$seurat_clusters)))
  dplyr::bind_rows(lapply(genes, function(g) {
    dplyr::bind_rows(lapply(clusters, function(cl) {
      cells <- Cells(combined)[combined$seurat_clusters == cl]
      tibble::tibble(
        gene = g,
        cluster = cl,
        avg_log_normalized_expr = mean(as.numeric(rna_data[g, cells])),
        pct_detected = mean(as.numeric(rna_counts[g, cells]) > 0) * 100
      )
    }))
  }))
}

droso_markers_raw <- read.delim(DROSO_MARKER_FILE, header = TRUE, quote = "",
                                check.names = FALSE)
droso_markers_mapped <- expand_old_id_table(
  droso_markers_raw, c("CPB Orthofinder", "CPB gene")
) |>
  dplyr::left_join(crosswalk, by = "old_gene") |>
  dplyr::mutate(
    mapping_status = dplyr::case_when(
      is.na(old_gene) ~ "no_old_cpb_id_in_marker_row",
      is.na(gene) ~ "old_id_not_found_in_v5_crosswalk",
      !(gene %in% rownames(combined)) ~ "mapped_but_not_detected_after_filtering",
      TRUE ~ "mapped_and_detected"
    )
  )

droso_marker_detected <- droso_markers_mapped |>
  dplyr::filter(mapping_status == "mapped_and_detected") |>
  dplyr::distinct(`Gene name`, `Annotation symbol`, `Flybase ID`,
                  `marker of`, ref, old_gene, gene)
droso_marker_expr_summary <- summarise_gene_expression(droso_marker_detected$gene)
droso_marker_expr <- droso_marker_detected |>
  dplyr::left_join(droso_marker_expr_summary, by = "gene",
                   relationship = "many-to-many") |>
  dplyr::group_by(`Gene name`, gene) |>
  dplyr::mutate(max_cluster = cluster[which.max(avg_log_normalized_expr)],
                max_pct_detected = max(pct_detected)) |>
  dplyr::ungroup()

write.table(droso_markers_mapped,
            file.path(ANN, "drosophila_marker_old_to_v5_mapping.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(droso_marker_expr,
            file.path(ANN, "drosophila_marker_expression_by_cluster.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cell_cycle_raw <- read.delim(CELL_CYCLE_FILE, header = TRUE, quote = "",
                             check.names = FALSE)
cell_cycle_mapped <- expand_old_id_table(cell_cycle_raw, "Geneid") |>
  dplyr::left_join(crosswalk, by = "old_gene") |>
  dplyr::mutate(
    mapping_status = dplyr::case_when(
      is.na(old_gene) ~ "no_old_cpb_id_in_cell_cycle_row",
      is.na(gene) ~ "old_id_not_found_in_v5_crosswalk",
      !(gene %in% rownames(combined)) ~ "mapped_but_not_detected_after_filtering",
      TRUE ~ "mapped_and_detected"
    )
  )

write.table(cell_cycle_mapped,
            file.path(ANN, "cell_cycle_old_to_v5_mapping.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cell_cycle_detected <- cell_cycle_mapped |>
  dplyr::filter(mapping_status == "mapped_and_detected") |>
  dplyr::distinct(phase, dm_orth, old_gene, gene)
cell_cycle_expr_summary <- summarise_gene_expression(cell_cycle_detected$gene)
cell_cycle_expr <- cell_cycle_detected |>
  dplyr::left_join(cell_cycle_expr_summary, by = "gene")
write.table(cell_cycle_expr,
            file.path(ANN, "cell_cycle_gene_expression_by_cluster.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

for (ph in sort(unique(cell_cycle_mapped$phase))) {
  genes_ph <- cell_cycle_mapped |>
    dplyr::filter(phase == ph, mapping_status == "mapped_and_detected") |>
    dplyr::pull(gene) |>
    unique() |>
    intersect(rownames(rna_data))
  score_name <- paste0("old_", gsub("[^A-Za-z0-9]+", "_", tolower(ph)), "_score")
  combined[[score_name]] <- if (length(genes_ph)) {
    Matrix::colMeans(rna_data[genes_ph, Cells(combined), drop = FALSE])
  } else {
    NA_real_
  }
}

cell_cycle_score_cols <- grep("^old_.*_score$", colnames(combined[[]]),
                              value = TRUE)
cell_cycle_scores <- combined[[]] |>
  dplyr::select(seurat_clusters, dplyr::all_of(cell_cycle_score_cols)) |>
  tidyr::pivot_longer(-seurat_clusters, names_to = "score",
                      values_to = "value") |>
  dplyr::group_by(seurat_clusters, score) |>
  dplyr::summarise(mean_score = mean(value, na.rm = TRUE),
                   median_score = median(value, na.rm = TRUE),
                   .groups = "drop")
write.table(cell_cycle_scores,
            file.path(ANN, "cell_cycle_scores_by_cluster.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cluster_annotation_support <- marker_strength |>
  dplyr::left_join(nearest_local |>
                     dplyr::filter(resolution == SELECTED_RES) |>
                     dplyr::select(cluster, nearest_neighbor = neighbor,
                                   n_meaningful_DE, n_rep_supported),
                   by = "cluster") |>
  dplyr::left_join(cluster_counts |>
                     dplyr::group_by(seurat_clusters) |>
                     dplyr::summarise(n_cells = sum(n_cells),
                                      .groups = "drop") |>
                     dplyr::rename(cluster = seurat_clusters),
                   by = "cluster") |>
  dplyr::mutate(
    annotation_confidence = dplyr::case_when(
      n_strong == 0 & n_meaningful_DE == 0 ~ "low_no_independent_marker_support",
      n_strong >= 5 & n_rep_consistent >= 5 ~ "high_marker_and_replicate_support",
      n_strong >= 5 ~ "moderate_marker_support",
      TRUE ~ "manual_review_required"
    )
  )
write.table(cluster_annotation_support,
            file.path(ANN, "cluster_annotation_support.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Marker-transfer and cell-cycle diagnostics.
if (nrow(droso_marker_expr)) {
  droso_features <- unique(droso_marker_expr$gene)
  p_droso <- DotPlot(combined, features = droso_features,
                     group.by = "seurat_clusters") +
    coord_flip() +
    labs(x = "Cluster", y = "Mapped Drosophila marker ortholog",
         title = "Old Drosophila marker table transferred to CPBv5 IDs") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 7),
          plot.margin = margin(5.5, 40, 5.5, 5.5))
  ggsave(file.path(FIG, "drosophila_marker_transfer_dotplot.pdf"), p_droso,
         width = 9, height = max(5, 0.22 * length(droso_features) + 2),
         limitsize = FALSE)
}

if (length(cell_cycle_score_cols)) {
  cc_cell_long <- combined[[]] |>
    dplyr::select(seurat_clusters, dplyr::all_of(cell_cycle_score_cols)) |>
    tidyr::pivot_longer(-seurat_clusters, names_to = "score",
                        values_to = "value")
  p_cc <- ggplot(cc_cell_long, aes(seurat_clusters, value,
                                   fill = seurat_clusters)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
    facet_wrap(~score, scales = "free_y") +
    labs(x = "Cluster", y = "Mean log-normalized expression",
         title = "Cell-cycle scores from old-ID Drosophila ortholog table") +
    theme_classic() + theme(legend.position = "none")
  ggsave(file.path(FIG, "cell_cycle_scores_by_cluster.pdf"), p_cc,
         width = 8, height = 4.8)
}

# Exploratory marker overview for the selected broad partition.
if (length(heatmap_features) && length(dotplot_features)) {
  p_heat <- DoHeatmap(combined, features = heatmap_features,
                         cells = heatmap_cells, group.by = "seurat_clusters",
                         raster = TRUE) +
    NoLegend() +
    theme(axis.text.y = element_text(size = 6),
          plot.margin = margin(5.5, 20, 5.5, 5.5))

  p_dot <- DotPlot(combined, features = dotplot_features,
                      group.by = "seurat_clusters", dot.scale = 5) +
    coord_flip() +
    labs(x = "Top marker gene", y = "Cluster") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 7),
          axis.text.x = element_text(size = 9),
          plot.margin = margin(5.5, 40, 5.5, 5.5))

  marker_overview <- (p_heat / p_dot) +
    patchwork::plot_annotation(tag_levels = "A")
  ggsave(file.path(FIG, "broad_cluster_marker_overview.pdf"), marker_overview,
         width = 12, height = 13, limitsize = FALSE)
}

saveRDS(combined, file.path(OUT, "combined_selected_resolution.rds"))
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))

message("\nRecommended resolution: ", SELECTED_RES)
message("Selection status: ", selection_status)
message("Resolution table: ", file.path(VAL, "resolution_selection_summary.tsv"))
message("Marker table:     ", file.path(MRK, "marker_genes_tables.xlsx"))
message("GO gene lists:    ", GO_DIR)
