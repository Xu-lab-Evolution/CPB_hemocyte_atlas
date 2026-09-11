# Export the CPB hemocyte atlas as a ShinyCell2 web portal.

suppressPackageStartupMessages({
  library(Seurat)
  library(ShinyCell2)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(hdf5r)
})

source(file.path("R", "config.R"))

INPUT_RDS <- file.path(OUT, "combined_annotated.rds")
APP_DIR <- file.path(OUT, "shinycell2_cpb")
PREFIX <- "cpb_"
SHINYCELL2_SOURCE <- "https://github.com/the-ouyang-lab/ShinyCell2"
SHINYCELL2_COMMIT <- "33bfc8ba232f0c829b6b23181cb83089d58e7879"

dir.create(APP_DIR, recursive = TRUE, showWarnings = FALSE)

obj <- JoinLayers(readRDS(INPUT_RDS))
DefaultAssay(obj) <- "RNA"

required_metadata <- c(
  "replicate", "res_0.2", "cluster_annotation", "short_label",
  "hemocyte_class", "hemocyte_state",
  "annotation_confidence", "cell_cycle_signal",
  "cell_cycle_s_score", "cell_cycle_g2m_score",
  "nCount_RNA", "nFeature_RNA", "percent.mt", "scDblFinder.score"
)
stopifnot(all(required_metadata %in% colnames(obj@meta.data)))
stopifnot(all(c("umap.harmony", "umap.unintegrated") %in% Reductions(obj)))
stopifnot(all(c("data", "counts") %in% Layers(obj[["RNA"]])))

cluster_levels <- as.character(1:8)
annotation_levels <- c(
  "Proliferating plasmatocytes",
  "Immune-active granulocytes",
  "Stress-responsive oenocytoids",
  "Adhesive/clotting plasmatocytes",
  "Prohemocytes (provisional)",
  "Differentiating oenocytoids",
  "Activated oenocytoids",
  "Germline/meiotic contaminant"
)
state_levels <- c(
  "Proliferating", "Immune-active", "Stress-responsive",
  "Adhesive/clotting", "RNA-low/provisional", "Differentiating",
  "Activated antimicrobial/stress", "Germline/meiotic contamination"
)

cluster_id <- as.character(obj$res_0.2)
portal_meta <- data.frame(
  sample = factor(
    obj$replicate,
    levels = c("rep1", "rep2"), labels = c("Sample 1", "Sample 2")
  ),
  cluster = factor(
    paste0("Cluster ", cluster_id), levels = paste0("Cluster ", cluster_levels)
  ),
  cell_annotation = factor(obj$cluster_annotation, levels = annotation_levels),
  hemocyte_class = factor(
    obj$hemocyte_class,
    levels = c(
      "Plasmatocyte", "Granulocyte", "Oenocytoid", "Prohemocyte",
      "Non-hemocyte"
    )
  ),
  hemocyte_state = factor(obj$hemocyte_state, levels = state_levels),
  annotation_confidence = factor(
    obj$annotation_confidence,
    levels = c(
      "high", "moderate-to-high", "moderate",
      "moderate from morphology; low from RNA"
    )
  ),
  analysis_status = factor(
    ifelse(
      cluster_id == "8",
      "Collection contaminant; excluded from hemocyte inference",
      "Hemocyte analysis set"
    ),
    levels = c(
      "Hemocyte analysis set",
      "Collection contaminant; excluded from hemocyte inference"
    )
  ),
  independent_recovery = factor(
    ifelse(
      cluster_id == "8", "Not independently recovered",
      "Independently recovered in both samples"
    ),
    levels = c(
      "Independently recovered in both samples",
      "Not independently recovered"
    )
  ),
  cell_cycle_state = factor(
    obj$cell_cycle_signal,
    levels = c(
      "Below upper-quartile reference", "S high only", "G2/M high only",
      "S and G2/M high"
    )
  ),
  s_phase_score = as.numeric(obj$cell_cycle_s_score),
  g2m_phase_score = as.numeric(obj$cell_cycle_g2m_score),
  umi_count = as.numeric(obj$nCount_RNA),
  detected_genes = as.numeric(obj$nFeature_RNA),
  mitochondrial_percent = as.numeric(obj$percent.mt),
  doublet_score = as.numeric(obj$scDblFinder.score),
  row.names = colnames(obj),
  check.names = FALSE
)
stopifnot(!anyNA(portal_meta$cell_annotation))

# Keep only curated, interpretable portal metadata in the export copy.
obj@meta.data <- portal_meta
Idents(obj) <- obj$cell_annotation

metadata_order <- colnames(portal_meta)
sc_conf <- createConfig(obj, meta.to.include = metadata_order, maxLevels = 50)
sc_conf <- modMetaName(
  sc_conf,
  meta.to.mod = metadata_order,
  new.name = c(
    "Sample", "Cluster", "Cell annotation", "Hemocyte class",
    "Hemocyte state", "Annotation confidence", "Analysis status",
    "Replicate recovery", "Cell-cycle state", "S-phase score",
    "G2/M-phase score", "UMI count", "Detected genes",
    "Mitochondrial reads (%)", "Doublet score"
  )
)
sc_conf <- reorderMeta(sc_conf, metadata_order)
sc_conf <- modDefault(
  sc_conf, default1 = "cell_annotation", default2 = "sample"
)

cluster_colours <- c(
  "Cluster 1" = "#E64B35", "Cluster 2" = "#4DBBD5",
  "Cluster 3" = "#00A087", "Cluster 4" = "#3C5488",
  "Cluster 5" = "#F39B7F", "Cluster 6" = "#91D1C2",
  "Cluster 7" = "#008B8B", "Cluster 8" = "#7F7F7F"
)
annotation_colours <- setNames(unname(cluster_colours), annotation_levels)
state_colours <- setNames(unname(cluster_colours), state_levels)
class_colours <- c(
  "Plasmatocyte" = "#E64B35", "Granulocyte" = "#4DBBD5",
  "Oenocytoid" = "#00A087", "Prohemocyte" = "#F39B7F",
  "Non-hemocyte" = "#7F7F7F"
)

apply_named_palette <- function(config, metadata, palette) {
  config_levels <- strsplit(
    config[config$ID == metadata, ]$fID, "\\|"
  )[[1]]
  stopifnot(all(config_levels %in% names(palette)))
  modColours(config, metadata, unname(palette[config_levels]))
}

sc_conf <- apply_named_palette(sc_conf, "cluster", cluster_colours)
sc_conf <- apply_named_palette(
  sc_conf, "cell_annotation", annotation_colours
)
sc_conf <- apply_named_palette(sc_conf, "hemocyte_class", class_colours)
sc_conf <- apply_named_palette(sc_conf, "hemocyte_state", state_colours)
sc_conf <- apply_named_palette(
  sc_conf, "sample",
  c("Sample 1" = "#3C5488", "Sample 2" = "#F39B7F")
)
sc_conf <- apply_named_palette(
  sc_conf, "analysis_status",
  c(
    "Hemocyte analysis set" = "#009E73",
    "Collection contaminant; excluded from hemocyte inference" = "#7F7F7F"
  )
)
sc_conf <- apply_named_palette(
  sc_conf, "independent_recovery",
  c(
    "Independently recovered in both samples" = "#009E73",
    "Not independently recovered" = "#7F7F7F"
  )
)
sc_conf <- apply_named_palette(
  sc_conf, "cell_cycle_state",
  c(
    "Below upper-quartile reference" = "#BDBDBD",
    "S high only" = "#56B4E9", "G2/M high only" = "#E69F00",
    "S and G2/M high" = "#CC79A7"
  )
)

default_gene1 <- "LDECv5g03007"
default_gene2 <- "LDECv5g15036"
default_multigene <- c(
  "LDECv5g04564", "LDECv5g10569", "LDECv5g13695",
  "LDECv5g03007", "LDECv5g08078", "LDECv5g15036",
  "LDECv5g00509", "LDECv5g08619", "LDECv5g06984",
  "LDECv5g00806"
)
stopifnot(all(c(default_gene1, default_gene2, default_multigene) %in% rownames(obj)))

makeShinyFiles(
  obj, sc_conf,
  assay = "RNA", assay.slot = "data",
  dimred.to.use = c("umap.harmony", "umap.unintegrated"),
  shiny.prefix = PREFIX, shiny.dir = APP_DIR,
  default.gene1 = default_gene1,
  default.gene2 = default_gene2,
  default.multigene = default_multigene,
  default.dimred = "umap.harmony",
  chunkSize = 250
)

makeShinyCodes(
  shiny.title = "Colorado potato beetle hemocyte single-cell atlas",
  shiny.footnotes = paste(
    "Leptinotarsa decemlineata hemocyte scRNA-seq atlas.",
    "Cluster 8 is retained as a collection-contamination control and excluded",
    "from hemocyte inference."
  ),
  shiny.prefix = PREFIX,
  shiny.headers = "CPB hemocytes",
  shiny.dir = APP_DIR,
  defPtSiz = 0.65
)

# ShinyCell2 updates the assay selector before the feature selector. Add guards
# so this brief transition waits silently instead of indexing an unavailable
# feature and emitting a transient HDF5 error.
insert_function_guard <- function(lines, function_name, guard_lines) {
  function_start <- grep(
    paste0("^", function_name, " <- function"), lines
  )
  stopifnot(length(function_start) == 1L)
  opening_line <- function_start - 1L + which(
    grepl("\\{[[:space:]]*$", lines[function_start:length(lines)])
  )[1]
  append(lines, guard_lines, after = opening_line)
}

helper_file <- file.path(APP_DIR, "shinyFunc.R")
helper_lines <- readLines(helper_file, warn = FALSE)
helper_lines <- insert_function_guard(
  helper_lines, "sc2Ddimr",
  c(
    "  req(inpdr, inp1, inpDtyp, inpsiz, inpord, inpcol, inpfsz, inpasp)",
    "  requestedAssay = gsub(\"^Assay: \", \"\", inpDtyp)",
    "  if(requestedAssay != \"Cell Information\"){",
    "    req(requestedAssay %in% names(inpGene))",
    "    req(inp1 %in% names(inpGene[[requestedAssay]]))",
    "  }"
  )
)
helper_lines <- insert_function_guard(
  helper_lines, "sc2Dnum",
  c(
    "  req(inpdr, inp1, inpDtyp, inpsplt)",
    "  requestedAssay = gsub(\"^Assay: \", \"\", inpDtyp)",
    "  if(requestedAssay != \"Cell Information\"){",
    "    req(requestedAssay %in% names(inpGene))",
    "    req(inp1 %in% names(inpGene[[requestedAssay]]))",
    "  }"
  )
)
writeLines(helper_lines, helper_file)

# Add searchable functional names without changing HDF5 row indexing.
functional_annotation <- read.delim(
  input_file("LdecV5_functional_annotation.txt"),
  check.names = FALSE, quote = "", fill = TRUE
) |>
  transmute(
    gene_id = Geneid,
    gene_name = case_when(
      !is.na(Uniprot_annotation) & Uniprot_annotation != "NA" &
        nzchar(Uniprot_annotation) ~ str_trim(str_split_fixed(
          Uniprot_annotation, ";", 2
        )[, 1]),
      TRUE ~ NA_character_
    )
  )

curated_aliases <- bind_rows(
  read.delim(
    file.path(OUT, "final_annotation", "selected_classification_markers.tsv"),
    check.names = FALSE
  ) |>
    transmute(gene_id = gene, gene_name = display_name, alias_priority = 1L),
  read.delim(
    file.path(OUT, "cell_cycle", "selected_cell_cycle_genes.tsv"),
    check.names = FALSE
  ) |>
    transmute(gene_id = gene, gene_name = display_name, alias_priority = 2L)
) |>
  filter(!is.na(gene_id), !is.na(gene_name), nzchar(gene_name)) |>
  arrange(alias_priority) |>
  distinct(gene_id, .keep_all = TRUE) |>
  select(gene_id, curated_gene_name = gene_name)

drosophila_orthologs <- read.delim(
  input_file("one-to-one-orthologs_CPB_Dm.tsv"),
  check.names = FALSE
) |>
  transmute(
    gene_id = Ldecemlineata,
    drosophila_gene_id = Dmelanogaster
  ) |>
  filter(!is.na(gene_id), !is.na(drosophila_gene_id)) |>
  distinct(gene_id, .keep_all = TRUE)

tribolium_orthologs <- read.delim(
  input_file("one-to-one-orthologs_CPB_Tc.tsv"),
  check.names = FALSE
) |>
  transmute(
    gene_id = Ldecemlineata,
    tribolium_gene_id = Tcastaneum
  ) |>
  filter(!is.na(gene_id), !is.na(tribolium_gene_id)) |>
  distinct(gene_id, .keep_all = TRUE)

make_portal_label <- function(gene_name, gene_id, drosophila_id, tribolium_id) {
  parts <- c(
    if (!is.na(gene_name) && nzchar(gene_name)) gene_name,
    gene_id,
    if (!is.na(drosophila_id) && nzchar(drosophila_id)) {
      paste0("Dm:", drosophila_id)
    },
    if (!is.na(tribolium_id) && nzchar(tribolium_id)) {
      paste0("Tc:", tribolium_id)
    }
  )
  paste(parts, collapse = " | ")
}

feature_aliases <- tibble(gene_id = rownames(obj)) |>
  left_join(functional_annotation, by = "gene_id") |>
  left_join(curated_aliases, by = "gene_id") |>
  left_join(drosophila_orthologs, by = "gene_id") |>
  left_join(tribolium_orthologs, by = "gene_id") |>
  mutate(
    gene_name = coalesce(curated_gene_name, gene_name),
    portal_label = mapply(
      make_portal_label,
      gene_name, gene_id, drosophila_gene_id, tribolium_gene_id,
      USE.NAMES = FALSE
    )
  ) |>
  select(
    gene_id, gene_name, drosophila_gene_id, tribolium_gene_id, portal_label
  )
stopifnot(!anyDuplicated(feature_aliases$portal_label))

alias_lookup <- setNames(feature_aliases$portal_label, feature_aliases$gene_id)
gene_file <- file.path(APP_DIR, paste0(PREFIX, "gene.rds"))
gene_map <- readRDS(gene_file)
for (assay_name in names(gene_map)) {
  original_ids <- names(gene_map[[assay_name]])
  replacement <- unname(alias_lookup[original_ids])
  replacement[is.na(replacement)] <- original_ids[is.na(replacement)]
  stopifnot(!anyDuplicated(replacement))
  names(gene_map[[assay_name]]) <- replacement
}
saveRDS(gene_map, gene_file)

default_file <- file.path(APP_DIR, paste0(PREFIX, "def.rds"))
portal_defaults <- readRDS(default_file)
portal_defaults$gene1$RNA <- unname(alias_lookup[default_gene1])
portal_defaults$gene2$RNA <- unname(alias_lookup[default_gene2])
portal_defaults$genes$RNA <- unname(alias_lookup[default_multigene])
saveRDS(portal_defaults, default_file)

write.table(
  feature_aliases,
  file.path(APP_DIR, "cpb_feature_aliases.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  feature_aliases |>
    filter(!is.na(drosophila_gene_id) | !is.na(tribolium_gene_id)),
  file.path(APP_DIR, "cpb_ortholog_lookup.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  as.data.frame(sc_conf),
  file.path(APP_DIR, "cpb_shinycell2_config.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

manifest <- c(
  paste("generated_at", format(Sys.time(), tz = "Europe/Berlin"), sep = "\t"),
  paste("input_rds", INPUT_RDS, sep = "\t"),
  paste("input_md5", unname(tools::md5sum(INPUT_RDS)), sep = "\t"),
  paste("cells", ncol(obj), sep = "\t"),
  paste("features", nrow(obj), sep = "\t"),
  paste(
    "features_with_Drosophila_ortholog",
    sum(!is.na(feature_aliases$drosophila_gene_id)), sep = "\t"
  ),
  paste(
    "features_with_Tribolium_ortholog",
    sum(!is.na(feature_aliases$tribolium_gene_id)), sep = "\t"
  ),
  paste("assay", "RNA/data", sep = "\t"),
  paste("dimension_reductions", "umap.harmony; umap.unintegrated", sep = "\t"),
  paste("default_dimension_reduction", "umap.harmony", sep = "\t"),
  paste("ShinyCell2_version", as.character(packageVersion("ShinyCell2")), sep = "\t"),
  paste("ShinyCell2_source", SHINYCELL2_SOURCE, sep = "\t"),
  paste("ShinyCell2_commit", SHINYCELL2_COMMIT, sep = "\t")
)
writeLines(manifest, file.path(APP_DIR, "export_manifest.tsv"))

writeLines(
  c(
    "# CPB hemocyte ShinyCell2 portal",
    "",
    "This directory is a self-contained ShinyCell2 application generated from",
    "the final CPB hemocyte Seurat object.",
    "",
    "## Install dependencies",
    "",
    "```r",
    "install.packages(c(",
    "  'shiny', 'shinyhelper', 'data.table', 'Matrix', 'DT', 'magrittr',",
    "  'ggplot2', 'ggrepel', 'hdf5r', 'ggdendro', 'gridExtra', 'ggpubr'",
    "))",
    "remotes::install_github(",
    "  'the-ouyang-lab/ShinyCell2@33bfc8ba232f0c829b6b23181cb83089d58e7879'",
    ")",
    "```",
    "",
    "## Run locally",
    "",
    "```r",
    "shiny::runApp('results/shinycell2_cpb')",
    "```",
    "",
    "The default view uses the Harmony UMAP and final cell annotations. The",
    "unintegrated UMAP remains available for checking sample-associated structure.",
    "Genes can be searched by CPB gene ID, functional name, Drosophila FBgn ID,",
    "or Tribolium TC gene ID. Ortholog searches resolve to the corresponding CPB",
    "feature; expression values are always from the CPB RNA assay.",
    "Cluster 8 is retained as a contamination control and is explicitly marked as",
    "excluded from hemocyte inference.",
    "",
    "## Deploy",
    "",
    "Deploy the complete directory without removing or renaming the HDF5/RDS files.",
    "For shinyapps.io, run `rsconnect::deployApp(appDir =",
    "'results/shinycell2_cpb')` after configuring an account.",
    "",
    "See `export_manifest.tsv` for source and version provenance."
  ),
  file.path(APP_DIR, "README.md")
)

capture.output(sessionInfo(), file = file.path(APP_DIR, "sessionInfo.txt"))

# Validate the generated data contract before returning.
meta_export <- readRDS(file.path(APP_DIR, paste0(PREFIX, "meta.rds")))
dimr_export <- readRDS(file.path(APP_DIR, paste0(PREFIX, "dimr.rds")))
defaults_export <- readRDS(default_file)
h5_file <- file.path(APP_DIR, paste0(PREFIX, "assay_RNA.h5"))
h5 <- H5File$new(h5_file, mode = "r")
h5_dims <- h5[["grp/data"]]$dims

search_examples <- bind_rows(
  feature_aliases |>
    filter(!is.na(drosophila_gene_id)) |>
    slice(1) |>
    transmute(
      species = "Drosophila melanogaster",
      query = drosophila_gene_id,
      gene_id, portal_label
    ),
  feature_aliases |>
    filter(!is.na(tribolium_gene_id)) |>
    slice(1) |>
    transmute(
      species = "Tribolium castaneum",
      query = tribolium_gene_id,
      gene_id, portal_label
    )
) |>
  mutate(
    hdf5_row = unname(gene_map$RNA[portal_label]),
    expressing_cells = vapply(
      hdf5_row,
      function(i) sum(h5[["grp/data"]][i, ] > 0),
      integer(1)
    )
  )
h5$close_all()

write.table(
  search_examples,
  file.path(APP_DIR, "ortholog_search_validation.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

stopifnot(
  nrow(meta_export) == ncol(obj),
  length(gene_map$RNA) == nrow(obj),
  identical(names(dimr_export), c("umap.harmony", "umap.unintegrated")),
  all(vapply(dimr_export, nrow, integer(1)) == ncol(obj)),
  all(as.integer(h5_dims) == c(nrow(obj), ncol(obj))),
  any(grepl("FBgn", names(gene_map$RNA), fixed = TRUE)),
  any(grepl("TC", names(gene_map$RNA), fixed = TRUE)),
  all(search_examples$expressing_cells > 0),
  defaults_export$dimrd[1] == "umap.harmony",
  file.exists(file.path(APP_DIR, "server.R")),
  file.exists(file.path(APP_DIR, "ui.R"))
)

message("ShinyCell2 portal written to: ", APP_DIR)
