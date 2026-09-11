# CPB hemocyte single-cell RNA-seq analysis

This repository contains the analysis used to process, cluster, validate, and
annotate the two CPB hemocyte scRNA-seq samples. The workflow starts from the
filtered Cell Ranger matrices and produces the analyzed Seurat object, marker
and GO-enrichment tables, manuscript figures, and an optional ShinyCell2 portal.

## Required input

Place the following files and directories in the repository root:

```text
CPBv5_single_mod_mito_prot_coding.gtf
LdecV5_functional_annotation.txt
table_LDECv5_vs_LdNA.txt
marker_genes_droso_cpb.txt
cpb_droso_cellcycle.tsv
one-to-one-orthologs_CPB_Dm.tsv
one-to-one-orthologs_CPB_Tc.tsv
```
The matrix directors (filtered_feature_bc_matrix_s1 and filtered_feature_bc_matrix_s2) can be download from:
10.6084/m9.figshare.33620566

The two matrix directories must contain the standard `barcodes.tsv.gz`,
`features.tsv.gz`, and `matrix.mtx.gz` files produced by Cell Ranger. The old
and current gene annotations are linked through `table_LDECv5_vs_LdNA.txt`
before Drosophila markers or cell-cycle genes are evaluated.

## Software

The workflow was written for R 4.x and Seurat 5. Install the required CRAN and
Bioconductor packages with:

```bash
Rscript install_dependencies.R
```

Each analysis stage records `sessionInfo()` in its output directory.

## Run the analysis

From the repository root:

```bash
Rscript run_analysis.R
```

Results are written to `results/`. Existing independent-clustering caches are
reused by stage 4. Delete `results/replicate_consensus_GO/independent_*rds` to
recompute those cached objects after changing clustering parameters.

Input and output locations can be set without editing the scripts:

```bash
CPB_INPUT_DIR=/path/to/input \
CPB_OUTPUT_DIR=/path/to/results \
Rscript run_analysis.R
```

The optional ShinyCell2 export is run separately after the main analysis:

```bash
Rscript R/08_export_shinycell2.R
```

ShinyCell2 is pinned in the export script to commit
`33bfc8ba232f0c829b6b23181cb83089d58e7879`.

## Workflow

1. `R/01_preprocess_and_cluster.R`: quality control, doublet detection,
   normalization, PCA, Harmony integration, graph clustering, resolution
   assessment, broad-cluster markers, and transfer of old annotation IDs.
2. `R/02_cell_cycle.R`: control-gene-adjusted S and G2/M module scores and
   expression summaries for selected mapped cell-cycle orthologs.
3. `R/03_cluster_validation.R`: broad-to-fine cluster hierarchy, replicate
   balance, graph-parameter sensitivity, low-RNA population assessment, and
   targeted collection-contamination screens.
4. `R/04_replicate_consensus_markers_GO.R`: independent clustering of each
   sample, recovery of integrated clusters, replicate-consistent marker
   selection, and GO enrichment.
5. `R/05_hemocyte_annotation.R`: evidence-based assignment of established
   insect hemocyte classes and transcriptional states.
6. `R/06_main_figures.R`: manuscript Figures 1-4.
7. `R/07_GO_and_annotation_markers.R`: top Biological Process enrichment
   results and annotation-marker expression dot plot (Figure 5).
8. `R/08_export_shinycell2.R`: optional ShinyCell2 application export.

The portal gene selectors accept CPB gene IDs, functional names, Drosophila
FlyBase gene IDs, and Tribolium gene IDs. Ortholog queries select the linked
CPB feature and display its expression in the CPB RNA assay.

## Clustering criteria

Candidate resolutions from 0.1 to 1.0 are compared using 80% cell subsampling,
cluster Jaccard stability, local differential-expression separation, marker
strength, and agreement of expression direction between samples. Resolution
0.1 is retained as the broad stable partition. Resolution 0.2 is used only as a
secondary biological partition: each fine cluster must contain at least 15
cells per sample and recover an independently clustered population in both
samples with at least 60% purity and 60% coverage. Fine clusters 1-7 satisfy the
replicate-recovery criterion. Fine cluster 8 contains 30 cells, only four from
sample 1, and expresses a coherent germline/meiotic program; it is retained as
a contamination control and excluded from hemocyte inference.

## Annotation policy

Cell labels require concordance among replicate-consistent marker genes,
functional annotation and GO enrichment, mapped Drosophila orthologs, the
cluster hierarchy, and established CPB morphology. Cell-cycle and stress
programs are treated as states rather than independent hemocyte classes.
Cluster 5 is reported as `Prohemocytes (provisional)` because morphology,
immunostaining, abundance, and recovery across both samples support this
interpretation, whereas its low RNA content and lack of a specific positive
marker prevent a definitive barcode-level assignment.

## Main outputs

```text
results/combined_annotated.rds
results/QC/
results/resolution_validation/
results/cell_cycle/
results/cluster_validation/
results/replicate_consensus_GO/
results/final_annotation/
results/figures/figure_1.pdf ... figure_5.pdf
results/shinycell2_cpb/                 # optional
```

## License

The analysis code is released under the MIT License. Input data and annotation
files should be distributed under their applicable data-access terms.
