# Assign established insect hemocyte terminology to validated clusters.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(openxlsx)
})

source(file.path("R", "config.R"))

CONSENSUS_OUT <- file.path(OUT, "replicate_consensus_GO")
ANNOTATION_OUT <- file.path(OUT, "final_annotation")
dir.create(ANNOTATION_OUT, recursive = TRUE, showWarnings = FALSE)

combined <- JoinLayers(readRDS(file.path(
  OUT, "combined_consensus.rds"
)))
DefaultAssay(combined) <- "RNA"

consensus <- read.delim(
  file.path(CONSENSUS_OUT, "cluster_annotations.tsv"),
  check.names = FALSE
)

# Classical names require agreement among CPB morphology, function/GO and
# orthologous markers. State qualifiers prevent activation or cycling from
# being mistaken for an additional classical cell type.
cluster_annotation <- tibble::tribble(
  ~fine_cluster, ~hemocyte_class, ~hemocyte_state, ~cluster_annotation,
  ~short_label, ~annotation_confidence, ~functional_GO_basis,
  ~ortholog_marker_basis, ~morphology_basis, ~limitation,
  "1", "Plasmatocyte", "Proliferating",
  "Proliferating plasmatocytes", "Proliferating plasmatocytes", "moderate",
  "Mitotic cell cycle, spindle, chromosome organization, DNA replication and cell division; this is the most actively dividing population.",
  "Hemolectin/hemocytin-like is detected in 82% and Nimrod/Eater-like in 53% of cells, with broad integrin expression; these support a plasmatocyte relationship beneath the dominant cycling program.",
  "Plasmatocytes are the predominant morphology in fourth-instar CPB hemolymph, but morphology was not linked to individual barcodes.",
  "Consensus markers are almost entirely cell-cycle genes, and the broad parent cluster also contains granulocytes; the lineage assignment is therefore less secure than the proliferative-state assignment.",
  "2", "Granulocyte", "Immune-active",
  "Immune-active granulocytes", "Granulocytes", "moderate-to-high",
  "Innate-immune and biotic-stimulus responses, integrin signaling, actin regulation and adhesion, with dual oxidase and PPO-activating-factor expression.",
  "Rac2-like, integrin beta-PS, integrin alpha-PS3, paxillin and vinculin agree with adhesive/phagocytic granulocyte criteria in other insects.",
  "Granulocytes are established in fourth-instar CPB and are classically granular, adhesive and phagocytic.",
  "A PPO paralog is strongly enriched, but PPO expression is not lineage-exclusive across insects; absence of strong Notch/Pebbled enrichment and the dominant adhesion/immune program argue against assigning this cluster to the oenocytoid branch.",
  "3", "Oenocytoid", "Stress-responsive",
  "Stress-responsive oenocytoids", "Stress-responsive oenocytoids", "high",
  "Heat response and transition-metal transport accompany a strong melanization-associated PPO program.",
  "The principal PPO paralogs, Notch and Pebbled-like expression agree with the conserved oenocytoid/crystal-cell differentiation program.",
  "Oenocytoids are an established, usually round and weakly adhesive CPB hemocyte class associated with phenoloxidase biology.",
  "Stress and metal-handling genes define the state within the class; they should not be interpreted as a separate lineage.",
  "4", "Plasmatocyte", "Adhesive/clotting",
  "Adhesive/clotting plasmatocytes", "Adhesive/clotting plasmatocytes", "high",
  "Cell adhesion, locomotion/motility regulation, receptor signaling, extracellular-matrix organization, calcium binding and clotting.",
  "Hemocytin/hemolectin-like, Nimrod/Eater-like, hemocyte transglutaminases, papilin, laminins and thrombospondin provide the strongest plasmatocyte ortholog support in the dataset.",
  "CPB plasmatocytes are the predominant established morphology and are adhesive/spreading cells involved in cellular encapsulation.",
  "Some adhesion and clotting functions also occur in granulocytes, but the combined Hml/Nim/ECM/transglutaminase evidence is strongest for plasmatocytes.",
  "5", "Prohemocyte", "RNA-low/provisional",
  "Prohemocytes (provisional)", "Prohemocytes (provisional)", "moderate from morphology; low from RNA",
  "No GO analysis is valid because no positive marker passes the replicate-consensus thresholds; RNA and feature counts are markedly lower than in mature clusters.",
  "Transferred Ance, Domeless, DE-cadherin and SPARC candidates are not specific, so no molecular marker establishes the identity.",
  "Independent CPB morphology and immunostaining demonstrate circulating prohemocytes; the cluster frequency (6.1%) is close to the published fourth-instar morphology estimate (6.5%).",
  "Low RNA also occurs in damaged or incomplete droplets. The label applies to a population enriched for prohemocytes and is not a definitive barcode-level assignment.",
  "6", "Oenocytoid", "Differentiating",
  "Differentiating oenocytoids", "Differentiating oenocytoids", "moderate-to-high",
  "Tissue-development regulation, cAMP metabolism and lipid-response terms define a differentiation-associated state within the PPO-rich branch.",
  "Very high Notch, Pebbled-like and the principal PPO paralogs, together with FoxO and Tob2, support an oenocytoid differentiation state.",
  "The class is consistent with established CPB oenocytoids; no matched morphology is available for this state.",
  "The developmental direction is inferred from conserved marker logic and graph proximity, not demonstrated by lineage tracing or time-course data.",
  "7", "Oenocytoid", "Activated antimicrobial/stress",
  "Activated oenocytoids", "Activated oenocytoids", "moderate",
  "Defense response to bacteria, general defense and stress response, with cathepsin, mannose-receptor-like, ficolin and redox genes.",
  "High principal PPO paralogs and a PPO2 ortholog place the rare state in the oenocytoid branch.",
  "Compatible with an activated oenocytoid state rather than a separate CPB morphology.",
  "Only 112 cells are present; activation during collection cannot be excluded.",
  "8", "Non-hemocyte", "Germline/meiotic contamination",
  "Germline/meiotic contaminant", "Germline contaminant", "high",
  "Coherent meiotic and germline program; excluded from hemocyte GO interpretation.",
  "No hemocyte ortholog support.",
  "Not compatible with any established CPB hemocyte morphology.",
  "Only four cells occur in replicate 1 and the branch fails independent recovery."
) |>
  left_join(
    consensus |>
      transmute(
        fine_cluster = as.character(fine_cluster), rep1, rep2, total_cells,
        n_robust_markers, replicate_effect_rho, independently_recovered
      ),
    by = "fine_cluster"
  )

class_definitions <- tibble::tribble(
  ~hemocyte_class, ~standard_morphology_function, ~molecular_criteria_used,
  ~CPB_interpretation,
  "Prohemocyte",
  "Small round cell, high nucleus-to-cytoplasm ratio, scant weakly granular cytoplasm; progenitor competence is possible but proliferation is not obligatory.",
  "Low RNA and depletion of mature effector programs are compatible but not sufficient; species-specific molecular markers require validation.",
  "Cluster 5 is assigned provisionally because morphology, immunostaining, frequency and independent recovery agree, despite absent positive RNA markers.",
  "Plasmatocyte",
  "Adhesive or spreading phagocytic/encapsulating cell, often spindle-shaped; contributes to extracellular matrix, wound repair and clotting.",
  "Hemolectin/hemocytin, Nimrod/Eater-like receptors, transglutaminase, ECM genes and adhesion/motility programs.",
  "Cluster 4 is the reference adhesive/clotting plasmatocyte; cluster 1 is assigned as a proliferating plasmatocyte state with lower lineage confidence.",
  "Granulocyte",
  "Granule-rich adhesive and phagocytic immune effector involved in recognition, spreading, nodulation and encapsulation.",
  "Rac-family activity, integrins, paxillin/vinculin, lectin/recognition genes, lysosomal or redox effectors and innate-immune GO terms.",
  "Cluster 2 meets the functional and ortholog criteria. Its PPO paralog does not override the stronger granulocyte evidence.",
  "Oenocytoid",
  "Usually large and round with limited adhesion; principal source or regulator of prophenoloxidase-dependent melanization in many insects.",
  "High principal PPO paralogs together with Notch/Pebbled/Lozenge-axis evidence and metal/copper handling; class can contain differentiation and activation states.",
  "Clusters 3, 6 and 7 are stress-responsive, differentiating and activated oenocytoid states, respectively.",
  "Spherulocyte",
  "Cell containing characteristic cytoplasmic spherules, often linked to cuticle components and wound repair.",
  "No conserved, validated cross-order marker set is available; morphology is essential.",
  "No replicate-supported CPB cluster can be assigned. The published 0.25% frequency is near the practical detection limit and does not justify relabeling cluster 8."
)

selected_markers <- tibble::tribble(
  ~gene, ~display_name, ~criterion,
  "LDECv5g04564", "Anillin", "Proliferation",
  "LDECv5g10569", "PCNA", "Proliferation",
  "LDECv5g09615", "Topoisomerase II", "Proliferation",
  "LDECv5g13695", "Rac2-like", "Granulocyte",
  "LDECv5g13501", "Integrin beta-PS", "Granulocyte",
  "LDECv5g09037", "Integrin alpha-PS3", "Granulocyte",
  "LDECv5g05220", "Dual oxidase", "Granulocyte",
  "LDECv5g13594", "PPO paralog C", "Granulocyte/Oenocytoid",
  "LDECv5g15036", "PPO paralog A", "Oenocytoid",
  "LDECv5g15038", "PPO paralog B", "Oenocytoid",
  "LDECv5g00509", "Notch", "Oenocytoid differentiation",
  "LDECv5g08619", "Pebbled-like", "Oenocytoid differentiation",
  "LDECv5g03007", "Hemocytin/Hml-like", "Plasmatocyte",
  "LDECv5g08078", "Nimrod/Eater-like", "Plasmatocyte",
  "LDECv5g15611", "Hemocyte transglutaminase", "Plasmatocyte",
  "LDECv5g10644", "Papilin", "Plasmatocyte",
  "LDECv5g06984", "SPARC-like", "Prohemocyte candidate (not specific)",
  "LDECv5g00806", "Ance-like", "Prohemocyte candidate (not specific)"
) |>
  filter(gene %in% rownames(combined))

rna <- GetAssayData(combined, assay = "RNA", layer = "data")
clusters <- as.character(combined$res_0.2)
marker_expression <- bind_rows(lapply(seq_len(nrow(selected_markers)), function(i) {
  bind_rows(lapply(sort(unique(clusters)), function(cluster) {
    values <- rna[selected_markers$gene[i], clusters == cluster]
    tibble(
      fine_cluster = cluster,
      gene = selected_markers$gene[i],
      display_name = selected_markers$display_name[i],
      criterion = selected_markers$criterion[i],
      average_log_normalized_expression = mean(values),
      percent_detected = 100 * mean(values > 0)
    )
  }))
}))

label_map <- setNames(cluster_annotation$cluster_annotation,
                      cluster_annotation$fine_cluster)
short_map <- setNames(cluster_annotation$short_label,
                      cluster_annotation$fine_cluster)
class_map <- setNames(cluster_annotation$hemocyte_class,
                      cluster_annotation$fine_cluster)
state_map <- setNames(cluster_annotation$hemocyte_state,
                      cluster_annotation$fine_cluster)
confidence_map <- setNames(cluster_annotation$annotation_confidence,
                           cluster_annotation$fine_cluster)

combined$cluster_annotation <- unname(label_map[as.character(combined$res_0.2)])
combined$short_label <- unname(short_map[as.character(combined$res_0.2)])
combined$hemocyte_class <- unname(class_map[as.character(combined$res_0.2)])
combined$hemocyte_state <- unname(state_map[as.character(combined$res_0.2)])
combined$annotation_confidence <- unname(
  confidence_map[as.character(combined$res_0.2)]
)

write.table(
  cluster_annotation,
  file.path(ANNOTATION_OUT, "cluster_annotations.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  class_definitions,
  file.path(ANNOTATION_OUT, "standard_hemocyte_classification_criteria.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  marker_expression,
  file.path(ANNOTATION_OUT, "selected_marker_expression_by_cluster.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  selected_markers,
  file.path(ANNOTATION_OUT, "selected_classification_markers.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

openxlsx::write.xlsx(
  list(
    cluster_annotations = cluster_annotation,
    class_definitions = class_definitions,
    selected_markers = selected_markers,
    marker_expression = marker_expression
  ),
  file.path(ANNOTATION_OUT, "standard_hemocyte_annotation.xlsx"),
  overwrite = TRUE
)

saveRDS(
  combined,
  file.path(OUT, "combined_annotated.rds")
)
capture.output(sessionInfo(), file = file.path(ANNOTATION_OUT, "sessionInfo.txt"))
message("Standard hemocyte annotations written to: ", ANNOTATION_OUT)
