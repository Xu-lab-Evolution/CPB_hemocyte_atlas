library(clusterProfiler)
library(dplyr)
library(readr)
library(tidyr)
library(reshape2) # to use decast
library(stringr) #for str_detect

cell_clusters <- read.table("cluster_markers_res05_pct010_logfc050.txt", sep="\t", quote="", stringsAsFactors=FALSE, header=TRUE)
colnames(cell_clusters)[7] <- "Geneid"

# increase the log fold change threshold
cell_clusters <- cell_clusters[cell_clusters$avg_log2FC >= 1 & cell_clusters$p_val_adj < 0.05,]

# universe of expressed genes in the single-cell dataset
expressed_genes <- read.table("expressed_genes.txt", sep="\t", stringsAsFactors=FALSE, header=FALSE)

go_anno <- read.table("GO_terms_ancestors_LDECv5.csv", sep="\t", quote="", stringsAsFactors=FALSE, header=TRUE) 
colnames(go_anno)[1] <- "GO_term"

gene2go <- read.table('LdecV5_functional_annotation.txt', sep="\t", quote="", 
                      stringsAsFactors=FALSE, header=TRUE)
gene2go <- gene2go[, c(1,4)]

gogene <- read.table('E:/new_assembly/duo_go_genes_LdecV5.txt', sep="\t", stringsAsFactors=FALSE, header=TRUE)

#remove obsolete terms
go_anno <- go_anno[!grepl("OBSOLETE", go_anno$Definition),]
gogene <- semi_join(gogene, go_anno, by="GO_term") 

go_molfunc <- go_anno[go_anno$Ontology=="molecular_function", ]
colnames(go_molfunc)[1] <- "GO_term"
go_biolproc <- go_anno[go_anno$Ontology=="biological_process", ]
colnames(go_biolproc)[1] <- "GO_term"
go_cellcomp <- go_anno[go_anno$Ontology=="cellular_component", ]

cell_types <- unique(cell_clusters$cluster)

# Associate the modules with the GO term ontologies they have genes in
gene_ontology <- left_join(gogene, go_anno, by="GO_term") # associate gene-go with go ontologies
gene_ontology <- gene_ontology[, -c(4:6)]
gene_ontology$Geneid <- gsub('.t[0-9]', '', gene_ontology$Geneid)

color_ont <- left_join(gene_ontology, cell_clusters, by="Geneid")
color_ontology <- color_ont %>% drop_na(cluster) %>% unique() 

vec_molfunc <- color_ontology$cluster[color_ontology$Ontology=="molecular_function"]
vec_biolproc <- color_ontology$cluster[color_ontology$Ontology=="biological_process"]

go2names <- go_anno[, c(1,3)] #get the name of the go

##### ORA -----
#### For molecular function ----

# extract association GO term-gene for GO terms corresponding to molecular functions
gogene_molfunc <- inner_join(go_molfunc, gogene, by="GO_term")
gogene_molfunc_red <- gogene_molfunc[, -c(2:5)]

go_enrichment <- function(cluster_name){
  clust_vector <- cell_clusters[cell_clusters$cluster==cluster_name, ]
  
  if (cluster_name %in% vec_molfunc){
    result <- enricher(as.character(clust_vector[,7]), # Geneid is in the 7th column
                       universe=expressed_genes$x,
                       TERM2GENE=gogene_molfunc_red,
                       TERM2NAME=go2names,
                       pvalueCutoff = 0.05)
    return(c(result, cluster_name))
  }
}

res_data <- lapply(cell_types, go_enrichment)

table_fin_molfunc <- data.frame()

enrich_tables <- lapply(res_data, function(x) {
  if (is.null(x[[1]]) || nrow(as.data.frame(x[[1]])) == 0) return(NULL)
  df <- as.data.frame(x[[1]])
  df$Module <- as.character(x[[2]])
  df
})

table_fin_molfunc <- bind_rows(enrich_tables)
table_fin_molfunc <- table_fin_molfunc %>% relocate(Module, .before = "ID")
table_fin_molfunc <- table_fin_molfunc %>% mutate(Above3genes = Count >= 3)

#write_delim(table_fin_molfunc, file = "GO_enrich_cluster_res05_pct010_logfc1_molfunct.txt", delim = "\t")

#### For biological process ----

# extract association GO term-gene for GO terms corresponding to biological process
gogene_biolproc <- inner_join(go_biolproc, gogene, by="GO_term")
gogene_biolproc_red <- gogene_biolproc[, -c(2:5)]

go_enrichment_biolproc <- function(cluster_name){
  clust_vector <- cell_clusters[cell_clusters$cluster == cluster_name, ]
  
  if (cluster_name %in% vec_biolproc){
    result <- enricher(as.character(clust_vector[,7]), 
                       universe = expressed_genes$x,
                       TERM2GENE = gogene_biolproc_red,
                       TERM2NAME = go2names,
                       pvalueCutoff = 0.05)
    return(list(result = result, cluster = cluster_name))
  } else {
    return(NULL)
  }
}

res_data_biolproc <- lapply(cell_types, go_enrichment_biolproc)

table_fin_biolproc <- data.frame()

enrich_tables_biolproc <- lapply(res_data_biolproc, function(x) {
  if (is.null(x$result) || nrow(as.data.frame(x$result)) == 0) return(NULL)
  df <- as.data.frame(x$result)
  df$Module <- as.character(x$cluster)
  df
})

table_fin_biolproc <- bind_rows(enrich_tables_biolproc)
table_fin_biolproc <- table_fin_biolproc %>% relocate(Module, .before = "ID")
table_fin_biolproc <- table_fin_biolproc %>% mutate(Above3genes = Count >= 3)

#write_delim(table_fin_biolproc, file = "GO_enrich_clusters_res05_pct010_logfc1_biolproc.txt", delim = "\t")
