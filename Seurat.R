library(dplyr)
library(readr)
library(tidyr)
library(Seurat)
library(patchwork)
library(stringr) 
library(ggplot2)
library(rtracklayer)
library(scDblFinder)
library(writexl)

s1.data <- Read10X(data.dir = "file1")
s2.data <- Read10X(data.dir = "file2")

s1 <- CreateSeuratObject(counts = s1.data, project="rep1")
s2 <- CreateSeuratObject(counts = s2.data, project="rep2")

gtf_data <- read.table("CPB_annotation.gtf", sep="\t", stringsAsFactors=FALSE, quote="")

sce <- as.SingleCellExperiment(s1)
sce_2 <- as.SingleCellExperiment(s2)

set.seed(1234)

# Run scDblFinder
sce <- scDblFinder(sce) # 762 (6.5%) doublets called 
sce_2 <- scDblFinder(sce_2) # 370 (5.2%) doublets called 

# Add back to Seurat
s1$scDblFinder <- sce$scDblFinder.class
s2$scDblFinder <- sce_2$scDblFinder.class

s1_singlet <- subset(s1, subset = scDblFinder == "singlet")
s2_singlet <- subset(s2, subset = scDblFinder == "singlet")

s1_singlet$replicate <- "rep1"
s2_singlet$replicate <- "rep2"
combined <- merge(s1_singlet, y = s2_singlet, add.cell.ids = c("rep1", "rep2"), project = "merged_reps")

#### removes cells with >= 15% mitochondrial genes ----
mt_genes <- gtf_data[gtf_data$V1=="MZ189364",]
mt_genes <- mt_genes[mt_genes$V3=="CDS",]
vec_mt_genes <- str_extract(mt_genes$V9, 'gene_id "\\w*"')
vec_mt_genes <- gsub('gene_id "',"",vec_mt_genes)
vec_mt_genes <- gsub('"', '',vec_mt_genes)

vec_mt_genes_red <- vec_mt_genes[which(vec_mt_genes %in% Features(combined))]
combined[["percent.mt"]] <- PercentageFeatureSet(combined, features = vec_mt_genes_red, assay = 'RNA')
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
hist(combined$percent.mt, breaks = 100, main = "Mitochondrial % Distribution")
mean(combined$percent.mt) # 9.5%
median(combined$percent.mt) # 4.3%
combined <- subset(combined, subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 15)

combined <- JoinLayers(combined) # layers need to be combined to run GetAssayData
expr_matrix_seurat_filt_combined <- GetAssayData(combined, assay = "RNA", layer = "counts")
genes_use_combined <- rowSums(expr_matrix_seurat_filt_combined > 0) >= 3
rowSums(expr_matrix_seurat_filt_combined) # 14736
sum(genes_use_combined) # 13651
expr_matrix_seurat_filt_combined <- expr_matrix_seurat_filt_combined[genes_use_combined, ]

combined <- CreateSeuratObject(counts = expr_matrix_seurat_filt_combined, meta.data = combined@meta.data)
combined[["RNA"]] <- split(combined[["RNA"]], f = combined$replicate)

#total number of cells after filtration
ncol(combined) # 15,477

set.seed(1234)

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 3000)
plot1 <- VariableFeaturePlot(combined)
top10 <- head(VariableFeatures(combined), 10)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

combined <- ScaleData(combined)
combined <- RunPCA(combined)

#### Integrate the replicates with Harmony ----
combined <- IntegrateLayers(
  object = combined, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.harmony",
  verbose = FALSE
)

#### Cluster cells ----
ElbowPlot(combined, ndims = 50, reduction = "pca")
combined <- FindNeighbors(combined, reduction = "integrated.harmony", dims = 1:10) 
combined <- FindClusters(combined, cluster.name = "harmony_clusters", resolution = 0.5)
table(Idents(combined))

combined <- RunUMAP(combined, reduction = "integrated.harmony", dims = 1:10)
DimPlot(combined, reduction = "umap", label = TRUE)

ncol(combined) # 15,477

combined <- BuildClusterTree(combined)
PlotClusterTree(combined)

combined <- JoinLayers(combined)

saveRDS(combined, "combined.RDS")

expr_matrix <- GetAssayData(combined, slot = "counts") 
genes_expressed <- rownames(expr_matrix)[Matrix::rowSums(expr_matrix) > 0] # 13,651 genes expressed
writeLines(genes_expressed, "expressed_genes.txt") # need the expressed genes as the background for the GO enrichment analysis

#### Determine marker genes ----
cluster_markers05 <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.5) 
cluster_markers05 <- dplyr::filter(cluster_markers05, p_val_adj < 0.05)

# Check function marker genes
cpb_anno <- read.table("LdecV5_functional_annotation.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE, quote="")
colnames(cpb_anno)[1] <- "gene"

anno_cluster_markers <- left_join(cluster_markers05, cpb_anno, by="gene")

write.table(anno_cluster_markers, "cluster_markers_res05_pct010_logfc050.txt", sep = "\t", quote=FALSE, row.names=FALSE)

top10_markers <- anno_cluster_markers %>%
  dplyr::filter(pct.1 >= 0.2) %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>% 
  arrange(cluster, desc(avg_log2FC))

#write_xlsx(top10_markers, "Top10_cluster_markers_pct010_res05.xlsx")
