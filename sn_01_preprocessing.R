setwd("")

library(Matrix)
library(dplyr)
library(Seurat)
library(SummarizedExperiment)
library(tidyverse)
library(ggplot2)
library(viridis)
library(scDblFinder)
library(harmony)
library(GEOquery)


my_color_1 <- c("red3", "dodgerblue3", "mediumpurple1", "gold3", "green4", "lightpink2", "chocolate4")


samples <- basename(list.files("", full.names = TRUE))

objs <- list()
for (sample in samples) {
  data <- Read10X_h5(paste0("", sample))
  sample_name <- substr(sample, 1, 10)
  colnames(data) <- paste(sample_name, colnames(data), sep = "-")
  data <- CreateSeuratObject(counts = data, project = sample_name, min.cells = 10, min.features = 300)
  
  data$sample <- sample
  objs[[sample]] <- data
}

data_sc <- Reduce(function(x, y) merge(x, y), objs)
data_sc <- JoinLayers(data_sc)

rm(data)
rm(objs)

data_sc$orig.ident <- substr(data_sc$orig.ident, 1, 10)
data_sc$sample <- substr(data_sc$sample, 1, 10)

gse <- getGEO("GSE168408", GSEMatrix = TRUE)
metadata_1 <- gse$`GSE168408-GPL21697_series_matrix.txt.gz`@phenoData@data
metadata_2 <- gse$`GSE168408-GPL24676_series_matrix.txt.gz`@phenoData@data
metadata_1 <- cbind(metadata_1$title, metadata_1$geo_accession) %>% as.data.frame()
metadata_2 <- cbind(metadata_2$title, metadata_2$geo_accession) %>% as.data.frame()
colnames(metadata_1) <- c("stage", "sample_id")
colnames(metadata_2) <- c("stage", "sample_id")
metadata <- rbind(metadata_1, metadata_2)
rm(gse)

samples <- substr(samples, 1, 10)
metadata <- subset(metadata, sample_id %in% samples)
metadata$stage <- str_sub(metadata$stage, 1, -10)
metadata$stage <- factor(metadata$stage, levels = c("ga22",  "ga24",  "ga34",  "1m",    "2m",    "3m",    
                                                    "4m",    "6m",    "10m",   "1y_a",  "1y_b",  "1y_c",  
                                                    "1y_d",  "2y_a",  "4y",    "6y",    "8y",    "10y_a", 
                                                    "10y_b", "12y",   "14y",   "16y",   "17y",   "20y_a", 
                                                    "20y_b", "25y"))
metadata$chemistry <- "v3"
metadata$chemistry[which(metadata$sample_id == "GSM5138516" | metadata$sample_id == "GSM5138535")] <- "v2"

write.csv(metadata, "", row.names = FALSE)

metadata <- metadata[match(data_sc$orig.ident, metadata$sample_id), c("stage", "chemistry")] %>% as.data.frame()
colnames(metadata) <- c("stage", "chemistry")
rownames(metadata) <- colnames(data_sc)

data_sc <- AddMetaData(data_sc, metadata = metadata)
head(data_sc@meta.data[, c("orig.ident", "stage", "chemistry")])

length(data_sc@meta.data$nFeature_RNA)

genes <- rownames(data_sc)
genes <- genes[!grepl("\\.", genes)]
data_sc <- subset(data_sc, features = genes)

gene_counts <- rowSums(data_sc[["RNA"]]$counts > 0)
genes <- names(gene_counts[gene_counts >= 300])

data_sc[["percent.mt"]] <- PercentageFeatureSet(data_sc, pattern = "^MT-")
VlnPlot(data_sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

data_sc <- subset(data_sc, features = genes, 
                  subset = nFeature_RNA > 800 & nFeature_RNA < 7000 & percent.mt < 5)

length(data_sc@meta.data$nFeature_RNA)


sce <- as.SingleCellExperiment(data_sc, assay = "RNA")

set.seed(123)
sce <- scDblFinder(sce, samples = sce$sample, clusters = TRUE, dbr = 0.2, verbose = FALSE)

res <- data.frame(
  cell = colnames(sce),
  scDblFinder.class  = sce$scDblFinder.class,
  scDblFinder.score  = sce$scDblFinder.score,
  scDblFinder.weighted = sce$scDblFinder.weighted,
  stringsAsFactors = FALSE
)

res <- res[match(colnames(data_sc), res$cell), ]
stopifnot(all(res$cell == colnames(data_sc)))
data_sc$scDblFinder.class <- res$scDblFinder.class
data_sc$scDblFinder.score <- res$scDblFinder.score
data_sc$scDblFinder.weighted<- res$scDblFinder.weighted

table(data_sc$scDblFinder.class)
# singlet doublet 
# 136678   15388 
table(sce$orig.ident, sce$scDblFinder.class)
# singlet doublet
# GSM5138509    8755    2211
# GSM5138511    9307     976
# GSM5138513    6575    1116
# GSM5138514    9625    1002
# GSM5138516    6242     651
# GSM5138517    6985     742
# GSM5138519    3470     342
# GSM5138520    6241     594
# GSM5138522    4887     509
# GSM5138524    4790     412
# GSM5138525    5033     520
# GSM5138527    1754     138
# GSM5138531    4141     402
# GSM5138533    7607     744
# GSM5138535    3288     327
# GSM5138537     789      80
# GSM5138538    3295     343
# GSM5138540    5857     560
# GSM5138541    2119     173
# GSM5138543    5457     538
# GSM5138545    5617     544
# GSM5138546    5229     591
# GSM5138547    6879     666
# GSM5138549    5490     511
# GSM6499429    3225     301
# GSM6499430    4021     395

data_sc <- subset(data_sc, subset = scDblFinder.class == "singlet")
length(data_sc@meta.data$scDblFinder.class)

rm(sce)
rm(res)

saveRDS(data_sc, file = "")


setwd("")

library(Matrix)
library(dplyr)
library(Seurat)
library(SummarizedExperiment)
library(tidyverse)
library(ggplot2)
library(viridis)
library(harmony)


data_sc <- readRDS("")

meta <- data_sc@meta.data$orig.ident %>% unique() %>% as.data.frame()
colnames(meta) <- "orig.ident"
meta$region <- c("dlPFC", "dlPFC", "dlPFC", "BA8",   "dlPFC", "dlPFC", "dlPFC", "dlPFC", 
                 "dlPFC", "BA8",   "dlPFC", "BA8",   "dlPFC", "dlPFC", "dlPFC", "dlPFC", 
                 "dlPFC", "dlPFC", "dlPFC", "BA10",  "BA8",   "BA8",   "BA8",   "BA10", 
                 "dlPFC", "dlPFC")
meta

meta <- meta[match(data_sc$orig.ident, meta$orig.ident), c("region")] %>% as.data.frame()
colnames(meta) <- c("region")
rownames(meta) <- colnames(data_sc)

data_sc <- AddMetaData(data_sc, metadata = meta)
head(data_sc@meta.data[, c("orig.ident", "region")])

data_sc <- SCTransform(data_sc, vars.to.regress = c("chemistry", "percent.mt"))

data_sc <- RunPCA(data_sc, assay = "SCT", features = VariableFeatures(object = data_sc))

p <- ElbowPlot(data_sc, ndims = 50)
print(p)
ggsave()

data_sc <- FindNeighbors(data_sc, dims = 1:8)
data_sc <- FindClusters(data_sc, resolution = 0.5)

data_sc <- RunUMAP(data_sc, dims = 1:8)

p <- DimPlot(data_sc, reduction = "umap", label = TRUE)
print(p)
ggsave()

p <- DimPlot(data_sc, group.by = "stage") + scale_color_viridis_d(option = "inferno")
print(p)
ggsave()

data_sc <- RunHarmony(data_sc, assay = "SCT", group.by.vars = c("chemistry", "region")) %>% 
  RunUMAP(reduction = "harmony", dims = 1:8) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:8) %>% 
  FindClusters(resolution = 0.5)

p <- DimPlot(data_sc, reduction = "umap", label = TRUE)
print(p)
ggsave()

p <- DimPlot(data_sc, group.by = "stage") + scale_color_viridis_d(option = "inferno")
print(p)
ggsave()

plot1 <- VlnPlot(data_sc, features = "nCount_RNA", pt.size = 0, group.by = "seurat_clusters") +
  NoLegend()
plot2 <- VlnPlot(data_sc, features = "nFeature_RNA", pt.size = 0, group.by = "seurat_clusters") +
  NoLegend()
p <- plot1 + plot2
print(p)
ggsave()

p <- FeaturePlot(data_sc, features = c("nFeature_RNA")) + scale_color_gradientn(colours = viridis::turbo(100))
print(p)
ggsave()

p <- FeaturePlot(data_sc, features = c("nCount_RNA")) + scale_color_gradientn(colours = viridis::turbo(100))
print(p)
ggsave()


saveRDS(data_sc, file = "")


setwd("")

library(Matrix)
library(dplyr)
library(Seurat)
library(SummarizedExperiment)
library(tidyverse)
library(ggplot2)
library(viridis)
library(harmony)


data_sc <- readRDS("")

data_sc.markers <- FindAllMarkers(data_sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data_sc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top20 <- data_sc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top10 <- data_sc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top3 <- data_sc.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
write.csv(top20, "")
write.csv(top10, "")
write.csv(top3, "")

pdf()
DoHeatmap(subset(data_sc, downsample = 500), features = top10$gene) + NoLegend()
dev.off()

genes <- top3$gene %>% as.character() %>% unique()
p <- DotPlot(data_sc, features = genes) + RotatedAxis() + scale_color_gradientn(colours = c("lightgray", "orange", "red", "darkred", "black"))
print(p)
ggsave()

markers_for_plot <- c("DCX",     "MEIS2",   "SOX6",    "PAX6",    "DLX2",    "RBFOX3",  "SLC17A7", "CUX2",    
                      "RORB",    "TLE4",    "GAD1",    "ADARB2",  "LAMP5",   "LHX6",    "RELN",    "VIP",     
                      "SST",     "PVALB",   "AQP4",    "DOCK8",   "MOBP",    "PDGFRA",  "CLDN5",   "LAMA2")

output_dir <- ""

for (marker in markers_for_plot) {
  p <- FeaturePlot(data_sc, features = marker) + scale_color_gradientn(colours = c("lightgray", "orange", "red", "darkred", "black"))
  ggsave()
}


cell_type_1 <- c("ExN", "ExN", "OPC", "Oli", "Ast", 
                 "Oli", "Ast", "InN", "ExN", "InN", 
                 "InN", "Mic", "ExN", "Ast", "ExN", 
                 "ExN", "InN", "End", "ExN", "Mic")
data_sc$cell_type_1 <- cell_type_1[data_sc$seurat_clusters]
data_sc$cell_type_1 <- factor(data_sc$cell_type_1, 
                              levels = c("ExN", "InN", "Ast", "Mic", "Oli", "OPC", "End"))
Idents(data_sc) <- data_sc$cell_type_1

p <- DimPlot(data_sc, reduction = "umap", label = FALSE, group.by = "cell_type_1") + 
  ggtitle("Cell type annotation on UMAP") + scale_color_manual(values = my_color_1)
print(p)
ggsave()

p <- DimPlot(data_sc, reduction = "umap", label = FALSE, group.by = "cell_type_1") + 
  NoLegend() + NoAxes() + 
  ggtitle("Cell type annotation") + scale_color_manual(values = my_color_1)
print(p)
ggsave()


data_sc.markers <- FindAllMarkers(data_sc, group.by = "cell_type_1", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data_sc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top20 <- data_sc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top10 <- data_sc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top3 <- data_sc.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
write.csv(top20, "")
write.csv(top10, "")
write.csv(top3, "")

pdf()
DoHeatmap(subset(data_sc, downsample = 500), features = top10$gene, 
          group.by = "cell_type_1", group.colors = my_color_1) + NoLegend()
dev.off()

genes <- top3$gene %>% as.character() %>% unique()
p <- DotPlot(data_sc, features = genes, group.by = "cell_type_1") + RotatedAxis() + scale_color_gradientn(colours = c("lightgray", "orange", "red", "darkred", "black"))
print(p)
ggsave()


cluster_sample_df <- as.data.frame(table(data_sc@meta.data$cell_type_1, data_sc@meta.data$stage))
colnames(cluster_sample_df) <- c("CellType", "Stage", "CellCount")

sample_totals <- aggregate(CellCount ~ Stage, data = cluster_sample_df, sum)

cluster_sample_df <- merge(cluster_sample_df, sample_totals, by = "Stage", suffixes = c("", "_Total"))
cluster_sample_df$Proportion <- cluster_sample_df$CellCount / cluster_sample_df$CellCount_Total

p <- ggplot(cluster_sample_df, aes(x = Stage, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) + 
  labs(title = "Cell proportion", fill = "Cell type") +
  scale_fill_manual(values = my_color_1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(p)
ggsave()


data_sc
# An object of class Seurat 
# 32858 features across 136678 samples within 2 assays 
# Active assay: SCT (16429 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, umap, harmony


saveRDS(data_sc, file = "")


p <- DimPlot(data_sc, cells.highlight = WhichCells(data_sc, expression = NPAS4 > 0), order = TRUE) + 
  NoLegend() + ggtitle("NPAS4 expression > 0") + theme(plot.title = element_text(hjust = 0.5))
print(p)
ggsave()

p <- DimPlot(data_sc, cells.highlight = WhichCells(data_sc, expression = NR4A1 > 0), order = TRUE) + 
  NoLegend() + ggtitle("NR4A1 expression > 0") + theme(plot.title = element_text(hjust = 0.5))
print(p)
ggsave()


data_sc_e <- subset(data_sc, cell_type_1 == "ExN")
dir.create("")

data_sc_e <- RunPCA(data_sc_e, assay = "SCT")

p <- ElbowPlot(data_sc_e, ndims = 50)
print(p)
ggsave()

data_sc_e <- FindNeighbors(data_sc_e, dims = 1:6)
data_sc_e <- FindClusters(data_sc_e, resolution = 0.5)

data_sc_e <- RunUMAP(data_sc_e, dims = 1:6)

data_sc_e <- RunHarmony(data_sc_e, assay = "SCT", group.by.vars = c("chemistry", "region")) %>% 
  RunUMAP(reduction = "harmony", dims = 1:6) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:6) %>% 
  FindClusters(resolution = 0.5)

p <- DimPlot(data_sc_e, reduction = "umap", label = TRUE)
print(p)
ggsave()

plot1 <- VlnPlot(data_sc_e, features = "nCount_RNA", pt.size = 0, group.by = "seurat_clusters") +
  NoLegend()
plot2 <- VlnPlot(data_sc_e, features = "nFeature_RNA", pt.size = 0, group.by = "seurat_clusters") +
  NoLegend()
p <- plot1 + plot2
print(p)
ggsave()

data_sc_e <- subset(data_sc_e, seurat_clusters != "12")

data_sc_e <- FindNeighbors(data_sc_e, dims = 1:4)
data_sc_e <- FindClusters(data_sc_e, resolution = 0.5)

data_sc_e <- RunUMAP(data_sc_e, dims = 1:4)

data_sc_e <- RunHarmony(data_sc_e, assay = "SCT", group.by.vars = c("chemistry", "region")) %>% 
  RunUMAP(reduction = "harmony", dims = 1:4) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:4) %>% 
  FindClusters(resolution = 0.5)

p <- DimPlot(data_sc_e, reduction = "umap", label = TRUE)
print(p)
ggsave()
