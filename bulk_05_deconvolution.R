setwd("")

library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(rstatix)
library(Matrix)
library(sva)
library(Seurat)
library(BisqueRNA)
library(SingleCellExperiment)
library(patchwork)
library(viridis)
library(scDblFinder)
library(GEOquery)


bulk_tpm <- read.csv("")
rownames(bulk_tpm) <- bulk_tpm[, 1]
bulk_tpm <- bulk_tpm[, 2:ncol(bulk_tpm)]
bulk_tpm <- as.data.frame(bulk_tpm)

metadata_bulk <- read.csv("")
metadata_bulk$group <- factor(metadata_bulk$group, levels = c("CTRL", "BD"))
metadata_bulk$batch <- as.character(metadata_bulk$batch)
str(metadata_bulk)

mean_expression <- apply(bulk_tpm, 2, mean, na.rm = TRUE)
barplot(mean_expression,
        xlab = "Sample", ylab = "Mean expression",
        main = "Mean expression across samples",
        names.arg = c(1:dim(bulk_tpm)[[2]]), cex.names = 0.7)

dim1 <- 1
dim2 <- 2

group <- metadata_bulk$group
batch <- metadata_bulk$batch

pca.object <- prcomp(t(bulk_tpm), scale. = TRUE)

pcaData <- as.data.frame(pca.object$x[, c(dim1, dim2)])
pcaData <- cbind(pcaData, group, batch)

percentVar <- round(100*summary(pca.object)$importance[2, c(dim1, dim2)], 0)

p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = batch, shape = group)) + 
  geom_point(size = 4) +
  geom_text_repel(aes(label = rownames(pcaData)), size = 2, max.overlaps = Inf) +
  xlab(paste0("PC", dim1, ": ",percentVar[1],"% variance")) +
  ylab(paste0("PC", dim2, ": ",percentVar[2],"% variance")) +
  ggtitle("PCA before batch correction") +
  coord_fixed(ratio=1.0) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)

print(p)
ggsave()


bulk_tpm <- ComBat(bulk_tpm, batch = metadata_bulk$batch)


pca.object <- prcomp(t(bulk_tpm), scale. = TRUE)

pcaData <- as.data.frame(pca.object$x[, c(dim1, dim2)])
pcaData <- cbind(pcaData, group, batch)

percentVar <- round(100*summary(pca.object)$importance[2, c(dim1, dim2)], 0)

p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = batch, shape = group)) + 
  geom_point(size = 4) +
  geom_text_repel(aes(label = rownames(pcaData)), size = 2, max.overlaps = Inf) +
  xlab(paste0("PC", dim1, ": ",percentVar[1],"% variance")) +
  ylab(paste0("PC", dim2, ": ",percentVar[2],"% variance")) +
  ggtitle("PCA after batch correction") +
  coord_fixed(ratio=1.0) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)

print(p)
ggsave()


mean_expression <- apply(bulk_tpm, 2, mean, na.rm = TRUE)
barplot(mean_expression,
        xlab = "Sample", ylab = "Mean expression",
        main = "Mean expression across samples",
        names.arg = c(1:dim(bulk_tpm)[[2]]), cex.names = 0.7)


gse <- getGEO("GSE213982", GSEMatrix = TRUE)
metadata <- gse$`GSE213982-GPL24676_series_matrix.txt.gz`@phenoData@data
metadata <- cbind(metadata$title, metadata$`group:ch1`, metadata$description) %>% as.data.frame()
colnames(metadata) <- c("donor_id", "group", "batch")

data_sc.data <- ReadMtx(
  mtx = "",
  features = "",
  cells = "",
  cell.column = 1,
  feature.column = 1,
  skip.cell = 1,
  skip.feature = 1
  )

data_sc <- CreateSeuratObject(counts = data_sc.data, project = "data_sc", min.cells = 10, min.features = 400)
head(data_sc@meta.data)
rm(data_sc.data)

cell_names <- Cells(data_sc)
control_cells <- grep("F7\\.|F10\\.|F13\\.|F21\\.|F22\\.|F30\\.|F31\\.|F33\\.|F34\\.|F35\\.|F36\\.|F38\\.", 
                      cell_names, value = TRUE)
substr(control_cells, 1, 3) %>% unique()
data_sc <- subset(data_sc, cells = control_cells)
data_sc
cell_names <- Cells(data_sc)

VlnPlot(data_sc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(data_sc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

genes <- rownames(data_sc)
genes <- genes[!grepl("\\.", genes)]
data_sc <- subset(data_sc, features = genes)

cm <- GetAssayData(data_sc, assay = "RNA", layer = "counts")
gene_counts <- Matrix::rowSums(cm > 0)
genes <- names(gene_counts[gene_counts >= 10])
rm(cm)
 
data_sc <- subset(data_sc, features = genes, 
                  subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & nCount_RNA < 30000)

VlnPlot(data_sc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(data_sc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

data_sc
# An object of class Seurat 
# 20270 features across 15537 samples within 1 assay 
# Active assay: RNA (20270 features, 0 variable features)
# 1 layer present: counts

data_sc <- NormalizeData(data_sc, normalization.method = "LogNormalize", scale.factor = 10000)

data_sc <- FindVariableFeatures(data_sc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(data_sc), 10)

plot1 <- VariableFeaturePlot(data_sc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(data_sc)
data_sc <- ScaleData(data_sc, features = all.genes)

data_sc <- RunPCA(data_sc, features = VariableFeatures(object = data_sc))

ElbowPlot(data_sc, ndims = 50)

data_sc <- FindNeighbors(data_sc, dims = 1:15)
data_sc <- FindClusters(data_sc, resolution = 0.2)

data_sc <- RunUMAP(data_sc,dims = 1:15)

p <- DimPlot(data_sc, reduction = "umap", label = TRUE)
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

data_sc
# An object of class Seurat 
# 20270 features across 15537 samples within 1 assay 
# Active assay: RNA (20270 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap

saveRDS(data_sc, "")


data_sc <- readRDS("")

sce <- as.SingleCellExperiment(data_sc, assay = "RNA")

set.seed(123)
sce <- scDblFinder(sce)

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
# 14315    1222 

p <- DimPlot(data_sc, reduction = "umap", group.by = "scDblFinder.class") +
  ggplot2::ggtitle("Doublets classification on UMAP")
print(p)
ggsave()

data_sc <- subset(data_sc, subset = scDblFinder.class == "singlet")
head(data_sc@meta.data)
length(data_sc@meta.data$scDblFinder.class)

rm(sce)

data_sc <- FindNeighbors(data_sc, dims = 1:15)
data_sc <- FindClusters(data_sc, resolution = 0.2)

data_sc <- RunUMAP(data_sc, dims = 1:15)

p <- DimPlot(data_sc, reduction = "umap", label = TRUE)
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

data_sc
# An object of class Seurat 
# 20270 features across 14315 samples within 1 assay 
# Active assay: RNA (20270 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap


data_sc.markers <- FindAllMarkers(data_sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data_sc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top10 <- data_sc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pdf()
DoHeatmap(subset(data_sc, downsample = 500), features = top10$gene) + NoLegend()
dev.off()


my_color_1 <- c("red3", "dodgerblue3", "mediumpurple1", "gold3", "green4", "lightpink2", "chocolate4", "grey")


cell_types <- c("Oli", "ExN", "Ast", "ExN", "OPC", "InN", "InN", "ExN", "InN", "End", "ExN", "ExN", "InN", "Mic", "End", "ExN")
cell_types <- factor(cell_types, levels = c("ExN", "InN", "Ast", "Mic", "Oli", "OPC", "End"))
data_sc$cell_type <- cell_types[data_sc$seurat_clusters] 
Idents(data_sc) <- data_sc$cell_type
cell_names <- Cells(data_sc)

p <- DimPlot(data_sc, reduction = "umap", label = FALSE, group.by = "cell_type") + 
  ggtitle("Cell type annotation") + scale_color_manual(values = my_color_1)
print(p)
ggsave()


all_markers <- FindAllMarkers(data_sc, group.by = "cell_type", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top10 <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

common_genes <- subset(all_markers, subset = rownames(all_markers) %in% rownames(bulk_tpm))
top40 <- all_markers %>% group_by(cluster) %>% top_n(n = 40, wt = avg_log2FC)
write.csv(top40, "")

pdf()
DoHeatmap(subset(data_sc, downsample = 500), features = top10$gene, 
          group.by = "cell_type", group.colors = my_color_1) + NoLegend()
dev.off()

markers <- top40$gene %>% as.character()


pie_plot <- function(data, title) {
  df <- as.data.frame(table(data_sc@meta.data$cell_type))
  colnames(df) <- c("Type", "Count")
  df <- df %>%
    mutate(Percent = 100 * Count / sum(Count),
           Label = sprintf("%.1f%%", Percent))
  
  ggplot(df, aes(x = "", y = Percent, fill = Type)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y") +
    geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 3) +
    ggtitle(title) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
}

p <- pie_plot(data_sc, "Cell Type Category")
print(p)


sc_counts <- GetAssayData(data_sc, layer = "counts")

sc_eset <- ExpressionSet(assayData = as.matrix(sc_counts))
pData(sc_eset)$cellType <- data_sc$cell_type
pData(sc_eset)$SubjectName <- cell_names

bulk_eset <- ExpressionSet(assayData = as.matrix(bulk_tpm))

res <- ReferenceBasedDecomposition(bulk.eset = bulk_eset, sc.eset = sc_eset, markers = markers, use.overlap = FALSE)

prop <- res$bulk.props
write.csv(prop, "")


prop <- prop[, metadata_bulk$SRR_id]
all(metadata_bulk$SRR_id == colnames(prop))
sample_to_group <- tibble(sample = colnames(prop), group = metadata_bulk$group)
prop <- as.data.frame(prop) %>% rownames_to_column("cell_type")

data_for_boxplot <- prop %>% 
  pivot_longer(cols = 2:ncol(prop), values_to = "proportion", names_to = "sample") %>% 
  left_join(sample_to_group, by = "sample")

data_for_boxplot$cell_type <- factor(data_for_boxplot$cell_type, 
                                     levels = c("ExN", "InN", "Ast", "Mic", "Oli", "OPC", "End"))

stat_test <- data_for_boxplot %>%
  group_by(cell_type) %>%
  t_test(proportion ~ group)

y_pos <- data_for_boxplot %>%
  group_by(cell_type) %>%
  summarise(y.position = max(proportion) + 0.05)

stat_test <- stat_test %>%
  left_join(y_pos, by = "cell_type")

p <- ggplot(data_for_boxplot) + 
  facet_wrap(~cell_type, scales = "free_y", ncol = 4) +
  geom_boxplot(aes(x = group, y = proportion, fill = group), alpha = 0.7) + 
  geom_jitter(aes(x = group, y = proportion, fill = group), alpha = 0.5, width = 0.2, size = 1) + 
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  stat_pvalue_manual(stat_test, label = "p", tip.length = 0, size = 3)
print(p)
ggsave()
