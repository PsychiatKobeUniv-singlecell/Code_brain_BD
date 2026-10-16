setwd("")

library(readr)
library(dplyr)
library(edgeR)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(GEOquery)


data <- read.csv("")

rownames(data) <- substr(data[, 1], 1, 15)
data <- data[, 3:ncol(data)]
data <- as.data.frame(data)

meta_id_match <- data.frame()

sraruntable <- read.csv("")
sraruntable <- as.data.frame(bind_cols(sraruntable$Run, sraruntable$GEO_Accession..exp.))
colnames(sraruntable) <- c("SRR_id", "GSE_id")
sraruntable <- sraruntable %>% dplyr::filter(SRR_id %in% colnames(data))
meta_id_match <- rbind(meta_id_match, sraruntable)

sraruntable <- read.csv("")
sraruntable <- as.data.frame(bind_cols(sraruntable$Run, sraruntable$GEO_Accession..exp.))
colnames(sraruntable) <- c("SRR_id", "GSE_id")
sraruntable <- sraruntable %>% dplyr::filter(SRR_id %in% colnames(data))
meta_id_match <- rbind(meta_id_match, sraruntable)

sraruntable <- read.csv("")
sraruntable <- as.data.frame(bind_cols(sraruntable$Run, sraruntable$GEO_Accession..exp.))
colnames(sraruntable) <- c("SRR_id", "GSE_id")
sraruntable <- sraruntable %>% dplyr::filter(SRR_id %in% colnames(data))
meta_id_match <- rbind(meta_id_match, sraruntable)

sraruntable <- read.csv("")
sraruntable <- as.data.frame(bind_cols(sraruntable$Run, sraruntable$GEO_Accession..exp.))
colnames(sraruntable) <- c("SRR_id", "GSE_id")
sraruntable <- sraruntable %>% dplyr::filter(SRR_id %in% colnames(data))
meta_id_match <- rbind(meta_id_match, sraruntable)

metadata_all <- data.frame()

gse <- getGEO("GSE53239", GSEMatrix = TRUE)
sample_info <- pData(phenoData(gse[[2]]))
head(sample_info)

metadata <- as.data.frame(bind_cols(sample_info$geo_accession, sample_info$`disease state:ch1`, 
                                    sample_info$`tissue:ch1`))
colnames(metadata) <- c("GSE_id", "group", "region")
metadata <- metadata %>% dplyr::filter(GSE_id %in% meta_id_match$GSE_id)
metadata$batch <- "1"
head(metadata)
metadata_all <- rbind(metadata_all, metadata)

gse <- getGEO("GSE78936", GSEMatrix = TRUE)
sample_info <- pData(phenoData(gse[[1]]))
head(sample_info)

metadata <- as.data.frame(bind_cols(sample_info$geo_accession, sample_info$`disease stage:ch1`, 
                                    sample_info$`brain region:ch1`))
colnames(metadata) <- c("GSE_id", "group", "region")
metadata <- metadata %>% dplyr::filter(GSE_id %in% meta_id_match$GSE_id)
metadata$batch <- "2"
head(metadata)
metadata_all <- rbind(metadata_all, metadata)

gse <- getGEO("GSE80655", GSEMatrix = TRUE)
sample_info <- pData(phenoData(gse[[1]]))
head(sample_info)

metadata <- as.data.frame(bind_cols(sample_info$geo_accession, sample_info$`clinical diagnosis:ch1`, 
                                    sample_info$`brain region:ch1`))
colnames(metadata) <- c("GSE_id", "group", "region")
metadata <- metadata %>% dplyr::filter(GSE_id %in% meta_id_match$GSE_id)
metadata$batch <- "3"
head(metadata)
metadata_all <- rbind(metadata_all, metadata)

gse <- getGEO("GSE112523", GSEMatrix = TRUE)
sample_info <- pData(phenoData(gse[[1]]))
head(sample_info)

metadata <- as.data.frame(bind_cols(sample_info$geo_accession, sample_info$`dist.dx:ch1`, 
                                    sample_info$source_name_ch1))
colnames(metadata) <- c("GSE_id", "group", "region")
metadata <- metadata %>% dplyr::filter(GSE_id %in% meta_id_match$GSE_id)
metadata$batch <- "4"
head(metadata)
metadata_all <- rbind(metadata_all, metadata)

metadata_all$group <- str_replace_all(metadata_all$group, c("Control" = "CTRL", "Bipolar Disorder" = "BD"))
metadata_all$group <- str_replace_all(metadata_all$group, c("Bipolar" = "BD"))
metadata_all$group <- factor(metadata_all$group, levels = c("CTRL", "BD"))
metadata_all$region <- "PFC"
metadata_all$batch <- as.character(metadata_all$batch)

metadata_all <- merge(meta_id_match, metadata_all, "GSE_id")
str(metadata_all)

write.csv(metadata_all, "", row.names = FALSE)

data <- data[, metadata_all$SRR_id]
all(metadata_all$SRR_id == colnames(data))


library(readr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(magrittr)
library(WGCNA)
library(ggplot2)
library(reshape2)
library(limma)
library(edgeR)
library(sva)


options(stringsAsFactors = FALSE)
allowWGCNAThreads()
set.seed(123)

dge <- DGEList(counts = data)

keep <- filterByExpr(dge)
table(keep)
dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge, method = "TMM") 
dge$samples
data <- cpm(dge$counts, log = TRUE, prior.count = 1)

mean_expression <- apply(data, 2, mean, na.rm = TRUE)
barplot(mean_expression,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples",
        names.arg = c(1:dim(data)[[2]]), cex.names = 0.7)

dim1 <- 1
dim2 <- 2

group <- metadata_all$group
batch <- metadata_all$batch

pca.object <- prcomp(t(data), scale. = TRUE)

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


data <- ComBat(data, batch = metadata_all$batch)


pca.object <- prcomp(t(data), scale. = TRUE)

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

mean_expression <- apply(data, 2, mean, na.rm = TRUE)
barplot(mean_expression,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples",
        names.arg = c(1:dim(data)[[2]]), cex.names = 0.7)

sampleTree <- hclust(dist(t(data)), method = "average")
dend <- as.dendrogram(sampleTree)
plot(dend)
cutHeight <- 155
abline(h = cutHeight, col = "red")

clust <- cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 1)
names(clust) <- sampleTree$labels
clust <- clust[labels(dend)]
table(clust)

keepSamples <- (clust[colnames(data)]==1)
data <- data[, keepSamples]
metadata_all <- subset(metadata_all, subset = SRR_id %in% colnames(data))
table(metadata_all$group)

gene_mad = apply(data, MARGIN = 1, FUN = mad)
hist(gene_mad)
gene_mad_rank = rank(-gene_mad, ties.method = "first")
keep = gene_mad_rank <= 3000
hist(gene_mad[keep])
data <- data[keep,]
dim(data)
gsg = goodSamplesGenes(data,verbose = 3)
gsg$allOK

data <- as.data.frame(t(data))
data[1:5,1:10]


powers = c(c(1:20))

sft = pickSoftThreshold(data, powerVector = powers, verbose = 5)

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

softpower <- 10
adjacency <- adjacency(data, power = softpower, type = "signed")
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1-TOM

k <- as.vector(apply(adjacency, 2, sum, na.rm=T))
hist(k)
scaleFreePlot(k, main = "Check scale free topology\n")


net <- blockwiseModules(data, power = softpower, TOMType = "signed",
                        minModuleSize = 20, reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, verbose = 3)

moduleNumbers <- net$colors
moduleColors <- labels2colors(net$colors)

plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)


MEs <- moduleEigengenes(data, colors = net$colors, excludeGrey = TRUE)$eigengenes

trait <- metadata_all[, c("SRR_id", "group", "batch")]
rownames(trait) <- trait[, 1]
trait <- trait[, 2:ncol(trait)]

trait$BD <- ifelse(trait$group == "BD", 1, 0)
trait$CTRL <- ifelse(trait$group == "CTRL", 1, 0)
trait$batch1 <- ifelse(trait$batch == "1", 1, 0)
trait$batch2 <- ifelse(trait$batch == "2", 1, 0)
trait$batch3 <- ifelse(trait$batch == "3", 1, 0)
trait$batch4 <- ifelse(trait$batch == "4", 1, 0)
trait <- trait[, 3:ncol(trait)]

moduleTraitCor <- cor(MEs, trait, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(trait))

corMat <- moduleTraitCor
pMat <- moduleTraitPvalue

mod_names <- rownames(corMat)
trait_names <- colnames(corMat)

cor_df <- melt(corMat)
p_df <- melt(pMat)
colnames(cor_df) <- c("Module", "Trait", "Correlation")
colnames(p_df) <- c("Module", "Trait", "Pvalue")

plot_df <- merge(cor_df, p_df, by = c("Module", "Trait"))
plot_df$label <- paste0("r = ", round(plot_df$Correlation, 2), 
                        "\np = ", formatC(plot_df$Pvalue, format = "e", digits = 2))

p <- ggplot(plot_df, aes(x = Trait, y = Module, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(low = "skyblue3", mid = "gray90", high = "gold3",
                       midpoint = 0, limit = c(-0.4, 0.4), name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Module-Trait Correlation")

print(p)
ggsave()


keep <- moduleColors != "grey"
genes <- colnames(data)[keep]
moduleColors <- moduleColors[keep]
moduleNumbers <- moduleNumbers[keep]

module_map_df <- data.frame(gene = genes, module_color  = moduleColors, module_number = moduleNumbers)
write.csv(module_map_df, "", row.names = FALSE)


library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)


module_map_df <- read.csv("")
module_map_df$module_number <- as.integer(module_map_df$module_number)

for (i in 1:length(unique(moduleNumbers))) {
  print(paste0("Start caculating: ME", i))
  
  genes <- subset(module_map_df, module_number == i)
  genes <- as.character(genes$gene)
  
  go_BP <- enrichGO(genes, universe = module_map_df$gene, OrgDb = "org.Hs.eg.db", 
                    keyType = "ENSEMBL", ont = "BP", readable = TRUE)
  
  if (!is.null(go_BP) && nrow(as.data.frame(go_BP)) > 0) {
    
    p <- dotplot(go_BP, showCategory = 10) + 
      labs(title = paste0("Enriched BP in ME", i)) + 
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    print(p)
    ggsave()
    
    p <- barplot(go_BP, showCategory = 10) + 
      labs(title = paste0("Enriched BP in ME", i)) + 
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    print(p)
    ggsave()
    
    p <- cnetplot(go_BP, node_label = "all")
    print(p)
    ggsave()
    
  } else {
    message("No significant GO terms for ME", i, " -> skip plotting")
  }
  
  go_MF <- enrichGO(genes, universe = module_map_df$gene, OrgDb = "org.Hs.eg.db", 
                    keyType = "ENSEMBL", ont = "MF", readable = TRUE)
  
  if (!is.null(go_MF) && nrow(as.data.frame(go_MF)) > 0) {
    
    p <- dotplot(go_MF, showCategory = 10) + 
      labs(title = paste0("Enriched MF in ME", i)) + 
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    print(p)
    ggsave()
    
    p <- barplot(go_MF, showCategory = 10) + 
      labs(title = paste0("Enriched MF in ME", i)) + 
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    print(p)
    ggsave()
    
  } else {
    message("No significant GO terms for ME", i, " -> skip plotting")
  }
}
