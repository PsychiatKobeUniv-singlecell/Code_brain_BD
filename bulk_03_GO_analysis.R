setwd("")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)


DEG_result <- read.csv("")
rownames(DEG_result) <- substr(DEG_result[, 1], 1, 15)
DEG_result <- DEG_result[, 2:ncol(DEG_result)]
DEG_result <- as.data.frame(DEG_result)

table(DEG_result$expression_change)

DEG_result <- DEG_result[order(DEG_result$logFC, decreasing = TRUE),]

geneList <- DEG_result[, 1]
names(geneList) <- rownames(DEG_result)
str(geneList)

DEGs_up <- subset(DEG_result, expression_change == "Up (FDR < 0.05)" | expression_change == "Up")
DEGs_up <- as.character(rownames(DEGs_up))

DEGs_down <- subset(DEG_result, expression_change == "Down (FDR < 0.05)" | expression_change == "Down")
DEGs_down <- as.character(rownames(DEGs_down))


go_BP <- enrichGO(DEGs_up, universe = names(geneList), OrgDb = "org.Hs.eg.db", 
                  keyType = "ENSEMBL", ont = "BP", readable = TRUE)

p <- dotplot(go_BP, showCategory = 10)
print(p)
ggsave()

p <- barplot(go_BP, showCategory = 10)
print(p)
ggsave()

p <- cnetplot(go_BP, foldChange = geneList, node_label = "all") + 
  scale_color_gradient2(low = "blue3", mid = "white", high = "red3", midpoint = 0)
print(p)
ggsave()

go_BP <- enrichGO(DEGs_down, universe = names(geneList), OrgDb = "org.Hs.eg.db", 
                  keyType = "ENSEMBL", ont = "BP", readable = TRUE)

p <- dotplot(go_BP, showCategory = 10)
print(p)
ggsave()

p <- barplot(go_BP, showCategory = 10)
print(p)
ggsave()

p <- cnetplot(go_BP, foldChange = geneList, node_label = "all") + 
  scale_color_gradient2(low = "blue3", mid = "white", high = "red3", midpoint = 0)
print(p)
ggsave()


go_MF <- enrichGO(DEGs_up, universe = names(geneList), OrgDb = "org.Hs.eg.db", 
                  keyType = "ENSEMBL", ont = "MF", readable = TRUE)

p <- dotplot(go_MF, showCategory = 10)
print(p)
ggsave()

p <- barplot(go_MF, showCategory = 10)
print(p)
ggsave()

p <- cnetplot(go_MF, foldChange = geneList, node_label = "all") + 
  scale_color_gradient2(low = "blue3", mid = "white", high = "red3", midpoint = 0)
print(p)
ggsave()

go_MF <- enrichGO(DEGs_down, universe = names(geneList), OrgDb = "org.Hs.eg.db", 
                  keyType = "ENSEMBL", ont = "MF", readable = TRUE)

p <- dotplot(go_MF, showCategory = 10)
print(p)
ggsave()

p <- barplot(go_MF, showCategory = 10)
print(p)
ggsave()

p <- cnetplot(go_MF, foldChange = geneList, node_label = "all") + 
  scale_color_gradient2(low = "blue3", mid = "white", high = "red3", midpoint = 0)
print(p)
ggsave()


set.seed(123)

gse_BP <- gseGO(geneList = geneList, ont = "BP", keyType = "ENSEMBL", seed = TRUE, 
                minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, OrgDb = org.Hs.eg.db)

p <- dotplot(gse_BP, showCategory = 10, split = ".sign") + facet_grid(.~.sign)
print(p)
ggsave()

p <- ridgeplot(gse_BP, showCategory = 10) + labs(x = "enrichment distribution")
print(p)
ggsave()

gse_MF <- gseGO(geneList = geneList, ont = "MF", keyType = "ENSEMBL", seed = TRUE, 
                minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, OrgDb = org.Hs.eg.db)

p <- dotplot(gse_MF, showCategory = 10, split = ".sign") + facet_grid(.~.sign)
print(p)
ggsave()

p <- ridgeplot(gse_MF, showCategory = 10) + labs(x = "enrichment distribution")
print(p)
ggsave()
