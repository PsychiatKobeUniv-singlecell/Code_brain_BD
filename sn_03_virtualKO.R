setwd("")

library(dplyr)
library(tidyverse)
library(Seurat)
library(scTenifoldKnk)

obj <- readRDS("")

obj$developmental_stage <- str_replace_all(obj$stage, 
                                           c("ga22|ga24|ga34" = "prenatal",
                                             "1m|2m|3m|4m|6m|10m" = "0y-1y",
                                             "1y_a|1y_b|1y_c|1y_d|2y_a|^4y|^6y|^8y" = "1y-10y",
                                             "10y_a|10y_b|12y|14y|16y|17y" = "10y-18y", 
                                             "20y_a|20y_b|25y" = "adult"))
any(is.na(obj$developmental_stage))
obj$developmental_stage <- factor(obj$developmental_stage, 
                                  levels = c("prenatal", "0y-1y", "1y-10y", "10y-18y", "adult"))

set.seed(123)
cell_counts <- obj@meta.data %>% count(developmental_stage, name = "n_cells")
cell_counts

cells_to_keep <- obj@meta.data %>% tibble::rownames_to_column("cell_id") %>%
  group_by(developmental_stage) %>% slice_sample(n = 1500) %>% pull(cell_id)
obj <- subset(obj, cells = cells_to_keep)
saveRDS(obj, "")


obj <- FindVariableFeatures(obj, nfeatures = 3000)
genes_use <- VariableFeatures(obj)

genes_use <- unique(c(genes_use, "NPAS4"))
length(genes_use)
# [1] 3001

stages <- c("prenatal", "0y-1y", "1y-10y", "10y-18y", "adult")

for (st in stages) {
  
  print(paste0("Processing: ", st))
  
  obj_sub <- subset(obj, developmental_stage == st)
  
  expr <- GetAssayData(obj_sub, assay = "SCT", slot = "data")
  expr <- expr[intersect(genes_use, rownames(expr)), ]
  
  result <- scTenifoldKnk(countMatrix = expr, gKO = "NPAS4", qc = FALSE, nCores = 16)
  res <- result$diffRegulation
  sig_genes <- subset(res, abs(Z) > 2)
  
  write.csv(sig_genes, , row.names = FALSE)
  
}


library(dplyr)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)


stages <- c("prenatal", "0y-1y", "1y-10y", "10y-18y", "adult")

for (st in stages) {
  
  sig_genes <- read.csv()
  sig_genes <- as.character(sig_genes$gene)
  
  go_BP <- enrichGO(sig_genes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
  
  if (!is.null(go_BP) && nrow(as.data.frame(go_BP)) > 0) {
    
    p <- dotplot(go_BP, showCategory = 5)
    print(p)
    ggsave()
    
    p <- cnetplot(go_BP, showCategory = 5, node_label = "all")
    print(p)
    ggsave()
    
  } else {
    message("No significant GO terms for ", st, " -> skip plotting")
  }

}


obj <- FindVariableFeatures(obj, nfeatures = 3000)
genes_use <- VariableFeatures(obj)

genes_use <- unique(c(genes_use, "NR4A1"))
length(genes_use)
# [1] 3000

stages <- c("prenatal", "0y-1y", "1y-10y", "10y-18y", "adult")

for (st in stages) {
  
  print(paste0("Processing: ", st))
  
  obj_sub <- subset(obj, developmental_stage == st)
  
  expr <- GetAssayData(obj_sub, assay = "SCT", slot = "data")
  expr <- expr[intersect(genes_use, rownames(expr)), ]
  
  result <- scTenifoldKnk(countMatrix = expr, gKO = "NR4A1", qc = FALSE, nCores = 16)
  res <- result$diffRegulation
  sig_genes <- subset(res, abs(Z) > 2)
  
  write.csv(sig_genes, , row.names = FALSE)
  
}


library(dplyr)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)


# 発達段階ごとに実行
stages <- c("prenatal", "0y-1y", "1y-10y", "10y-18y", "adult")

for (st in stages) {
  
  sig_genes <- read.csv()
  sig_genes <- as.character(sig_genes$gene)
  
  go_BP <- enrichGO(sig_genes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
  
  if (!is.null(go_BP) && nrow(as.data.frame(go_BP)) > 0) {
    
    p <- dotplot(go_BP, showCategory = 5)
    print(p)
    ggsave()
    
    p <- cnetplot(go_BP, showCategory = 5, node_label = "all")
    print(p)
    ggsave()
    
  } else {
    message("No significant GO terms for ", st, " -> skip plotting")
  }
  
}


library(tidyverse)
library(ggplot2)
library(ggrepel)
library(Seurat)


obj <- readRDS("")
mean_expression <- AverageExpression(obj, assay = "SCT", group.by = "developmental_stage")
mean_expression <- as.data.frame(mean_expression$SCT)
colnames(mean_expression) <- c("prenatal", "0y-1y", "1y-10y", "10y-18y", "adult")

stages <- c("prenatal", "0y-1y", "1y-10y", "10y-18y", "adult")

for (st in stages) {
  
  sig_genes <- read.csv()
  sig_genes <- sig_genes[-1,]
  sig_genes$abs.Z <- abs(sig_genes$Z)
  
  mean_expression_sub <- mean_expression %>% 
    dplyr::filter(rownames(mean_expression) %in% sig_genes$gene) %>% 
    dplyr::select(all_of(st)) %>% 
    dplyr::rename(mean_expression = paste0(st)) %>% 
    rownames_to_column("gene")
  
  sig_genes <- left_join(sig_genes, mean_expression_sub, by = "gene")
  sig_genes$log2MeanExpression <- log2(sig_genes$mean_expression)
  sig_genes <- subset(sig_genes, log2MeanExpression != "-Inf")
  
  p <- ggplot(sig_genes, aes(x = abs.Z, y = log2MeanExpression)) +
    geom_point(color = "deepskyblue") +
    geom_text_repel(data = sig_genes[1:5, ], aes(label = gene),
                    size = 4, max.overlaps = Inf, show.legend = FALSE) +
    scale_x_continuous(breaks = c(seq(2, 6, 1)), limits = c(2, 6)) +
    scale_y_continuous(breaks = c(seq(-4, 6, 2)), limits = c(-4, 6)) +
    ggtitle(paste0("NPAS4 virtual KO: ", st)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)
  
  print(p)
  ggsave()
  
}


stages <- c("prenatal", "0y-1y", "1y-10y", "10y-18y", "adult")

for (st in stages) {
  
  sig_genes <- read.csv()
  sig_genes <- sig_genes[-1,]
  sig_genes$abs.Z <- abs(sig_genes$Z)
  
  mean_expression_sub <- mean_expression %>% 
    dplyr::filter(rownames(mean_expression) %in% sig_genes$gene) %>% 
    dplyr::select(all_of(st)) %>% 
    dplyr::rename(mean_expression = paste0(st)) %>% 
    rownames_to_column("gene")
  
  sig_genes <- left_join(sig_genes, mean_expression_sub, by = "gene")
  sig_genes$log2MeanExpression <- log2(sig_genes$mean_expression)
  sig_genes <- subset(sig_genes, log2MeanExpression != "-Inf")
  
  p <- ggplot(sig_genes, aes(x = abs.Z, y = log2MeanExpression)) +
    geom_point(color = "limegreen") +
    geom_text_repel(data = sig_genes[1:5, ], aes(label = gene),
                    size = 4, max.overlaps = Inf, show.legend = FALSE) +
    scale_x_continuous(breaks = c(seq(2, 6, 1)), limits = c(2, 6)) +
    scale_y_continuous(breaks = c(seq(-4, 6, 2)), limits = c(-4, 6)) +
    ggtitle(paste0("NR4A1 virtual KO: ", st)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)
  
  print(p)
  ggsave()
  
}


library(tidyverse)
library(eulerr)


gene_overlap_test <- function(gene_set1, gene_set2, universe) {
  
  gene_set1 <- intersect(unique(gene_set1), universe)
  gene_set2 <- intersect(unique(gene_set2), universe)
  
  overlap_genes <- intersect(gene_set1, gene_set2)
  
  a <- length(overlap_genes)
  b <- length(setdiff(gene_set1, gene_set2))
  c <- length(setdiff(gene_set2, gene_set1))
  d <- length(universe) - a - b - c
  
  tbl <- matrix(
    c(a, b,
      c, d),
    nrow = 2,
    byrow = TRUE
  )
  
  fisher_res <- fisher.test(tbl, alternative = "greater")
  
  data.frame(
    GeneSet1_size = length(gene_set1),
    GeneSet2_size = length(gene_set2),
    Overlap = a,
    Odds_ratio = unname(fisher_res$estimate),
    P_value = fisher_res$p.value
  )
}


stages <- c("prenatal", "0y-1y", "1y-10y", "10y-18y", "adult")

for (st in stages) {
  
  NPAS4_genes <- read.csv()
  NPAS4_genes <- as.character(NPAS4_genes$gene)
  
  NR4A1_genes <- read.csv()
  NR4A1_genes <- as.character(NR4A1_genes$gene)
  
  venn_list <- euler(list(NPAS4_genes = NPAS4_genes, NR4A1_genes = NR4A1_genes))
  
  p <- plot(venn_list, quantities = TRUE, fill = c("deepskyblue", "limegreen"), alpha = 0.8, main = st)
  print(p)
  ggsave()
  
  universe <- unique(c(genes_use, "NPAS4", "NR4A1"))
  overlap <- gene_overlap_test(gene_set1 = NPAS4_genes, gene_set2 = NR4A1_genes, universe = universe)
  print(overlap)
}

# GeneSet1_size GeneSet2_size Overlap Odds_ratio      P_value
#            97            96      83   1246.234 3.20483e-132
# GeneSet1_size GeneSet2_size Overlap Odds_ratio      P_value
#            62            71      58   3341.159 2.71434e-104
# GeneSet1_size GeneSet2_size Overlap Odds_ratio       P_value
#            83            97      71   648.2492 6.227638e-108
# GeneSet1_size GeneSet2_size Overlap Odds_ratio      P_value
#            40            90      37   665.7985 8.929144e-57
# GeneSet1_size GeneSet2_size Overlap Odds_ratio       P_value
#            63            69      60   5496.175 7.770059e-112
