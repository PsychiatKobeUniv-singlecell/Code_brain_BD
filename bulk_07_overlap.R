setwd("")

library(tidyverse)
library(eulerr)


genes_GWAS <- read_tsv("")
genes_GWAS <- genes_GWAS$MAPPED_GENE %>%
  as.character() %>%
  str_split("\\s-\\s|,\\s|;\\s") %>%
  unlist() %>%
  str_trim() %>%
  na_if("") %>%
  na.omit() %>%
  unique()

DEGs <- read.csv("")
DEGs <- unique(as.character(DEGs$hgnc_symbol))

DETs <- read.csv("")
DETs <- unique(as.character(DETs$hgnc_symbol))

venn_list <- euler(list(DEGs = DEGs, DETs = DETs, genes_GWAS = genes_GWAS))

p <- plot(venn_list, quantities = TRUE, fill = c("lightskyblue", "plum", "wheat"), alpha = 0.8)
print(p)
ggsave()

Reduce(intersect, list(genes_GWAS, DEGs, DETs))
# [1] "FADS2"   "RGMA"    "PIM1"    "ARAP3"   "CHST3"   "ATP13A4" "EGFR"    "ZNF385D" "FBLN1"  
# [10] "TMEM179" "GFRA2"   "CABLES1" "EPB41"   "GNA12"   "ADAMTS2" "TPBG"    "SLC44A2" "DLGAP2" 
# [19] "GNA14"   "GRTP1"  


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

all_genes <- read.csv("")
all_genes <- unique(na.omit(as.character(all_genes$hgnc_symbol)))
all_genes <- all_genes[all_genes != ""]

all_transcripts <- read.csv("")
all_transcripts <- unique(na.omit(as.character(all_transcripts$hgnc_symbol)))
all_transcripts <- all_transcripts[all_transcripts != ""]

universe <- unique(c(all_genes, all_transcripts, genes_GWAS))

overlap_DEG_GWAS <- gene_overlap_test(gene_set1 = DEGs, gene_set2 = genes_GWAS, universe = universe)
overlap_DEG_GWAS
# GeneSet1_size GeneSet2_size Overlap Odds_ratio   P_value
#          1166          1558      58   1.341831 0.0223552

overlap_DET_GWAS <- gene_overlap_test(gene_set1 = DETs, gene_set2 = genes_GWAS, universe = universe)
overlap_DET_GWAS
# GeneSet1_size GeneSet2_size Overlap Odds_ratio      P_value
#           921          1558      67   2.036555 3.528979e-07

overlap_DEG_DET <- gene_overlap_test(gene_set1 = DEGs, gene_set2 = DETs, universe = universe)
overlap_DEG_DET
# GeneSet1_size GeneSet2_size Overlap Odds_ratio P_value
#          1166           921     381   35.40852       0
