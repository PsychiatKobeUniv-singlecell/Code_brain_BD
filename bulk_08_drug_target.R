setwd("")

library(tidyverse)
library(ggplot2)
library(viridis)


genes_GWAS <- read_tsv("")
genes_GWAS <- genes_GWAS$MAPPED_GENE %>%
  as.character() %>%
  str_split("\\s-\\s|,\\s|;\\s") %>%
  unlist() %>%
  str_trim() %>%
  na_if("") %>%
  na.omit() %>%
  unique()

DEG_result <- read.csv("")
DET_result <- read.csv("")

up_DEG <- subset(DEG_result, expression_change == "Up (FDR < 0.05)" | expression_change == "Up")
up_DEG <- unique(as.character(up_DEG$hgnc_symbol))
up_DET <- subset(DET_result, expression_change == "Up (FDR < 0.05)" | expression_change == "Up")
up_DET <- unique(as.character(up_DET$hgnc_symbol))

down_DEG <- subset(DEG_result, expression_change == "Down (FDR < 0.05)" | expression_change == "Down")
down_DEG <- unique(as.character(down_DEG$hgnc_symbol))
down_DET <- subset(DET_result, expression_change == "Down (FDR < 0.05)" | expression_change == "Down")
down_DET <- unique(as.character(down_DET$hgnc_symbol))

up_genes <- unique(c(up_DEG, up_DET))
down_genes <- unique(c(down_DEG, down_DET))

all_genes <- read.csv("")
all_genes <- unique(na.omit(as.character(all_genes$hgnc_symbol)))
all_genes <- all_genes[all_genes != ""]

all_transcripts <- read.csv("")
all_transcripts <- unique(na.omit(as.character(all_transcripts$hgnc_symbol)))
all_transcripts <- all_transcripts[all_transcripts != ""]

universe_genes <- unique(c(all_genes, all_transcripts))


run_drug_enrichment <- function(
    drug_target_file,
    up_genes,
    down_genes,
    genes_GWAS,
    drug_col = "drug",
    gene_col = "gene",
    universe_genes = NULL,
    min_overlap = 1,
    p_adjust_method = "BH"
) {
  
  drug_target <- read_tsv(drug_target_file)
  
  if (!all(c(drug_col, gene_col) %in% colnames(drug_target))) {
    stop("drug_col or gene_col is not found in drug_target_file.")
  }
  
  drug_target <- drug_target %>%
    dplyr::select(
      drug = all_of(drug_col),
      gene = all_of(gene_col)
    ) %>%
    filter(!is.na(drug), !is.na(gene)) %>%
    mutate(
      drug = as.character(drug),
      gene = as.character(gene)
    ) %>%
    distinct()
  
  drug_target$drug <- str_replace_all(drug_target$drug, c(" LAUROXIL" = "",
                                                          " HYDROCHLORIDE" = "",
                                                          " FUMARATE" = "",
                                                          " DECANOATE" = ""))
  
  gene_sets <- list(
    "Upregulated genes" = unique(na.omit(as.character(up_genes))),
    "Downregulated genes" = unique(na.omit(as.character(down_genes))),
    "GWAS genes" = unique(na.omit(as.character(genes_GWAS)))
  )
  
  if (is.null(universe_genes)) {
    universe_genes <- unique(c(
      drug_target$gene,
      unlist(gene_sets)
    ))
  } else {
    universe_genes <- unique(na.omit(as.character(universe_genes)))
  }
  
  drug_target <- drug_target %>%
    filter(gene %in% universe_genes)
  
  gene_sets <- lapply(gene_sets, function(x) {
    unique(x[x %in% universe_genes])
  })
  
  N <- length(universe_genes)
  
  drug_gene_list <- split(drug_target$gene, drug_target$drug)
  drug_gene_list <- lapply(drug_gene_list, unique)
  
  res <- purrr::map_dfr(names(gene_sets), function(set_name) {
    
    query_genes <- gene_sets[[set_name]]
    K <- length(query_genes)
    
    purrr::map_dfr(names(drug_gene_list), function(drug_name) {
      
      target_genes <- drug_gene_list[[drug_name]]
      M <- length(target_genes)
      
      overlap_genes <- intersect(query_genes, target_genes)
      x <- length(overlap_genes)
      
      mat <- matrix(
        c(
          x,
          M - x,
          K - x,
          N - M - K + x
        ),
        nrow = 2,
        byrow = TRUE
      )
      
      pval <- fisher.test(mat, alternative = "greater")$p.value
      
      tibble(
        gene_set = set_name,
        drug = drug_name,
        overlap_count = x,
        drug_target_count = M,
        gene_set_count = K,
        universe_count = N,
        pvalue = pval,
        overlap_genes = paste(overlap_genes, collapse = ";")
      )
    })
  }) %>%
    group_by(gene_set) %>%
    mutate(
      padj = p.adjust(pvalue, method = p_adjust_method),
      minus_log10_padj = -log10(padj)
    ) %>%
    ungroup() %>%
    filter(overlap_count >= min_overlap) %>%
    arrange(gene_set, padj, desc(overlap_count))
  
  return(res)
}


plot_drug_enrichment_dotplot <- function(
    enrichment_result,
    use_padj = TRUE,
    p_cutoff = 0.05,
    dot_size_range = c(2, 10)
) {
  
  plot_df <- enrichment_result
  
  if (use_padj) {
    plot_df <- plot_df %>%
      mutate(
        significance = -log10(padj),
        p_label = padj
      )
  } else {
    plot_df <- plot_df %>%
      mutate(
        significance = -log10(pvalue),
        p_label = pvalue
      )
  }
  
  plot_df <- plot_df %>%
    mutate(
      gene_set = factor(
        gene_set,
        levels = c("Upregulated genes", "Downregulated genes", "GWAS genes")
      ),
      drug = factor(drug, levels = c("LITHIUM", "VALPROIC ACID", "CARBAMAZEPINE", "LAMOTRIGINE", 
                                     "HALOPERIDOL", "OLANZAPINE", "QUETIAPINE", "ASENAPINE", 
                                     "ARIPIPRAZOLE", "PALIPERIDONE", "RISPERIDONE", "LURASIDONE"))
    )
  
  ggplot(plot_df, aes(x = drug, y = gene_set)) +
    geom_point(
      aes(size = overlap_count, color = significance)
    ) +
    scale_size_continuous(
      name = "Gene count",
      range = dot_size_range
    ) +
    scale_color_gradient(
      name = ifelse(use_padj, "-log10(adj. P)", "-log10(P)"),
      low = "grey80",
      high = "red3"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
}


plot_drug_enrichment_heatmap <- function(
    enrichment_result,
    use_padj = TRUE,
    top_n_drugs = NULL,
    fill_na = 0,
    show_gene_count = TRUE
) {
  
  plot_df <- enrichment_result
  
  if (use_padj) {
    plot_df <- plot_df %>%
      mutate(
        significance = -log10(padj),
        p_label = padj
      )
  } else {
    plot_df <- plot_df %>%
      mutate(
        significance = -log10(pvalue),
        p_label = pvalue
      )
  }
  
  plot_df <- plot_df %>%
    mutate(
      significance = ifelse(is.infinite(significance), NA, significance)
    )
  
  plot_df <- plot_df %>%
    complete(
      gene_set,
      drug,
      fill = list(
        overlap_count = 0,
        significance = fill_na
      )
    )
  
  plot_df <- plot_df %>%
    mutate(
      gene_set = factor(
        gene_set,
        levels = c("Upregulated genes", "Downregulated genes", "GWAS genes")
      )
    )
  
  drug_order <- plot_df %>%
    group_by(drug) %>%
    summarise(max_sig = max(significance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_sig)) %>%
    pull(drug)
  
  plot_df <- plot_df %>%
    mutate(
      drug = factor(drug, levels = drug_order)
    )
  
  p <- ggplot(plot_df, aes(x = drug, y = gene_set, fill = significance)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_viridis_c(
      name = ifelse(use_padj, "-log10(adj. P)", "-log10(P)"),
      option = "viridis",
      na.value = "grey90"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank()
    )
  
  if (show_gene_count) {
    p <- p +
      geom_text(
        aes(label = ifelse(overlap_count > 0, overlap_count, "")),
        size = 4
      )
  }
  
  return(p)
}


enr <- run_drug_enrichment(
  drug_target_file = "",
  up_genes = up_genes,
  down_genes = down_genes,
  genes_GWAS = genes_GWAS,
  universe_genes = universe_genes,
  drug_col = "drug",
  gene_col = "gene"
)

write.csv(enr, "", row.names = FALSE)

p <- plot_drug_enrichment_dotplot(
  enrichment_result = enr,
  use_padj = TRUE,
  )
print(p)
ggsave()

p <- plot_drug_enrichment_heatmap(
  enrichment_result = enr,
  use_padj = TRUE,
  show_gene_count = TRUE
)
print(p)
ggsave()
