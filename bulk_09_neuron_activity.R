setwd("")

library(tidyverse)
library(ggplot2)
library(viridis)
library(readxl)


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

up_TF <- read_tsv("")
up_TF <- up_TF[1:100, ]$TF
down_TF <- read_tsv("")
down_TF <- down_TF[1:100, ]$TF

universe_genes <- unique(c(all_genes, all_transcripts))


read_neuron_activity_gene_sets <- function(
    excel_file,
    gene_col = "Gene"
) {
  
  sheet_names <- excel_sheets(excel_file)
  
  gene_sets <- lapply(sheet_names, function(sh) {
    
    df <- read_excel(excel_file, sheet = sh)
    
    if (!gene_col %in% colnames(df)) {
      stop(paste0("Column '", gene_col, "' is not found in sheet: ", sh))
    }
    
    genes <- df[[gene_col]] %>% as.character() %>% unique()
    
    return(genes)
    
  })
  
  names(gene_sets) <- sheet_names %>% str_replace_all(c(".45min" = "_ERG", ". 4hr" = "_LRG"))
  
  return(gene_sets)
}


run_neuron_activity_enrichment <- function(
    excel_file,
    up_genes,
    down_genes,
    up_TF,
    down_TF,
    genes_GWAS,
    gene_col = "Gene",
    universe_genes = NULL,
    min_overlap = 1,
    p_adjust_method = "BH"
) {
  
  neuron_activity_gene_sets <- read_neuron_activity_gene_sets(
    excel_file = excel_file,
    gene_col = gene_col
  )
  
  gene_sets <- list(
    "Upregulated genes" = unique(na.omit(as.character(up_genes))),
    "Downregulated genes" = unique(na.omit(as.character(down_genes))),
    "Up_DEG_enriched_TF" = unique(na.omit(as.character(up_TF))),
    "Down_DEG_enriched_TF" = unique(na.omit(as.character(down_TF))),
    "GWAS genes" = unique(na.omit(as.character(genes_GWAS)))
  )
  
  if (is.null(universe_genes)) {
    universe_genes <- unique(c(
      unlist(neuron_activity_gene_sets),
      unlist(gene_sets)
    ))
  } else {
    universe_genes <- unique(na.omit(as.character(universe_genes)))
  }
  
  neuron_activity_gene_sets <- lapply(neuron_activity_gene_sets, function(x) {
    unique(x[x %in% universe_genes])
  })
  
  gene_sets <- lapply(gene_sets, function(x) {
    unique(x[x %in% universe_genes])
  })
  
  N <- length(universe_genes)
  
  res <- purrr::map_dfr(names(gene_sets), function(set_name) {
    
    query_genes <- gene_sets[[set_name]]
    K <- length(query_genes)
    
    purrr::map_dfr(names(neuron_activity_gene_sets), function(neuron_activity_name) {
      
      target_genes <- neuron_activity_gene_sets[[neuron_activity_name]]
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
        neuron_activity = neuron_activity_name,
        overlap_count = x,
        neuron_activity_count = M,
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


plot_neuron_activity_enrichment_heatmap <- function(
    enrichment_result,
    use_padj = TRUE,
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
      neuron_activity,
      fill = list(
        overlap_count = 0,
        significance = fill_na
      )
    )
  
  plot_df <- plot_df %>%
    mutate(
      gene_set = factor(
        gene_set,
        levels = c("Upregulated genes", "Downregulated genes", 
                   "Up_DEG_enriched_TF", "Down_DEG_enriched_TF", 
                   "GWAS genes")
      )
    )
  
  neuron_activity_order <- plot_df %>%
    group_by(neuron_activity) %>%
    summarise(max_sig = max(significance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_sig)) %>%
    pull(neuron_activity)
  
  plot_df <- plot_df %>%
    mutate(neuron_activity = factor(neuron_activity, levels = neuron_activity_order))
  
  p <- ggplot(plot_df, aes(x = neuron_activity, y = gene_set, fill = significance)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_viridis_c(
      name = ifelse(use_padj, "-log10(adj. P)", "-log10(P)"),
      option = "viridis",
      na.value = "grey90"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank()
    )
  
  if (show_gene_count) {
    p <- p + geom_text(
      aes(label = ifelse(overlap_count > 0, overlap_count, "")),
      size = 4
      )
  }
  return(p)
}


enr <- run_neuron_activity_enrichment(
  excel_file = "",
  up_genes = up_genes,
  down_genes = down_genes,
  up_TF = up_TF,
  down_TF = down_TF,
  genes_GWAS = genes_GWAS,
  universe_genes = universe_genes,
  gene_col = "Gene"
)

write.csv(enr, "", row.names = FALSE)

p <- plot_neuron_activity_enrichment_heatmap(
  enrichment_result = enr,
  use_padj = TRUE,
  show_gene_count = TRUE
)
print(p)
ggsave()


activity_genes <- read_neuron_activity_gene_sets(
  excel_file = "",
  gene_col = "Gene"
) %>% unlist() %>% unique()

tf_file <- read_tsv("", show_col_types = FALSE)

selected_tf <- c("NPAS4", "NR4A1")

tf_targets <- tf_file %>%
  filter(TF %in% selected_tf) %>%
  mutate(gene = str_split(Overlapping_Genes, pattern = "\\s*,\\s*")) %>%
  select(TF, gene) %>%
  unnest(gene) %>%
  distinct(TF, gene) %>%
  group_by(TF) %>%
  summarise(genes = list(sort(unique(gene))), .groups = "drop") %>%
  deframe()

plot_tf_target_categories <- function(
    tf_targets,
    genes_GWAS,
    activity_genes,
    bar_width = 0.8,
    label_wrap_width = 20,
    label_size = 4,
    show_count = TRUE,
    tf_order = NULL
) {
  
  target_long <- imap_dfr(tf_targets, function(genes, tf_name) {
    tibble(TF = tf_name, gene = genes)
  }) %>% distinct(TF, gene)
  
  category_levels <- c("GWAS gene", "activity-dependent", "both", "none")
  classification_table <- target_long %>%
    mutate(
      is_gwas = gene %in% genes_GWAS,
      is_activity = gene %in% activity_genes,
      category = case_when(
        is_gwas & is_activity ~ "both",
        is_gwas ~ "GWAS gene",
        is_activity ~ "activity-dependent",
        TRUE ~ "none"
      ),
      category = factor(category, levels = category_levels)
    ) %>% select(TF, gene, category)
  
  gene_label_table <- classification_table %>%
    filter(category != "none") %>%
    group_by(TF, category) %>%
    summarise(genes = paste(unique(gene), collapse = ", "), .groups = "drop")
  
  plot_data <- classification_table %>%
    count(TF, category, name = "n") %>%
    complete(TF, 
             category = factor(category_levels, levels = category_levels),
             fill = list(n = 0)) %>%
    left_join(gene_label_table, by = c("TF", "category")) %>%
    mutate(category = factor(category, levels = category_levels))
  
  if (is.null(tf_order)) {
    
    tf_order <- names(tf_targets)
    
  } else {
    tf_order <- as.character(tf_order)
    
    missing_tf <- setdiff(unique(plot_data$TF), tf_order)
    
    if (length(missing_tf) > 0) {
      tf_order <- c(tf_order, missing_tf)
    }
  }
  
  plot_data <- plot_data %>%
    mutate(TF = factor(TF, levels = rev(tf_order)))
  
  classification_table <- classification_table %>%
    mutate(TF = factor(TF, levels = tf_order),
           category = factor(category, levels = category_levels)
           ) %>% 
    arrange(TF, category, gene)
  
  plot_data <- plot_data %>%
    mutate(
      gene_label = case_when(
        category == "none" ~ "",
        n == 0 ~ "",
        is.na(genes) ~ "",
        TRUE ~ str_wrap(genes, width = label_wrap_width)
      ),
      plot_label = case_when(
        category == "none" ~ "",
        n == 0 ~ "",
        show_count & gene_label != "" ~ paste0("n = ", n, "\n", gene_label),
        !show_count & gene_label != "" ~ gene_label,
        TRUE ~ ""
      )
    )
  
  p <- ggplot(plot_data, aes(x = TF, y = n, fill = category)) +
    geom_col(width = bar_width, color = "white", linewidth = 0.2) +
    geom_text(aes(label = plot_label), 
              position = position_stack(vjust = 0.5),
              size = label_size) +
    coord_flip(clip = "off") +
    scale_fill_manual(
      values = c(
        "GWAS gene" = "pink",
        "activity-dependent" = "lightblue1",
        "both" = "lightgoldenrod1",
        "none" = "gray80"
      ),
      breaks = category_levels, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Number of target genes", fill = "Category") +
    theme_bw() +
    theme(axis.text.y = element_text(face = "bold", size = 14), legend.position = "right")
  
  return(list(plot = p, classification = classification_table, summary = plot_data))
}

result <- plot_tf_target_categories(
  tf_targets = tf_targets,
  genes_GWAS = genes_GWAS,
  activity_genes = activity_genes,
  show_count = TRUE
)

p <- result$plot
print(p)
ggsave()
