setwd("")

library(tidyverse)

plot_tf_overlap <- function(input_file, top_n = 10) {
  
  df <- read_tsv(input_file, show_col_types = FALSE)
  
  top_tf <- df %>%
    arrange(Rank) %>%
    slice_head(n = top_n) %>%
    mutate(TF = as.character(TF))
  
  plot_df <- top_tf %>%
    select(Rank, TF, Overlapping_Genes) %>%
    separate_rows(Overlapping_Genes, sep = "\\s*,\\s*") %>%
    rename(Gene = Overlapping_Genes) %>%
    mutate(Gene = str_trim(Gene)) %>%
    filter(Gene != "")
  
  all_genes <- sort(unique(plot_df$Gene))
  
  matrix_df <- expand_grid(
    TF = top_tf$TF,
    Gene = all_genes
  ) %>%
    left_join(
      plot_df %>% mutate(has_overlap = TRUE),
      by = c("TF", "Gene")
    ) %>%
    left_join(
      top_tf %>% select(TF, Rank),
      by = "TF"
    ) %>%
    mutate(
      has_overlap = replace_na(has_overlap, FALSE),
      mark = case_when(
        TF == Gene ~ "x",
        has_overlap ~ "dot",
        TRUE ~ "empty"
      ),
      TF = factor(TF, levels = rev(top_tf$TF)),
      Gene = factor(Gene, levels = all_genes)
    )
  
  p <- ggplot(matrix_df, aes(x = Gene, y = TF)) +
    geom_tile(fill = "white", color = "grey85") +
    geom_point(
      data = matrix_df %>% filter(mark == "dot"),
      size = 3
    ) +
    geom_text(
      data = matrix_df %>% filter(mark == "x"),
      aes(label = "x"),
      size = 4
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )
  
  return(p)
  
}


p <- plot_tf_overlap(input_file = "", top_n = 3)
print(p)
ggsave()


p <- plot_tf_overlap(input_file = "", top_n = 10)
print(p)
ggsave()
