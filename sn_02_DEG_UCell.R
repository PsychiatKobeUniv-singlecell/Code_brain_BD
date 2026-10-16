setwd("")

library(Matrix)
library(dplyr)
library(tidyr)
library(Seurat)
library(UCell)
library(SummarizedExperiment)
library(tidyverse)
library(ggplot2)
library(viridis)
library(patchwork)


data_sc <- readRDS("")

data_sc$developmental_stage <- str_replace_all(data_sc$stage, 
                                               c("ga22|ga24|ga34" = "prenatal",
                                                 "1m|2m|3m|4m|6m|10m" = "0y-1y",
                                                 "1y_a|1y_b|1y_c|1y_d|2y_a|^4y|^6y|^8y" = "1y-10y",
                                                 "10y_a|10y_b|12y|14y|16y|17y" = "10y-18y", 
                                                 "20y_a|20y_b|25y" = "adult"))
any(is.na(data_sc$developmental_stage))
data_sc$developmental_stage <- factor(data_sc$developmental_stage, 
                                      levels = c("prenatal", "0y-1y", "1y-10y", "10y-18y", "adult"))


deg <- read.csv("")

deg_up <- subset(deg, expression_change == "Up (FDR < 0.05)" | expression_change == "Up")
deg_up <- unique(as.character(deg_up$hgnc_symbol))

deg_down <- subset(deg, expression_change == "Down (FDR < 0.05)" | expression_change == "Down")
deg_down <- unique(as.character(deg_down$hgnc_symbol))

deg_up <- intersect(deg_up, rownames(data_sc))
length(deg_up)

deg_down <- intersect(deg_down, rownames(data_sc))
length(deg_down)

data_sc <- AddModuleScore_UCell(
  data_sc,
  features = list(DEG_up = deg_up, DEG_down = deg_down)
)

colnames(data_sc@meta.data)[grepl("DEG", colnames(data_sc@meta.data))]


vals <- FetchData(data_sc, vars = "DEG_up_UCell")[,1]
v_min <- suppressWarnings(min(vals, na.rm = TRUE))
v_max <- suppressWarnings(max(vals, na.rm = TRUE))

p_list <- FeaturePlot(data_sc, features = "DEG_up_UCell", split.by = "developmental_stage", 
                      pt.size = 0.5, order = TRUE, combine = FALSE)

p_list <- lapply(p_list, function(p) {
  p +
    scale_color_gradientn(colours = c("lavenderblush", "pink1", "orchid3", "purple4"),
                          limits = c(v_min, v_max)) +
    theme(legend.position = "right", 
          plot.title = element_blank(), 
          panel.border = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.line = element_blank()
    )
})

gene_plot <- (p_list[[1]] | p_list[[2]] | p_list[[3]] | p_list[[4]] | p_list[[5]]) + 
  plot_layout(heights = c(1, 1), guides = "collect") & 
  theme(legend.position = "bottom", legend.title = element_blank()) &
  guides(color = guide_colorbar(direction = "horizontal",
                                barwidth  = unit(50, "mm"), barheight = unit(5, "mm")))

p <- wrap_plots(gene_plot, nrow = 1) + plot_layout(guides = "keep") + 
  plot_annotation(title = "DEG UCell score - Upregulated", 
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
print(p)
ggsave()


vals <- FetchData(data_sc, vars = "DEG_down_UCell")[,1]
v_min <- suppressWarnings(min(vals, na.rm = TRUE))
v_max <- suppressWarnings(max(vals, na.rm = TRUE))

p_list <- FeaturePlot(data_sc, features = "DEG_down_UCell", split.by = "developmental_stage", 
                      pt.size = 0.5, order = TRUE, combine = FALSE)

p_list <- lapply(p_list, function(p) {
  p +
    scale_color_gradientn(colours = c("azure2", "palegreen3", "steelblue3", "blue4"),
                          limits = c(v_min, v_max)) +
    theme(legend.position = "right", 
          plot.title = element_blank(), 
          panel.border = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.line = element_blank()
    )
})

gene_plot <- (p_list[[1]] | p_list[[2]] | p_list[[3]] | p_list[[4]] | p_list[[5]]) + 
  plot_layout(heights = c(1, 1), guides = "collect") & 
  theme(legend.position = "bottom", legend.title = element_blank()) &
  guides(color = guide_colorbar(direction = "horizontal",
                                barwidth  = unit(50, "mm"), barheight = unit(5, "mm")))

p <- wrap_plots(gene_plot, nrow = 1) + plot_layout(guides = "keep") + 
  plot_annotation(title = "DEG UCell score - Downregulated", 
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
print(p)
ggsave()


data_sc$DEG_up_UCell_z <- as.numeric(scale(data_sc@meta.data[["DEG_up_UCell"]]))
data_sc$DEG_down_UCell_z <- as.numeric(scale(data_sc@meta.data[["DEG_down_UCell"]]))

get_sig_label <- function(padj) {
  case_when(
    padj < 0.001 ~ "***",
    padj < 0.01  ~ "**",
    padj < 0.05  ~ "*",
    TRUE ~ ""
  )
}

meta <- data_sc@meta.data %>%
  filter(
    !is.na(DEG_up_UCell_z),
    !is.na(DEG_down_UCell_z),
    !is.na(cell_type_1),
    !is.na(developmental_stage)
  )


celltype_stat <- meta %>%
  group_by(cell_type_1) %>%
  summarise(
    mean_z = mean(DEG_up_UCell_z),
    n_cells = n(),
    .groups = "drop"
  )

p_celltype <- sapply(celltype_stat$cell_type_1, function(ct) {
  idx <- meta$cell_type_1 == ct
  
  if (sum(idx) < 10) return(NA)
  
  wilcox.test(
    meta$DEG_up_UCell_z[idx],
    meta$DEG_up_UCell_z[!idx],
    alternative = "greater"
  )$p.value
})

celltype_stat$p_val <- p_celltype
celltype_stat$p_adj <- p.adjust(celltype_stat$p_val, method = "BH")
celltype_stat$sig <- get_sig_label(celltype_stat$p_adj)

celltype_stat <- celltype_stat %>%
  arrange(desc(mean_z)) %>%
  mutate(cell_type_1 = factor(cell_type_1, levels = cell_type_1))

p <- ggplot(celltype_stat, aes(x = cell_type_1, y = mean_z)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = sig),
    vjust = ifelse(celltype_stat$mean_z >= 0, -0.4, 1.2),
    size = 5
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "black") +
  theme_classic(base_size = 14) +
  labs(x = "Cell type", y = "Mean Z-scored UCell", 
       title = "DEG UCell score by cell type - Upregulated") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

print(p)
ggsave()


celltype_stat <- meta %>%
  group_by(cell_type_1) %>%
  summarise(
    mean_z = mean(DEG_down_UCell_z),
    n_cells = n(),
    .groups = "drop"
  )

p_celltype <- sapply(celltype_stat$cell_type_1, function(ct) {
  idx <- meta$cell_type_1 == ct
  
  if (sum(idx) < 10) return(NA)
  
  wilcox.test(
    meta$DEG_down_UCell_z[idx],
    meta$DEG_down_UCell_z[!idx],
    alternative = "greater"
  )$p.value
})

celltype_stat$p_val <- p_celltype
celltype_stat$p_adj <- p.adjust(celltype_stat$p_val, method = "BH")
celltype_stat$sig <- get_sig_label(celltype_stat$p_adj)

celltype_stat <- celltype_stat %>%
  arrange(desc(mean_z)) %>%
  mutate(cell_type_1 = factor(cell_type_1, levels = cell_type_1))

p <- ggplot(celltype_stat, aes(x = cell_type_1, y = mean_z)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = sig),
    vjust = ifelse(celltype_stat$mean_z >= 0, -0.4, 1.2),
    size = 5
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "black") +
  theme_classic(base_size = 14) +
  labs(x = "Cell type", y = "Mean Z-scored UCell", 
       title = "DEG UCell score by cell type - Downregulated") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

print(p)
ggsave()


stage_stat <- meta %>%
  group_by(developmental_stage) %>%
  summarise(
    mean_z = mean(DEG_up_UCell_z),
    n_cells = n(),
    .groups = "drop"
  )

p_stage <- sapply(stage_stat$developmental_stage, function(st) {
  idx <- meta$developmental_stage == st
  
  if (sum(idx) < 10) return(NA)
  
  wilcox.test(
    meta$DEG_up_UCell_z[idx],
    meta$DEG_up_UCell_z[!idx],
    alternative = "greater"
  )$p.value
})

stage_stat$p_val <- p_stage
stage_stat$p_adj <- p.adjust(stage_stat$p_val, method = "BH")
stage_stat$sig <- get_sig_label(stage_stat$p_adj)

stage_stat <- stage_stat %>%
  arrange(desc(mean_z)) %>%
  mutate(developmental_stage = factor(developmental_stage, levels = developmental_stage))

p <- ggplot(stage_stat, aes(x = developmental_stage, y = mean_z)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = sig),
    vjust = ifelse(stage_stat$mean_z >= 0, -0.4, 1.2),
    size = 5
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "black") +
  theme_classic(base_size = 14) +
  labs(x = "Developmental stage", y = "Mean Z-scored UCell", 
       title = "DEG UCell score by developmental stage - Upregulated") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

print(p)
ggsave()


stage_stat <- meta %>%
  group_by(developmental_stage) %>%
  summarise(
    mean_z = mean(DEG_down_UCell_z),
    n_cells = n(),
    .groups = "drop"
  )

p_stage <- sapply(stage_stat$developmental_stage, function(st) {
  idx <- meta$developmental_stage == st
  
  if (sum(idx) < 10) return(NA)
  
  wilcox.test(
    meta$DEG_down_UCell_z[idx],
    meta$DEG_down_UCell_z[!idx],
    alternative = "greater"
  )$p.value
})

stage_stat$p_val <- p_stage
stage_stat$p_adj <- p.adjust(stage_stat$p_val, method = "BH")
stage_stat$sig <- get_sig_label(stage_stat$p_adj)

stage_stat <- stage_stat %>%
  arrange(desc(mean_z)) %>%
  mutate(developmental_stage = factor(developmental_stage, levels = developmental_stage))

p <- ggplot(stage_stat, aes(x = developmental_stage, y = mean_z)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = sig),
    vjust = ifelse(stage_stat$mean_z >= 0, -0.4, 1.2),
    size = 5
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "black") +
  theme_classic(base_size = 14) +
  labs(x = "Developmental stage", y = "Mean Z-scored UCell", 
       title = "DEG UCell score by developmental stage - Downregulated") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

print(p)
ggsave()
