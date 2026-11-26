setwd("")

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)


G_2018a <- read.csv("")
head(G_2018a)
up_DEG_array <- subset(G_2018a, beta > 0 & p < 0.05)
up_DEG_array <- unique(as.character(up_DEG_array$X))

down_DEG_array <- subset(G_2018a, beta < 0 & p < 0.05)
down_DEG_array <- unique(as.character(down_DEG_array$X))

G_2018b_DEG <- read.csv("")
head(G_2018b_DEG)
unique(G_2018b_DEG$Disorder.DGE_RegulationDirection)
G_2018b_DEG <- subset(G_2018b_DEG, Disorder.DGE_RegulationDirection == "BD.DGE_up" | Disorder.DGE_RegulationDirection == "BD.DGE_down")

up_DEG_RNAseq <- subset(G_2018b_DEG, Disorder.DGE_RegulationDirection == "BD.DGE_up")
up_DEG_RNAseq <- unique(as.character(up_DEG_RNAseq$Ensembl_Name))

down_DEG_RNAseq <- subset(G_2018b_DEG, Disorder.DGE_RegulationDirection == "BD.DGE_down")
down_DEG_RNAseq <- unique(as.character(down_DEG_RNAseq$Ensembl_Name))


G_2018b_DET <- read.csv("")
head(G_2018b_DET)
unique(G_2018b_DET$Disorder.DTE_RegulationDirection)
G_2018b_DET <- subset(G_2018b_DET, Disorder.DTE_RegulationDirection == "BD.DTE_up" | Disorder.DTE_RegulationDirection == "BD.DTE_down")

up_DET_RNAseq <- subset(G_2018b_DET, Disorder.DTE_RegulationDirection == "BD.DTE_up")
up_DET_RNAseq <- unique(as.character(up_DET_RNAseq$Ensembl_Transcript))

down_DET_RNAseq <- subset(G_2018b_DET, Disorder.DTE_RegulationDirection == "BD.DTE_down")
down_DET_RNAseq <- unique(as.character(down_DET_RNAseq$Ensembl_Transcript))


DEG_result <- read.csv("")
DET_result <- read.csv("")

up_DEG <- subset(DEG_result, expression_change == "Up (FDR < 0.05)" | expression_change == "Up")
up_DEG <- unique(as.character(up_DEG$X))
up_DET <- subset(DET_result, expression_change == "Up (FDR < 0.05)" | expression_change == "Up")
up_DET <- unique(as.character(up_DET$X))

down_DEG <- subset(DEG_result, expression_change == "Down (FDR < 0.05)" | expression_change == "Down")
down_DEG <- unique(as.character(down_DEG$X))
down_DET <- subset(DET_result, expression_change == "Down (FDR < 0.05)" | expression_change == "Down")
down_DET <- unique(as.character(down_DET$X))


## DEGs vs DEGs

row_sets <- list(
  Gandal_etal_RNAseq_DEG_up   = up_DEG_RNAseq,
  Gandal_etal_RNAseq_DEG_down = down_DEG_RNAseq,
  Gandal_etal_array_DEG_up    = up_DEG_array,
  Gandal_etal_array_DEG_down  = down_DEG_array
)

col_sets <- list(
  This_study_DEG_up   = up_DEG,
  This_study_DEG_down = down_DEG
)

clean_genes <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  
  x <- substr(x, 1, 15)
  
  unique(x[!is.na(x) & x != ""])
}

row_sets <- lapply(row_sets, clean_genes)
col_sets <- lapply(col_sets, clean_genes)

overlap_df <- expand.grid(
  Row = names(row_sets),
  Column = names(col_sets),
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    Overlapping_Genes = paste(
      sort(intersect(
        row_sets[[Row]],
        col_sets[[Column]]
      )),
      collapse = ";"
    ),
    Overlap_Count = length(
      intersect(
        row_sets[[Row]],
        col_sets[[Column]]
      )
    ),
    Column_Total = length(col_sets[[Column]]),
    Overlap_Percent = ifelse(
      Column_Total > 0,
      Overlap_Count / Column_Total * 100,
      NA_real_
    )
  ) %>%
  ungroup()

overlap_df <- overlap_df %>%
  mutate(
    Label = ifelse(
      is.na(Overlap_Percent),
      "NA",
      sprintf(
        "%d\n(%.1f%%)",
        Overlap_Count,
        Overlap_Percent
      )))

overlap_df$Row <- factor(overlap_df$Row, levels = rev(names(row_sets)))
overlap_df$Column <- factor(overlap_df$Column, levels = names(col_sets))

p <- ggplot(overlap_df, aes(x = Column, y = Row, fill = Overlap_Percent)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = Label),
    size = 3.5,
    color = "white",
    lineheight = 1.1
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(0, 50),
    name = "Overlap (%)",
  ) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right"
    )

print(p)
ggsave()

write.csv(overlap_df, "", row.names = FALSE)


## DETs vs DETs

row_sets <- list(
  Gandal_etal_RNAseq_DET_up   = up_DET_RNAseq,
  Gandal_etal_RNAseq_DET_down = down_DET_RNAseq
)

col_sets <- list(
  This_study_DET_up   = up_DET,
  This_study_DET_down = down_DET
)

clean_genes <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  
  x <- substr(x, 1, 15)
  
  unique(x[!is.na(x) & x != ""])
}

row_sets <- lapply(row_sets, clean_genes)
col_sets <- lapply(col_sets, clean_genes)

overlap_df <- expand.grid(
  Row = names(row_sets),
  Column = names(col_sets),
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    Overlapping_Genes = paste(
      sort(intersect(
        row_sets[[Row]],
        col_sets[[Column]]
      )),
      collapse = ";"
    ),
    Overlap_Count = length(
      intersect(
        row_sets[[Row]],
        col_sets[[Column]]
      )
    ),
    Column_Total = length(col_sets[[Column]]),
    Overlap_Percent = ifelse(
      Column_Total > 0,
      Overlap_Count / Column_Total * 100,
      NA_real_
    )
  ) %>%
  ungroup()

overlap_df <- overlap_df %>%
  mutate(
    Label = ifelse(
      is.na(Overlap_Percent),
      "NA",
      sprintf(
        "%d\n(%.1f%%)",
        Overlap_Count,
        Overlap_Percent
      )))

overlap_df$Row <- factor(overlap_df$Row, levels = rev(names(row_sets)))
overlap_df$Column <- factor(overlap_df$Column, levels = names(col_sets))

p <- ggplot(overlap_df, aes(x = Column, y = Row, fill = Overlap_Percent)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = Label),
    size = 3.5,
    color = "white",
    lineheight = 1.1
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(0, 50),
    name = "Overlap (%)",
  ) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right"
  )

print(p)
ggsave()

write.csv(overlap_df, "", row.names = FALSE)
