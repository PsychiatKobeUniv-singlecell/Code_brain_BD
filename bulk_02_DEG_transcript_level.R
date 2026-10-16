setwd("")

library(readr)
library(dplyr)
library(edgeR)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(GEOquery)


data <- read.csv("")

zero_rows <- apply(data[, 2:ncol(data)], 1, function(x) all(as.numeric(x) == 0))
sum(zero_rows)
data <- data[!zero_rows, ]
nrow(data)

rownames(data) <- data[, 1]
data <- data[, 2:ncol(data)]
data <- as.data.frame(data)


metadata_all <- read.csv("")

metadata_all$group <- factor(metadata_all$group, levels = c("CTRL", "BD"))
metadata_all$batch <- as.character(metadata_all$batch)
metadata_all$age <- as.numeric(metadata_all$age)
metadata_all$pmi <- as.numeric(metadata_all$pmi)
metadata_all$sex <- factor(metadata_all$sex, levels = c("M", "F"))

str(metadata_all)

data <- data[, metadata_all$SRR_id]
all(metadata_all$SRR_id == colnames(data))

mean_expression <- apply(data, 2, mean, na.rm = T)
barplot(mean_expression,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples",
        names.arg = c(1:dim(data)[[2]]), cex.names = 0.7)


dge <- DGEList(counts = data)

cpm_vals <- cpm(dge)
hist(log2(rowMeans(cpm_vals) + 1), main = "log2 average CPM", xlab = "log2(CPM + 1)")

keep <- rowSums(cpm(dge$counts) > 1) >= 0.8 * ncol(dge$counts)
table(keep)
dge <- dge[keep, , keep.lib.sizes = FALSE]

cpm_vals <- cpm(dge)
hist(log2(rowMeans(cpm_vals) + 1), main = "log2 average CPM", xlab = "log2(CPM + 1)")

dge <- calcNormFactors(dge, method = "TMM") 
dge$samples
dge_lcpm <- cpm(dge$counts, log = TRUE)

mean_expression <- apply(dge_lcpm, 2, mean, na.rm = T)
barplot(mean_expression,
        xlab = "Sample", ylab = "Mean expression (log2cpm)",
        main ="Mean expression across samples",
        names.arg = c(1:dim(dge_lcpm)[[2]]), cex.names = 0.7)

boxplot(dge_lcpm)

write.csv(dge_lcpm, "", row.names = TRUE)


design <- model.matrix(~ group + age + sex + pmi + batch, data = metadata_all)
design

dge <- estimateDisp(dge, design)

fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = "groupBD")
top_tag <- topTags(qlf, n = Inf, sort.by = "logFC")
head(top_tag)
top_tag$comparison

DEG_result <- top_tag$table
DEG_result$log10PValue <- -log10(DEG_result$PValue)
DEG_result$expression_change <- "Not Significant"

FC_threshold <- 1.2
FDR_threshold <- 0.05
PValue_threshold <- 0.05

DEG_result$expression_change[which(DEG_result$logFC >  log2(FC_threshold) & DEG_result$PValue < PValue_threshold)] <- "Up"
DEG_result$expression_change[which(DEG_result$logFC < -log2(FC_threshold) & DEG_result$PValue < PValue_threshold)] <- "Down"

DEG_result$expression_change[which(DEG_result$logFC >  log2(FC_threshold) & DEG_result$FDR < FDR_threshold)] <- "Up (FDR < 0.05)"
DEG_result$expression_change[which(DEG_result$logFC < -log2(FC_threshold) & DEG_result$FDR < FDR_threshold)] <- "Down (FDR < 0.05)"

write.csv(DEG_result, "", row.names = TRUE)


DEG_result <- read.csv("")
rownames(DEG_result) <- DEG_result[, 1]
DEG_result <- DEG_result[, 2:ncol(DEG_result)]
DEG_result <- as.data.frame(DEG_result)

tmp <- DEG_result[ , c("logFC", "log10PValue", "expression_change"), ]
tmp <- as.data.frame(tmp)
tmp <- tmp[which(!(is.na(tmp$logFC)|is.na(tmp$log10PValue))),]

tmp <- tmp %>% 
  rownames_to_column(var = "ensembl_transcript_id_version") %>% 
  mutate(expression_change = factor(expression_change, 
                                    levels = c("Up (FDR < 0.05)", "Down (FDR < 0.05)", "Up", "Down", "Not Significant")))

DEGs <- subset(tmp, subset = expression_change == "Up (FDR < 0.05)" | expression_change == "Down (FDR < 0.05)" | 
                 expression_change == "Up" | expression_change == "Down")


result <- read.csv("")
match_indices <- match(DEGs$ensembl_transcript_id_version, result$ensembl_transcript_id_version)

extract_matched_data <- function(index) {
  if (is.na(index)) {
    return(data.frame(ensembl_transcript_id_version = NA, hgnc_symbol = NA, stringsAsFactors = FALSE))
  } else {
    return(result[index, , drop = FALSE])
  }
}

table <- do.call(rbind, lapply(match_indices, extract_matched_data))

DEGs <- cbind(table, DEGs)

na_ensembl_idx <- is.na(DEGs[, 1])
DEGs[na_ensembl_idx, 1] <- DEGs[na_ensembl_idx, 3]
DEGs[na_ensembl_idx, 2] <- DEGs[na_ensembl_idx, 3]

empty_symbol_idx <- !is.na(DEGs[, 1]) & DEGs[, 2] == ""
DEGs[empty_symbol_idx, 2] <- DEGs[empty_symbol_idx, 3]

mismatch_idx <- DEGs[, 1] != DEGs[, 3]
if (any(mismatch_idx, na.rm = TRUE)) {
  warning("")
} else {
  DEGs <- DEGs[, -3]
}

rownames(DEGs) <- DEGs[, 1]
DEGs <- DEGs[, 2:ncol(DEGs)]

write.csv(DEGs, "", row.names = TRUE)
write.table(rownames(DEGs), "", row.names = FALSE, col.names = FALSE, sep = ",")


gene_map <- read.csv("", header = TRUE)
tmp <- merge(gene_map, tmp, by = "ensembl_transcript_id_version", all = TRUE)
tmp <- tmp[!is.na(tmp$logFC), ]
tmp[is.na(tmp)] <- ""
tmp <- tmp[order(tmp$log10PValue, decreasing = TRUE),]

p <- tmp %>% ggplot(aes(x = logFC, y = log10PValue, color = expression_change)) +
  geom_point(size = 1) +
  geom_text_repel(data = tmp[1:15, ], aes(label = hgnc_symbol),
                  size = 4, max.overlaps = Inf, show.legend = FALSE) +
  xlab("log2 (Fold Change)") +
  ylab("-log10 (pvalue)") +
  ggtitle("Transcript-level gene expression changes in BD vs CTRL") +
  scale_x_continuous(breaks = c(seq(-5, 5, 1)), limits = c(-5, 5)) +
  scale_y_continuous(breaks = c(seq(0, 14, 2)), limits = c(0, 14)) +
  scale_color_manual(name = "Expression Change", values = c(
    "Up (FDR < 0.05)" = "red3",
    "Down (FDR < 0.05)" = "blue3",
    "Up" = "lightpink2",
    "Down" = "skyblue1",
    "Not Significant" = "gray"
  )) +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.line = element_line(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)
        )

print(p)
ggsave()
