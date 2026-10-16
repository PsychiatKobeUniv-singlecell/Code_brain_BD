setwd("")

library(readr)
library(dplyr)
library(biomaRt)


data <- read.csv("")

rownames(data) <- data[, 1]
data <- data[, 2:ncol(data)]
data <- as.data.frame(data)

ncol(data)
nrow(data)

database <- useMart("ensembl")
listDatasets(database)

dataset <- useDataset("hsapiens_gene_ensembl", mart = database)
dataset

filters <- listFilters(dataset)
filters[grep("ensembl_gene_id", filters[,1]),]

attributes <- c("ensembl_gene_id_version", "hgnc_symbol")
filters <- "ensembl_gene_id_version"

result <- getBM(attributes = attributes, filters = filters, values = rownames(data), mart = dataset, useCache = FALSE)
write.csv(result, "", row.names = FALSE)


match_indices <- match(rownames(data), result$ensembl_gene_id_version)

extract_matched_data <- function(index) {
  if (is.na(index)) {
    return(data.frame(ensembl_gene_id_version = NA, hgnc_symbol = NA, stringsAsFactors = FALSE))
  } else {
    return(result[index, , drop = FALSE])
  }
}

table <- do.call(rbind, lapply(match_indices, extract_matched_data))

data <- cbind(table, data)

empty_symbol_idx <- !is.na(data[, 1]) & data[, 2] == ""
data[empty_symbol_idx, 2] <- data[empty_symbol_idx, 1]

ncol(data)
nrow(data)

zero_rows <- apply(data[, 3:ncol(data)], 1, function(x) all(as.numeric(x) == 0))
sum(zero_rows)
data <- data[!zero_rows, ]
nrow(data)

write.csv(data, "", row.names = FALSE)


setwd("")

library(readr)
library(dplyr)
library(tidyverse)


data <- read.csv("")

zero_rows <- apply(data[, 2:ncol(data)], 1, function(x) all(as.numeric(x) == 0))
sum(zero_rows)
data <- data[!zero_rows, ]
nrow(data)

rownames(data) <- data[, 1]
data <- data[, 2:ncol(data)]
data <- as.data.frame(data)

gene_map <- read.csv("", header = TRUE)

data$ensembl_gene_id_version <- rownames(data)
merged <- merge(gene_map, data, by = "ensembl_gene_id_version")
merged <- merged[!merged$hgnc_symbol == "", ]
merged <- merged %>% group_by(hgnc_symbol) %>% summarise(across(where(is.numeric), sum))
merged <- column_to_rownames(as.data.frame(merged), var = "hgnc_symbol")

ncol(merged)
nrow(merged)

write.csv(merged, "")
