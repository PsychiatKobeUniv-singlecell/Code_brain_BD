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
filters[grep("transcript", filters[,1]),]

attributes <- c("ensembl_transcript_id_version", "hgnc_symbol")
filters <- "ensembl_transcript_id_version"

result <- getBM(attributes = attributes, filters = filters, values = rownames(data), mart = dataset, useCache = FALSE)
write.csv(result, "", row.names = FALSE)
