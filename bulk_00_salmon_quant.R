setwd("")

library(dplyr)
library(readr)
library(tximport)
library(txdbmaker)
library(GenomicFeatures)


txdb <- txdbmaker::makeTxDbFromGFF("", format = "gtf")

tx2gene <- select(txdb,
                  keys = keys(txdb, keytype = "TXNAME"),
                  columns = c("TXNAME", "GENEID"),
                  keytype = "TXNAME")

write.table(tx2gene, file = "", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)


tx2gene <- read.delim("", header = FALSE, stringsAsFactors = FALSE)
colnames(tx2gene) <- c("TXNAME", "GENEID")

sample_dirs <- list.dirs(path = c("./GSE53239_salmon", "./GSE78936_salmon", "./GSE80655_salmon", "./GSE112523_salmon"), full.names = TRUE, recursive = FALSE)
sample_names <- basename(sample_dirs)

files <- file.path(sample_dirs, "quant.sf")
names(files) <- sample_names
head(read_tsv(files[1]))


txi_tx <- tximport(files, type = "salmon", txOut = TRUE, ignoreAfterBar = TRUE)

txi_gene <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = TRUE)

write.csv(as.data.frame(txi_tx$abundance), "")
write.csv(as.data.frame(txi_tx$counts), "")

write.csv(as.data.frame(txi_gene$abundance), "")
write.csv(as.data.frame(txi_gene$counts), "")


txi_tx_tpm <- tximport(files, type = "salmon", txOut = TRUE, 
                         ignoreAfterBar = TRUE, countsFromAbundance = "lengthScaledTPM")

write.csv(as.data.frame(txi_tx_tpm$abundance), "")


txi_gene_tpm <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                         ignoreAfterBar = TRUE, countsFromAbundance = "lengthScaledTPM")

write.csv(as.data.frame(txi_gene_tpm$abundance), "")


sample_dirs <- list.dirs(path = c("./GSE80655_salmon", "./GSE112523_salmon"), full.names = TRUE, recursive = FALSE)
sample_names <- basename(sample_dirs)

files <- file.path(sample_dirs, "quant.sf")
names(files) <- sample_names
head(read_tsv(files[1]))


txi_tx <- tximport(files, type = "salmon", txOut = TRUE, ignoreAfterBar = TRUE)

txi_gene <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = TRUE)

write.csv(as.data.frame(txi_tx$abundance), "")
write.csv(as.data.frame(txi_tx$counts), "")

write.csv(as.data.frame(txi_gene$abundance), "")
write.csv(as.data.frame(txi_gene$counts), "")


txi_tx_tpm <- tximport(files, type = "salmon", txOut = TRUE, 
                       ignoreAfterBar = TRUE, countsFromAbundance = "lengthScaledTPM")

write.csv(as.data.frame(txi_tx_tpm$abundance), "")


txi_gene_tpm <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                         ignoreAfterBar = TRUE, countsFromAbundance = "lengthScaledTPM")

write.csv(as.data.frame(txi_gene_tpm$abundance), "")
