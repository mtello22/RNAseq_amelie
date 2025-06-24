library(data.table)
library(readr)
library(tximport)
library(tximportData)
library(tximeta)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(tidyverse)


# This will receive the stdin
exp_dirs <- list.dirs(path = "~/GitHub/RNAseq_amelie/data/salmon_quant", full.names = TRUE, 
                      recursive = FALSE)
exp_files <- list.files(path = exp_dirs, pattern = "quant.sf", full.names = TRUE)

# Extract sample names
exp_name <- exp_files %>% 
  str_replace(pattern = ".+/data/salmon_quant/", replacement = "") %>%
  str_replace(pattern = "/quant.sf", replacement = "") %>%
  str_replace(pattern = "-HAZ.+", replacement = "") 

# Code from https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# Takes the gtf file to make a table reference between ENST and ENSG
txdb <- makeTxDbFromGFF(file = "~/GitHub/RNAseq_amelie/data/genome_annot/Homo_sapiens.GRCh38.114.gtf.gz")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# Reads files and merge them 
txi <- tximport(files = exp_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# Generates metadata info
condition <- exp_name %>%
  str_replace("[1-3]", "")
deseq_coldat <- data.frame(condition = condition, row.names = exp_name, stringsAsFactors = TRUE)
deseq_coldat$condition <- relevel(deseq_coldat$condition, ref = "CTRL")


# Outputs estiamted counts summarized by ENSG
est_counts <- txi$counts
est_counts <- data.table(est_counts, keep.rownames = TRUE)
names(est_counts) <- c("ENSG", row.names(deseq_coldat))


# Outputs full txi and coldat objects
fwrite(est_counts, "~/GitHub/RNAseq_amelie/data/salmon_quant/Salmon_EstCount_ENSG.tsv", 
       quote = FALSE, sep = '\t', append = FALSE,
       row.names = FALSE, col.names = TRUE)
save(txi, deseq_coldat, file = "~/GitHub/RNAseq_amelie/data/salmon_quant/txi_coldat_object.rds")

