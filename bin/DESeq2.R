library(DESeq2)
library(data.table)
library(biomaRt)

load(file = "~/GitHub/RNAseq_amelie/data/salmon_quant/txi_coldat_object.rds")

DESeqDataset <- DESeqDataSetFromTximport(txi = txi, colData = deseq_coldat, design = ~ condition)

# Filter lowly expressed genes
min_reads <- 4
min_samples <- 10
index <- rowSums(counts(DESeqDataset) >= min_reads) > min_samples
DESeqDataset <- DESeqDataset[index, ]

# Perform DESeq2
DESeqDataset <- DESeq(DESeqDataset)

# resultsNames(DESeqDataset)
# "condition_FA_vs_CTRL" "condition_PNS_vs_CTRL" "condition_PTS_vs_CTRL"

DEG_results <- results(DESeqDataset, 
                        name="condition_FA_vs_CTRL",
                        pAdjustMethod = "BH", 
                        alpha = 0.05)

DEG_results <- as.data.table(keep.rownames = TRUE, DEG_results)
setnames(DEG_results, "rn", "ENSG")
DEG_results[, Group := "FA_vs_Ctrl"]

####################

deseq_coldat$condition <- relevel(deseq_coldat$condition, ref = "FA")
DESeqDataset <- DESeqDataSetFromTximport(txi = txi, colData = deseq_coldat, design = ~ condition)
# Filter lowly expressed genes
DESeqDataset <- DESeqDataset[index, ]

# Perform DESeq2
DESeqDataset <- DESeq(DESeqDataset)

# resultsNames(DESeqDataset)
# "condition_CTRL_vs_FA" "condition_PNS_vs_FA"  "condition_PTS_vs_FA"

DEGs <- results(DESeqDataset, 
                name="condition_PNS_vs_FA",
                pAdjustMethod = "BH", 
                alpha = 0.05)
DEGs <- as.data.table(DEGs, keep.rownames = TRUE)
setnames(DEGs, "rn", "ENSG")
DEGs[, Group := "PNS_vs_FA"]
DEG_results <- rbind(DEG_results, DEGs)

DEGs <- results(DESeqDataset, 
                name="condition_PTS_vs_FA",
                pAdjustMethod = "BH", 
                alpha = 0.05)
DEGs <- as.data.table(DEGs, keep.rownames = TRUE)
setnames(DEGs, "rn", "ENSG")
DEGs[, Group := "PTS_vs_FA"]
DEG_results <- rbind(DEG_results, DEGs)



## ENSEMBL Homo_sapiens.GRCh38.114.gtf.gz
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "hsapiens_gene_ensembl", 
                      version = 114)
gene_IDs <- getBM(filters= "ensembl_gene_id", 
                  attributes= c("ensembl_gene_id",
                                "ensembl_gene_id_version",
                                "external_gene_name", 
                                "entrezgene_id", 
                                "entrezgene_accession"),
                  values = unique(DEG_results$ENSG), 
                  mart= ensembl)
gene_IDs <- as.data.table(na.omit(gene_IDs))


DEG_results <- merge.data.table(y = DEG_results, 
                                x= gene_IDs[, .SD, 
                                            .SDcols = c("ensembl_gene_id",
                                                        "external_gene_name",
                                                        "entrezgene_id")], 
                                by.y = "ENSG", 
                                by.x = "ensembl_gene_id", 
                                all.x = FALSE)

fwrite(DEG_results, file = "~/GitHub/RNAseq_amelie/output/DEG_results.tsv", 
       quote = FALSE, append = FALSE, sep = '\t', 
       row.names = FALSE, col.names = TRUE)
