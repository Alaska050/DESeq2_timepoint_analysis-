load("/Users/disha/Downloads/sim_counts_22.RData")

head(count_table)
is.matrix(count_table)
is.numeric(count_table)

pheno <- data.frame(tp = c("T1", "T1", "T1", "T24", "T24", "T24"))
deseq_dataset <- DESeqDataSetFromMatrix(countData = count_table, colData = pheno, design = ~tp)
colData(deseq_dataset)
deseq_dataset <- estimateSizeFactors(deseq_dataset)
sizeFactors(deseq_dataset)
counts(deseq_dataset, normalized=TRUE)
deseq_dataset <- estimateDispersions(deseq_dataset)
deseq_dataset <- nbinomWaldTest(deseq_dataset)
