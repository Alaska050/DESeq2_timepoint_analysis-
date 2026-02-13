#Differentail Gene Expression analysis using DESeq2

#Load count data
load("/Users/disha/Downloads/sim_counts_22.RData")

head(count_table)
is.matrix(count_table)
is.numeric(count_table)
dim(count_table)

library(DESeq2)
#Creating the DSeq2 object
pheno <- data.frame(tp = c("T1", "T1", "T1", "T24", "T24", "T24"))
deseq_dataset <- DESeqDataSetFromMatrix(countData = count_table, colData = pheno, design = ~tp)
colData(deseq_dataset)

#Normalizing library size differences
deseq_dataset <- estimateSizeFactors(deseq_dataset)
sizeFactors(deseq_dataset)
counts(deseq_dataset, normalized=TRUE)

#Seperating biological noise with techinal noise 
deseq_dataset <- estimateDispersions(deseq_dataset)
plotDispEsts(deseq_dataset)

#Differential expression testing
deseq_dataset <- nbinomWaldTest(deseq_dataset)
results_table <- results (deseq_dataset, contrast = c("tp", "T24", "T1"))
summary(results_table)
results_table <- results_table[order(results_table$padj),]
plotMA(deseq_dataset, alpha=0.01)
plotDispEsts(deseq_dataset)

#Plasmid gene exppression
plasmid_ids <- paste0("PROKKA_00", 950:957)
plasmid_res <- results_table[rownames(results_table) %in% plasmid_ids, ]
plasmid_res
plasmid_res[order(plasmid_res$padj), ]
counts(deseq_dataset, normalized=TRUE)[plasmid_ids, ]


#Data Transformation
rlog_data <- rlogTransformation(deseq_dataset, blind=TRUE)
dist_rl = dist(t(assay(rlog_data)))
dist_mat = as.matrix (dist_rl)
heatmap.2(dist_mat, trace = "none")
png("Sample_distance_heatmap.png", width=900, height=800)
heatmap.2(dist_mat, trace="none")
dev.off()
graphics.off()


plotPCA(rlog_data, intgroup="tp")
png("PCA_plot.png", width=900, height=700)
plotPCA(rlog_data, intgroup="tp")
dev.off()

# Load annotation map
ann_name <- load("/Users/disha/Downloads/annotation_map_22.RData")
annotation_table <- as.data.frame(get(ann_name), stringsAsFactors = FALSE)

# DESeq2 results to df
res_df <- as.data.frame(results_table)
res_df$prokka_id <- trimws(rownames(res_df))

# Create matching key column in annotation table
annotation_table$prokka_id <- trimws(as.character(annotation_table$PROKKA_ID))

# Merge
res_annot <- merge(res_df, annotation_table, by = "prokka_id", all.x = TRUE)

# Filter significant
sig_hits <- subset(res_annot, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)

# Check
nrow(sig_hits)
head(sig_hits)
sum(is.na(res_annot$NCBI_ID))
sum(duplicated(res_annot$prokka_id))

top10 <- res_annot[order(res_annot$padj), 
                   c("prokka_id","NCBI_Gene","NCBI_Description",
                     "log2FoldChange","padj","baseMean")]

top10 <- top10[1:10, ]
top10

plasmid_ids <- paste0("PROKKA_00", 950:957)

plasmid_table <- res_annot[res_annot$prokka_id %in% plasmid_ids,
                           c("prokka_id","NCBI_Gene","NCBI_Description",
                             "baseMean","log2FoldChange","padj")]

plasmid_table <- plasmid_table[order(plasmid_table$padj), ]
plasmid_table
write.csv(plasmid_table, "Plasmid_DESeq2_results.csv", row.names = FALSE)

norm_plasmid <- counts(deseq_dataset, normalized = TRUE)[plasmid_ids, ]

#Plamid expression barplot
png("Plasmid_normalised_counts.png", width=1100, height=700)
barplot(as.matrix(norm_plasmid),
        beside=TRUE, las=2,
        ylab="Normalised read counts",
        main="Plasmid gene expression across samples (DESeq2 normalised)")
legend("topright", legend=rownames(norm_plasmid), cex=0.65)
dev.off()

