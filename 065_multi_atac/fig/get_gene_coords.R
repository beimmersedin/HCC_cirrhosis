
library(EnsDb.Hsapiens.v86)

genes <- read.csv("/data1/project/yeonu/065_multi_atac/fig/target_genes_list.csv", header=TRUE)$gene
gene_data <- genes(EnsDb.Hsapiens.v86, filter=GeneNameFilter(genes), columns=c("gene_name","seq_name","gene_seq_start","gene_seq_end","gene_biotype","seq_strand"))
gene_df <- as.data.frame(gene_data)
gene_df$chrom <- paste0("chr", gene_df$seq_name)
write.csv(gene_df, "/data1/project/yeonu/065_multi_atac/fig/target_gene_coords.csv", row.names=FALSE)
cat("Done:", nrow(gene_df), "genes
")
