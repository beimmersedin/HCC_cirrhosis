# 라이브러리 로드
library(Seurat)
library(MAST)
library(Matrix)

# 출력 경로
out_dir <- "Deg/results_output"

# 1. 데이터 불러오기
# MTX가 cells x genes 형식(Python/AnnData 기준)이므로 transpose 필요
cat("Loading data...\n")
counts_raw <- readMM("data/counts_sparse.mtx")
cell_meta <- read.csv("data/cell_metadata.csv", row.names = 1)
gene_meta <- read.csv("data/gene_metadata.csv")

# genes x cells로 전치
cat("Transposing matrix...\n")
counts <- t(counts_raw)
rownames(counts) <- gene_meta$gene
colnames(counts) <- rownames(cell_meta)
rm(counts_raw)

# Seurat 객체 생성
cat("Creating Seurat object...\n")
seurat_obj <- CreateSeuratObject(counts = counts, meta.data = cell_meta)
rm(counts)

# Normalize
cat("Normalizing data...\n")
seurat_obj <- NormalizeData(seurat_obj)

# 2. Hepatitis + Mesenchymal만 추출 (predicted_celltype 사용)
cat("Subsetting Hepatitis Mesenchymal...\n")
seurat_sub <- subset(seurat_obj, subset = status == "Hepatitis" & predicted_celltype == "Mesenchymal")
rm(seurat_obj)

# nFeature_RNA NA 제거
seurat_sub <- subset(seurat_sub, cells = colnames(seurat_sub)[!is.na(seurat_sub$nFeature_RNA)])
cat(paste0("Hepatitis Mesenchymal after NA removal: ", ncol(seurat_sub), "\n"))

# 3. 비교군 설정: NT vs TC (group 컬럼)
Idents(seurat_sub) <- "group"
cat(paste0("NT: ", sum(Idents(seurat_sub) == "NT"),
           ", TC: ", sum(Idents(seurat_sub) == "TC"), "\n"))

# 4. MAST 실행 (NT vs TC)
cat("Running MAST DEG analysis (NT vs TC in Hepatitis Mesenchymal)...\n")
mast_results <- FindMarkers(
  object = seurat_sub,
  ident.1 = "TC",
  ident.2 = "NT",
  test.use = "MAST",
  latent.vars = "nFeature_RNA",
  logfc.threshold = 0.25
)

# 결과 정리: 컬럼명을 gene, log2FC, pvalue, padj로 변환
cat("Saving results...\n")
out_df <- data.frame(
  gene = rownames(mast_results),
  log2FC = mast_results$avg_log2FC,
  pvalue = mast_results$p_val,
  padj = mast_results$p_val_adj,
  stringsAsFactors = FALSE
)
write.csv(out_df, file.path(out_dir, "MAST_DEGs_Hepatitis_NTvsTC_Mesenchymal.csv"), row.names = FALSE)
cat(paste0("Total DEGs found: ", nrow(out_df), "\n"))
cat(paste0("Significant DEGs (padj < 0.05): ",
           sum(out_df$padj < 0.05), "\n"))

cat("Done! Results saved to Deg/results_output/\n")
