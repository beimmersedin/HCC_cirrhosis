# 1. 필수 라이브러리 로드
library(Signac)
library(Seurat)
library(Matrix)
library(data.table)
library(ggplot2)

# -----------------------------------------------------------------------------
# 2. 경로 설정 및 폴더 생성
# -----------------------------------------------------------------------------
# 데이터 위치 (입력)
data_dir <- "/data1/project2/yeonu/cirrhosis/multiome/"

# 결과 저장 위치 (출력: DA 폴더 내에 'results_output' 폴더 생성)
output_dir <- "/data1/project/yeonu/065_multi_atac/DA/results_output/"

# 폴더가 없으면 생성 (recursive = TRUE로 중간 경로까지 생성 가능)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("결과 저장 폴더가 생성되었습니다: ", output_dir)
}

# -----------------------------------------------------------------------------
# 3. 데이터 로드 및 객체 생성
# -----------------------------------------------------------------------------
# 메타데이터 로드 — rna_cell_metadata.csv 사용 (celltype_assign, group, status 포함)
metadata <- fread(paste0(data_dir, "rna_cell_metadata.csv"), data.table = FALSE)
rownames(metadata) <- metadata[, 1]

# Sparse Count 로드 및 Matrix 변환
sparse_data <- fread(paste0(data_dir, "atac_counts_sparse.csv"), data.table = FALSE)
counts <- sparseMatrix(
  i = as.numeric(as.factor(sparse_data$peak)),
  j = as.numeric(as.factor(sparse_data$cell_id)),
  x = sparse_data$value,
  dimnames = list(levels(as.factor(sparse_data$peak)), levels(as.factor(sparse_data$cell_id)))
)

# Signac 객체 생성 [cite: 14, 59, 83]
chrom_assay <- CreateChromatinAssay(counts = counts, sep = c("-", "-"), genome = 'hg38')
multiome_atac <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = metadata)

# -----------------------------------------------------------------------------
# 4. TC 샘플 서브셋팅 및 분석 [cite: 91, 94, 96]
# -----------------------------------------------------------------------------
# Tumor Core(TC)만 추출
tc_atac <- subset(multiome_atac, subset = group == "TC")

# 정규화 및 차원 축소 [cite: 10, 40]
tc_atac <- RunTFIDF(tc_atac)
tc_atac <- FindTopFeatures(tc_atac, min.cutoff = 'q0')
tc_atac <- RunSVD(tc_atac)


# # -----------------------------------------------------------------------------
# # 5. 세포 유형별 DA 분석 (예: Hepatocytes) [cite: 93, 111, 190, 192]
# # -----------------------------------------------------------------------------
# cell_type <- "Hepatocytes"
# tc_subset <- subset(tc_atac, subset = celltype_assign == cell_type)
# Idents(tc_subset) <- "status"


# -----------------------------------------------------------------------------
# 5. 전체 세포 유형 DA 분석 (for loop)
# -----------------------------------------------------------------------------
cell_types <- unique(tc_atac$celltype_assign)
message("분석 대상 세포 유형 (", length(cell_types), "개): ", paste(cell_types, collapse = ", "))

for (cell_type in cell_types) {

  message("\n>>> [", which(cell_types == cell_type), "/", length(cell_types), "] ", cell_type, " 분석 시작...")

  # 해당 세포 유형 서브셋
  tc_subset <- subset(tc_atac, subset = celltype_assign == cell_type)
  n_cirr <- sum(tc_subset$status == "Cirrhosis")
  n_hep  <- sum(tc_subset$status == "Hepatitis")
  message("    Cells: ", ncol(tc_subset), " (Cirrhosis=", n_cirr, ", Hepatitis=", n_hep, ")")

  Idents(tc_subset) <- "status"

  # DA 분석 수행 (세포 수 부족 시 에러 방지)
  da_results <- tryCatch({
    FindMarkers(
      object = tc_subset,
      ident.1 = "Cirrhosis",
      ident.2 = "Hepatitis",
      test.use = 'LR',
      latent.vars = 'nCount_ATAC',
      min.pct = 0.05
    )
  }, error = function(e) {
    message("    ⚠ 분석 실패 (", e$message, ") → SKIP")
    return(NULL)
  })

  if (is.null(da_results)) next

  # 결과 정리: 컬럼명 변경 (peak, log2FC, pct.1, pct.2, pvalue, padj)
  da_df <- data.frame(
    peak = rownames(da_results),
    log2FC = da_results$avg_log2FC,
    pct.1 = da_results$pct.1,
    pct.2 = da_results$pct.2,
    pvalue = da_results$p_val,
    padj = da_results$p_val_adj,
    stringsAsFactors = FALSE
  )

  out_file <- paste0(output_dir, "DA_results_", cell_type, ".csv")
  write.csv(da_df, out_file, row.names = FALSE)

  n_sig <- sum(da_df$padj < 0.05, na.rm = TRUE)
  message("    완료: ", nrow(da_df), " peaks 테스트, ", n_sig, " significant → ", out_file)
}

message("\n모든 세포 유형 DA 분석 완료. 결과: ", output_dir)