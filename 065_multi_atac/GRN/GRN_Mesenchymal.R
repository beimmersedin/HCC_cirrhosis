# =============================================================================
#  GRN (Gene Regulatory Network) — TC Mesenchymal
#  chromVAR TF Activity와 RNA 발현의 cell-level 상관관계로 TF-Target 교차 검증
# =============================================================================

# 1. 필수 라이브러리 로드
library(data.table)
library(Matrix)
library(dplyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# 2. 경로 설정
# -----------------------------------------------------------------------------
data_dir    <- "/data1/project2/yeonu/cirrhosis/multiome/"
tf_dir      <- "/data1/project/yeonu/065_multi_atac/TFmotif/results_output/"
linkage_dir <- "/data1/project/yeonu/065_multi_atac/Linkage/results_output/"
output_dir  <- "/data1/project/yeonu/065_multi_atac/GRN/results_output/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 3. 데이터 로드
# -----------------------------------------------------------------------------
message("=== 데이터 로드 시작 ===")

# (A) TF-Peak-Gene Triplets (Linkage 결과)
triplets <- fread(paste0(linkage_dir, "Mesenchymal_TF_Peak_Gene_Links.csv"), data.table = FALSE)
message("Triplets 로드 완료: ", nrow(triplets), " rows")

# (B) 세포별 TF Activity Score (chromVAR 결과, motif × cell)
tf_scores <- fread(paste0(tf_dir, "Mesenchymal_TF_activity.csv"), data.table = FALSE)
rownames(tf_scores) <- tf_scores[, 1]
tf_scores <- tf_scores[, -1]
message("TF scores 로드 완료: ", nrow(tf_scores), " motifs × ", ncol(tf_scores), " cells")

# (C) RNA 발현 데이터 로드 (tc_subset 객체 대신 직접 구축)
metadata <- fread(paste0(data_dir, "rna_cell_metadata.csv"), data.table = FALSE)
rownames(metadata) <- metadata[, 1]

message("RNA sparse count 로드 중 (시간 소요)...")
rna_sparse <- fread(paste0(data_dir, "rna_counts_sparse.csv"), data.table = FALSE)
rna_counts <- sparseMatrix(
  i = as.numeric(as.factor(rna_sparse$gene)),
  j = as.numeric(as.factor(rna_sparse$cell_id)),
  x = rna_sparse$value,
  dimnames = list(levels(as.factor(rna_sparse$gene)),
                  levels(as.factor(rna_sparse$cell_id)))
)
rm(rna_sparse); gc()
message("RNA matrix 구축 완료: ", nrow(rna_counts), " genes × ", ncol(rna_counts), " cells")

# -----------------------------------------------------------------------------
# 4. 공통 세포 정렬 (tf_scores와 rna_counts 모두에 존재하는 TC Mesenchymal)
# -----------------------------------------------------------------------------
tc_cells <- rownames(metadata[metadata$group == "TC" & metadata$celltype_assign == "Mesenchymal", ])
common_cells <- Reduce(intersect, list(colnames(tf_scores), colnames(rna_counts), tc_cells))
message("공통 세포 수: ", length(common_cells))

tf_scores <- tf_scores[, common_cells]

# -----------------------------------------------------------------------------
# 5. RNA 정규화 (타겟 유전자만 추출 후 log-normalization)
# -----------------------------------------------------------------------------
# 전체 library size 계산 (정규화 정확도 위해)
lib_size <- colSums(rna_counts[, common_cells])

# 타겟 유전자만 추출하여 메모리 절약
target_genes <- unique(triplets$gene)
available_genes <- intersect(target_genes, rownames(rna_counts))
rna_sub <- rna_counts[available_genes, common_cells]
message("타겟 유전자: ", length(target_genes), "개 중 ", length(available_genes), "개 매칭")

# Log-normalization (TPM-like: counts / library_size * 10000 → log1p)
rna_matrix <- log1p(t(t(rna_sub) / lib_size * 10000))
rna_matrix <- as.matrix(rna_matrix)
rm(rna_counts, rna_sub); gc()

# -----------------------------------------------------------------------------
# 6. TF Activity – Gene Expression 상관관계 계산 (Spearman)
# -----------------------------------------------------------------------------
unique_pairs <- triplets %>% select(tf_name, motif_id, gene) %>% distinct()
message("\n상관관계 계산 시작 (총 ", nrow(unique_pairs), " 쌍)...")

calc_cor <- function(tf_id, gene_name) {
  if (tf_id %in% rownames(tf_scores) && gene_name %in% rownames(rna_matrix)) {
    return(cor(as.numeric(tf_scores[tf_id, ]),
               as.numeric(rna_matrix[gene_name, ]),
               method = "spearman"))
  } else {
    return(NA)
  }
}

unique_pairs$grn_weight <- mapply(calc_cor, unique_pairs$motif_id, unique_pairs$gene)

n_computed <- sum(!is.na(unique_pairs$grn_weight))
message("계산 완료: ", n_computed, " / ", nrow(unique_pairs), " 쌍 (NA 제외)")

# -----------------------------------------------------------------------------
# 7. 전체 GRN 결과 저장 (필터 전) — causal evidence 분석용
# -----------------------------------------------------------------------------
all_grn <- unique_pairs %>%
  filter(!is.na(grn_weight)) %>%
  arrange(desc(grn_weight))

write.csv(all_grn, paste0(output_dir, "Mesenchymal_All_GRN_Weights.csv"), row.names = FALSE)
message("전체 GRN 저장 완료: ", nrow(all_grn), " edges (필터 전)")

# -----------------------------------------------------------------------------
# 8. 필터링된 GRN 리스트 생성 및 Master Regulator 식별
# -----------------------------------------------------------------------------
final_grn <- all_grn %>%
  filter(grn_weight > 0.1)

write.csv(final_grn, paste0(output_dir, "Mesenchymal_Final_Weighted_GRN.csv"), row.names = FALSE)
message("Weighted GRN 저장 완료: ", nrow(final_grn), " edges (weight > 0.1)")

# Master Regulator 식별 (가장 많은 타겟 유전자를 거느린 TF)
master_regulators <- final_grn %>%
  count(tf_name) %>%
  arrange(desc(n))

write.csv(master_regulators, paste0(output_dir, "Mesenchymal_Master_Regulators_Summary.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 완료
# -----------------------------------------------------------------------------
message("\n", strrep("=", 70))
message("GRN 구축 완료!")
message(strrep("=", 70))
message("Top 5 Master Regulators:")
print(head(master_regulators, 5))
message("\n결과 위치: ", output_dir)
message("  Mesenchymal_All_GRN_Weights.csv              — 전체 TF-Gene 상관 (causal evidence 입력용)")
message("  Mesenchymal_Final_Weighted_GRN.csv           — TF-Gene edge list (weight > 0.1)")
message("  Mesenchymal_Master_Regulators_Summary.csv    — TF별 타겟 유전자 수")
