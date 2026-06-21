# =============================================================================
#  TF-Target Gene 매칭 — TC Macrophages
#  기존 Peak-Gene Links (LinkPeaks 결과)와 TFmotif Peak-TF Mapping을 조인
#  ※ LinkPeaks는 LinkPeaks_Macrophages.R에서 별도 수행 (이미 완료된 경우 재실행 불필요)
# =============================================================================

library(data.table)

# ─── 경로 설정 ───────────────────────────────────────────────────────────────
tf_dir     <- "/data1/project/yeonu/065_multi_atac/TFmotif/results_output/"
linkage_dir <- "/data1/project/yeonu/065_multi_atac/Linkage/results_output/"
output_dir <- "/data1/project/yeonu/065_multi_atac/Linkage/results_output/"

# ─── 데이터 로드 ─────────────────────────────────────────────────────────────
message("=== TF-Target Gene 매칭 시작 (Macrophages) ===")

# Peak-Gene Links (LinkPeaks 결과, 이미 저장된 CSV)
all_links <- fread(paste0(linkage_dir, "Macrophages_Peak_Gene_Links.csv"), data.table = FALSE)
message("Peak-Gene Links 로드 완료: ", nrow(all_links), " links")

# Peak-TF Mapping (TFmotif 결과)
peak_tf_mapping <- fread(paste0(tf_dir, "Macrophages_Peak_TF_Match_Map.csv"), data.table = FALSE)
message("Peak-TF Mapping 로드 완료: ", nrow(peak_tf_mapping), " mappings")

# =============================================================================
# TF-Target Gene 매칭 (TFmotif 결과와 연결)
# =============================================================================
message("\n", strrep("=", 70))
message("TF-Target Gene 매칭")
message(strrep("=", 70))

# Peak-TF Mapping과 Peak-Gene Link를 조인하여 TF → Peak → Gene 관계 생성
tf_gene_links <- merge(
  peak_tf_mapping[, c("peak", "motif_id", "tf_name", "avg_log2FC", "p_val_adj")],
  all_links[, c("peak", "gene", "score", "zscore", "pvalue")],
  by = "peak"
)
tf_gene_links <- tf_gene_links[order(tf_gene_links$pvalue), ]

write.csv(tf_gene_links, paste0(output_dir, "Macrophages_TF_Peak_Gene_Links.csv"), row.names = FALSE)
message("저장 완료: Macrophages_TF_Peak_Gene_Links.csv (", nrow(tf_gene_links), " TF-Peak-Gene triplets)")

# =============================================================================
# 완료
# =============================================================================
message("\n", strrep("=", 70))
message("TF-Target Gene 매칭 완료! 결과 위치: ", output_dir)
message(strrep("=", 70))
message("\n출력 파일:")
message("  Macrophages_TF_Peak_Gene_Links.csv   — TF → Peak → Gene 통합 매핑")
