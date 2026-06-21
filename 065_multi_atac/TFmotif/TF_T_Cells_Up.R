# =============================================================================
#  TF Motif Analysis Pipeline — TC T_Cells (Cirrhosis vs Hepatitis)
#
#  Step 1: 후보 TF 선정 (Motif Enrichment)       → T_Cells_Motif_Enrichment.csv
#  Step 2: 결합 부위 확정 (Peak-TF Mapping)       → T_Cells_Peak_TF_Match_Map.csv
#  Step 3: 세포별 활성 검증 (chromVAR)             → T_Cells_TF_activity.csv
# =============================================================================

# ─── 라이브러리 ──────────────────────────────────────────────────────────────
library(Signac)
library(Seurat)
library(Matrix)
library(data.table)
library(JASPAR2022)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)
library(motifmatchr)
library(reshape2)

# ─── 경로 설정 ───────────────────────────────────────────────────────────────
data_dir   <- "/data1/project2/yeonu/cirrhosis/multiome/"
da_dir     <- "/data1/project/yeonu/065_multi_atac/DA/results_output/"
output_dir <- "/data1/project/yeonu/065_multi_atac/TFmotif/results_output/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ─── 데이터 로드 & 객체 생성 ─────────────────────────────────────────────────
message("=== 데이터 로드 시작 ===")
metadata    <- fread(paste0(data_dir, "rna_cell_metadata.csv"), data.table = FALSE)
rownames(metadata) <- metadata[, 1]
sparse_data <- fread(paste0(data_dir, "atac_counts_sparse.csv"), data.table = FALSE)

counts <- sparseMatrix(
  i = as.numeric(as.factor(sparse_data$peak)),
  j = as.numeric(as.factor(sparse_data$cell_id)),
  x = sparse_data$value,
  dimnames = list(levels(as.factor(sparse_data$peak)),
                  levels(as.factor(sparse_data$cell_id)))
)

chrom_assay   <- CreateChromatinAssay(counts = counts, sep = c("-", "-"), genome = "hg38")
multiome_atac <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = metadata)
tc_subset     <- subset(multiome_atac, subset = group == "TC" & celltype_assign == "T_Cells")

tc_subset <- RunTFIDF(tc_subset) %>% FindTopFeatures(min.cutoff = "q0") %>% RunSVD()
message("객체 생성 완료: ", ncol(tc_subset), " cells")

# ─── JASPAR Motif DB 로드 & Peak에 Motif 매핑 ────────────────────────────────
pfm <- getMatrixSet(x = JASPAR2022, opts = list(species = 9606, all_versions = FALSE))
tc_subset <- AddMotifs(object = tc_subset, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)

# DA 결과 로드
sig_dars <- read.csv(paste0(da_dir, "Significant_DARs_T_Cells.csv"), row.names = 1)
up_peaks <- rownames(sig_dars[sig_dars$avg_log2FC > 0, ])
message("Up-regulated DARs (Cirrhosis): ", length(up_peaks), " peaks")

# ─── Significant DARs에 대한 PWM Score 계산 ──────────────────────────────────
message("PWM Score 계산 중 (Significant DARs only)...")
dar_ranges <- StringToGRanges(rownames(sig_dars), sep = c("-", "-"))
pwm_results <- matchMotifs(pfm,
                           subject = dar_ranges,
                           genome = BSgenome.Hsapiens.UCSC.hg38,
                           out = "scores")
pwm_score_matrix <- motifScores(pwm_results)
pwm_match_matrix <- motifMatches(pwm_results)
rownames(pwm_score_matrix) <- rownames(sig_dars)
rownames(pwm_match_matrix) <- rownames(sig_dars)
message("PWM Score 완료: ", nrow(pwm_score_matrix), " peaks × ", ncol(pwm_score_matrix), " motifs")

# Motif별 PWM score 임계값 계산 (p < 5e-05 통과 최소 score)
pwm_threshold_vec <- sapply(seq_len(ncol(pwm_score_matrix)), function(j) {
  matched <- pwm_match_matrix[, j]
  if (any(matched)) min(pwm_score_matrix[matched, j]) else NA
})
names(pwm_threshold_vec) <- colnames(pwm_score_matrix)
message("PWM Threshold 계산 완료: ", sum(!is.na(pwm_threshold_vec)), " motifs with matches")

# Motif별 이론적 최대 score 계산 (consensus 서열 완벽 매칭 시)
pwm_max_vec <- sapply(pfm, function(m) {
  pwm <- toPWM(m, type = "log2probratio")
  sum(apply(as.matrix(Matrix(pwm)), 2, max))
})
names(pwm_max_vec) <- sapply(pfm, function(x) x@ID)
message("PWM Max Score 계산 완료: ", length(pwm_max_vec), " motifs")


# =============================================================================
# Step 1. 후보 TF 선정 (Candidate Selection)
# =============================================================================
message("\n", strrep("=", 70))
message("Step 1: 후보 TF 선정 (Motif Enrichment)")
message(strrep("=", 70))

enriched_motifs <- FindMotifs(object = tc_subset, features = up_peaks)
write.csv(enriched_motifs, paste0(output_dir, "T_Cells_Motif_Enrichment.csv"))
message("저장 완료: T_Cells_Motif_Enrichment.csv (", nrow(enriched_motifs), " motifs)")


# =============================================================================
# Step 2. 결합 부위 확정 (Peak-TF Mapping)
# =============================================================================
message("\n", strrep("=", 70))
message("Step 2: 결합 부위 확정 (Peak-TF Mapping)")
message(strrep("=", 70))

# Motif match matrix (peak x motif, 0/1)
motif_matrix <- GetMotifData(object = tc_subset, assay = "peaks", slot = "data")

# Significant DARs만 필터링
common_peaks <- intersect(rownames(sig_dars), rownames(motif_matrix))
dar_motifs   <- motif_matrix[common_peaks, ]

# Long format 변환
peak_tf_list <- as.data.frame(as.matrix(dar_motifs))
peak_tf_list$peak <- rownames(peak_tf_list)
peak_tf_mapping <- melt(peak_tf_list, id.vars = "peak",
                        variable.name = "motif_id", value.name = "is_bound")
peak_tf_mapping <- peak_tf_mapping[peak_tf_mapping$is_bound == 1, ]

# Motif ID -> TF 이름 변환
motif_id_vec   <- sapply(pfm, function(x) x@ID)
motif_name_vec <- sapply(pfm, function(x) x@name)
names(motif_name_vec) <- motif_id_vec
peak_tf_mapping$tf_name    <- motif_name_vec[as.character(peak_tf_mapping$motif_id)]
peak_tf_mapping$avg_log2FC <- sig_dars[peak_tf_mapping$peak, "avg_log2FC"]
peak_tf_mapping$p_val_adj  <- sig_dars[peak_tf_mapping$peak, "p_val_adj"]

# PWM Score 추가 (각 peak-motif 조합의 최대 매칭 점수)
pwm_idx <- cbind(match(peak_tf_mapping$peak, rownames(pwm_score_matrix)),
                 match(as.character(peak_tf_mapping$motif_id), colnames(pwm_score_matrix)))
peak_tf_mapping$pwm_score <- pwm_score_matrix[pwm_idx]

# Motif별 PWM 임계값 추가 (p < 5e-05 통과 기준 score)
peak_tf_mapping$pwm_threshold <- pwm_threshold_vec[as.character(peak_tf_mapping$motif_id)]

# Motif별 이론적 최대 score 추가
peak_tf_mapping$pwm_max <- pwm_max_vec[as.character(peak_tf_mapping$motif_id)]

write.csv(peak_tf_mapping, paste0(output_dir, "T_Cells_Peak_TF_Match_Map.csv"), row.names = FALSE)
message("저장 완료: T_Cells_Peak_TF_Match_Map.csv (", nrow(peak_tf_mapping), " pairs)")


# =============================================================================
# Step 3. 세포별 활성 검증 (chromVAR)
# =============================================================================
message("\n", strrep("=", 70))
message("Step 3: 세포별 활성 검증 (chromVAR)")
message(strrep("=", 70))

register(MulticoreParam(workers = 4))
tc_subset <- RunChromVAR(object = tc_subset, genome = BSgenome.Hsapiens.UCSC.hg38)

tf_scores <- GetAssayData(tc_subset, assay = "chromvar", layer = "data")
write.csv(as.matrix(tf_scores), paste0(output_dir, "T_Cells_TF_activity.csv"))
message("저장 완료: T_Cells_TF_activity.csv")


# =============================================================================
# 완료
# =============================================================================
message("\n", strrep("=", 70))
message("TF 분석 완료! 결과 위치: ", output_dir)
message(strrep("=", 70))
message("\n출력 파일:")
message("  [Step 1] T_Cells_Motif_Enrichment.csv")
message("  [Step 2] T_Cells_Peak_TF_Match_Map.csv")
message("  [Step 3] T_Cells_TF_activity.csv")
