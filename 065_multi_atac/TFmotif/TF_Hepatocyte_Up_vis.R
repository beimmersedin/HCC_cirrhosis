# =============================================================================
#  TF Motif Analysis — 시각화 (TF_Hepatocyte_Up.R 결과 기반)
#
#  입력: Hepatocyte_Motif_Enrichment.csv, Hepatocyte_TF_activity.csv
#  출력: 01_Top12_Motif_Logos.png
#        02_Top3_TF_Activity_UMAP.png
#        03_Top3_TF_Violin_CirrhosisVsHepatitis.png
#
#  주의: TF_Hepatocyte_Up.R 실행 후 tc_subset 객체가 메모리에 있는 상태에서
#        실행하거나, 아래 데이터 로드 섹션을 활성화하여 독립 실행 가능
# =============================================================================

# ─── 라이브러리 ──────────────────────────────────────────────────────────────
library(Signac)
library(Seurat)
library(Matrix)
library(data.table)
library(JASPAR2022)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(BiocParallel)
library(motifmatchr)
library(ggplot2)

# ─── 경로 설정 ───────────────────────────────────────────────────────────────
data_dir   <- "/data1/project2/yeonu/cirrhosis/multiome/"
output_dir <- "/data1/project/yeonu/065_multi_atac/TFmotif/results_output/"

# =============================================================================
# 데이터 로드 (독립 실행 시 활성화)
# TF_Hepatocyte_Up.R 직후 실행하면 tc_subset, enriched_motifs, tf_scores가
# 이미 메모리에 있으므로 이 섹션은 건너뛸 수 있음
# =============================================================================
if (!exists("tc_subset")) {
  message("=== tc_subset 없음 → 데이터 재로드 ===")

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
  tc_subset     <- subset(multiome_atac, subset = group == "TC" & celltype_assign == "Hepatocytes")

  tc_subset <- RunTFIDF(tc_subset) %>% FindTopFeatures(min.cutoff = "q0") %>% RunSVD()

  pfm <- getMatrixSet(x = JASPAR2022, opts = list(species = 9606, all_versions = FALSE))
  tc_subset <- AddMotifs(object = tc_subset, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)

  register(MulticoreParam(workers = 4))
  tc_subset <- RunChromVAR(object = tc_subset, genome = BSgenome.Hsapiens.UCSC.hg38)

  message("데이터 재로드 완료: ", ncol(tc_subset), " cells")
}

# CSV 로드
enriched_motifs <- read.csv(paste0(output_dir, "Hepatocyte_Motif_Enrichment.csv"), row.names = 1)
tf_scores       <- read.csv(paste0(output_dir, "Hepatocyte_TF_activity.csv"), row.names = 1)

message("Enriched motifs: ", nrow(enriched_motifs))
message("TF scores: ", nrow(tf_scores), " motifs × ", ncol(tf_scores), " cells")


# =============================================================================
# 1. Top 12 Motif Logo
# =============================================================================
message("\n", strrep("=", 70))
message("시각화 1: Top 12 Motif Logos")
message(strrep("=", 70))

p1 <- MotifPlot(tc_subset, motifs = head(rownames(enriched_motifs), 12))
ggsave(paste0(output_dir, "01_Top12_Motif_Logos.png"), plot = p1,
       width = 10, height = 12, dpi = 300)
message("01_Top12_Motif_Logos.png 저장 완료")


# =============================================================================
# 2. Top 3 TF Activity UMAP
# =============================================================================
message("\n", strrep("=", 70))
message("시각화 2: Top 3 TF Activity UMAP")
message(strrep("=", 70))

top_tfs <- head(enriched_motifs[order(enriched_motifs$p.adjust), ], 3)
tc_subset <- RunUMAP(tc_subset, reduction = "lsi", dims = 2:30)
DefaultAssay(tc_subset) <- "chromvar"

plots <- list()
for (i in seq_len(nrow(top_tfs))) {
  plots[[i]] <- FeaturePlot(tc_subset, features = top_tfs$motif[i],
                            min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.3) +
    scale_color_gradientn(colors = c("navy", "white", "firebrick")) +
    ggtitle(paste0(top_tfs$motif.name[i], " Activity"))
}
p2 <- wrap_plots(plots, ncol = 3)
ggsave(paste0(output_dir, "02_Top3_TF_Activity_UMAP.png"), plot = p2,
       width = 15, height = 5, dpi = 300)
message("02_Top3_TF_Activity_UMAP.png 저장 완료")


# =============================================================================
# 3. Cirrhosis vs Hepatitis Violin Plot (Top 3 TFs)
# =============================================================================
message("\n", strrep("=", 70))
message("시각화 3: Top 3 TF Violin (Cirrhosis vs Hepatitis)")
message(strrep("=", 70))

plots <- list()
for (i in seq_len(nrow(top_tfs))) {
  tf_id   <- top_tfs$motif[i]
  tf_name <- top_tfs$motif.name[i]
  plot_df <- data.frame(
    score  = as.numeric(tf_scores[tf_id, ]),
    status = tc_subset$status
  )
  plots[[i]] <- ggplot(plot_df, aes(x = status, y = score, fill = status)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
    scale_fill_manual(values = c("Cirrhosis" = "#C44E52", "Hepatitis" = "#8172B3")) +
    labs(title = paste0(tf_name, " Activity"), x = "Status", y = "chromVAR Deviation Score") +
    theme_classic() + theme(legend.position = "none")
}
p3 <- wrap_plots(plots, ncol = 1)
ggsave(paste0(output_dir, "03_Top3_TF_Violin_CirrhosisVsHepatitis.png"), plot = p3,
       width = 10, height = 12, dpi = 300)
message("03_Top3_TF_Violin_CirrhosisVsHepatitis.png 저장 완료")


# =============================================================================
# 완료
# =============================================================================
message("\n", strrep("=", 70))
message("시각화 완료! 결과 위치: ", output_dir)
message(strrep("=", 70))
message("\n출력 파일:")
message("  [시각화] 01_Top12_Motif_Logos.png")
message("  [시각화] 02_Top3_TF_Activity_UMAP.png")
message("  [시각화] 03_Top3_TF_Violin_CirrhosisVsHepatitis.png")
