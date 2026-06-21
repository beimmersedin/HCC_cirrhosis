# =============================================================================
#  Peak-to-Gene LinkPeaks — TC NK_Cells
#  ATAC Peak과 RNA 발현 간 상관관계를 계산 (시간이 오래 걸리는 단계)
#  결과: NK_Cells_Peak_Gene_Links.csv
# =============================================================================

# ─── 라이브러리 ──────────────────────────────────────────────────────────────
library(Signac)
library(Seurat)
library(Matrix)
library(data.table)
library(ggplot2)
library(magrittr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

# ─── 경로 설정 ───────────────────────────────────────────────────────────────
data_dir   <- "/data1/project2/yeonu/cirrhosis/multiome/"
output_dir <- "/data1/project/yeonu/065_multi_atac/Linkage/results_output/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ─── 데이터 로드 & ATAC 객체 생성 ────────────────────────────────────────────
message("=== 데이터 로드 시작 ===")

# 메타데이터
metadata <- fread(paste0(data_dir, "rna_cell_metadata.csv"), data.table = FALSE)
rownames(metadata) <- metadata[, 1]

# ATAC sparse count
sparse_data <- fread(paste0(data_dir, "atac_counts_sparse.csv"), data.table = FALSE)
atac_counts <- sparseMatrix(
  i = as.numeric(as.factor(sparse_data$peak)),
  j = as.numeric(as.factor(sparse_data$cell_id)),
  x = sparse_data$value,
  dimnames = list(levels(as.factor(sparse_data$peak)),
                  levels(as.factor(sparse_data$cell_id)))
)

# RNA sparse count
rna_sparse <- fread(paste0(data_dir, "rna_counts_sparse.csv"), data.table = FALSE)
rna_counts <- sparseMatrix(
  i = as.numeric(as.factor(rna_sparse$gene)),
  j = as.numeric(as.factor(rna_sparse$cell_id)),
  x = rna_sparse$value,
  dimnames = list(levels(as.factor(rna_sparse$gene)),
                  levels(as.factor(rna_sparse$cell_id)))
)

# ─── Gene Annotation 설정 (EnsDb.Hsapiens.v86) ──────────────────────────────
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"

# ─── Seurat 객체 생성 ────────────────────────────────────────────────────────
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts, sep = c("-", "-"), genome = "hg38"
)
Annotation(chrom_assay) <- annotations

multiome_atac <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = metadata)

# TC NK_Cells 서브셋
tc_subset <- subset(multiome_atac, subset = group == "TC" & celltype_assign == "NK_Cells")
message("ATAC 객체 생성 완료: ", ncol(tc_subset), " cells")

# ATAC 정규화
tc_subset <- RunTFIDF(tc_subset) %>% FindTopFeatures(min.cutoff = "q0") %>% RunSVD()

# ─── RNA Assay 추가 ──────────────────────────────────────────────────────────
common_cells <- intersect(colnames(tc_subset), colnames(rna_counts))
tc_subset <- subset(tc_subset, cells = common_cells)
rna_sub   <- rna_counts[, common_cells]

tc_subset[["RNA"]] <- CreateAssayObject(counts = rna_sub)

# RNA 정규화
DefaultAssay(tc_subset) <- "RNA"
tc_subset <- SCTransform(tc_subset, verbose = FALSE)
message("RNA 통합 완료: ", ncol(tc_subset), " cells, ", nrow(rna_sub), " genes")

# =============================================================================
# Peak-to-Gene Linkage 수행
# =============================================================================
message("\n", strrep("=", 70))
message("Peak-to-Gene Linkage 수행 (distance = 500kb)")
message(strrep("=", 70))

DefaultAssay(tc_subset) <- "peaks"

# DNA sequence 정보 계산 (LinkPeaks 필수 선행 단계)
tc_subset <- RegionStats(tc_subset, genome = BSgenome.Hsapiens.UCSC.hg38)

tc_subset <- LinkPeaks(
  object = tc_subset,
  peak.assay = "peaks",
  expression.assay = "SCT",
  distance = 5e+05
)

# Link 정보 추출 및 저장
all_links <- as.data.frame(Links(tc_subset))
write.csv(all_links, paste0(output_dir, "NK_Cells_Peak_Gene_Links.csv"))
message("저장 완료: NK_Cells_Peak_Gene_Links.csv (", nrow(all_links), " links)")

message("\n", strrep("=", 70))
message("LinkPeaks 완료! 결과 위치: ", output_dir)
message(strrep("=", 70))
