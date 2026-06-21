import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

# ==========================================
# 1. 경로 및 파라미터 설정
# ==========================================
BASE_DIR = "/data1/project/yeonu/065_multi_rna"
DATA_PATH = os.path.join(BASE_DIR, "data/combined_final_annotated.h5ad")
GENE_80_PATH = os.path.join(BASE_DIR, "data/gene_metadata_80.csv")
OUT_DIR = os.path.join(BASE_DIR, "WGCNA_80/results_output")

MIN_CELLS_PER_BIN = 50     # 한 빈에 묶을 세포 수 (Pseudobulk 1개 생성 단위)

os.makedirs(OUT_DIR, exist_ok=True)

def make_pseudobulk(adata_full, celltype, out_dir, bin_size):
    """단일 cell type에 대해 pseudobulk를 생성하고 저장한다."""
    print(f"\n{'='*50}")
    print(f"--- Processing: {celltype} ---")
    print(f"{'='*50}")

    # 1) 특정 Cell Type 및 NT/TC 그룹만 추출
    keep = (adata_full.obs["predicted_celltype"] == celltype) & (adata_full.obs["group"].isin(["NT", "TC"]))
    adata = adata_full[keep].copy()
    print(f"Genes: {adata.n_vars}, Cells: {adata.n_obs}")

    if adata.n_obs < bin_size:
        print(f"--- Skipped: {celltype} (only {adata.n_obs} cells, need >= {bin_size}) ---")
        return

    # Raw count 확보
    X_raw = adata.X
    if not sparse.issparse(X_raw):
        X_raw = sparse.csr_matrix(X_raw)

    # 2) Pseudobulk 생성 (Random Binning)
    print(f"--- Generating Pseudobulk (Bin size: {bin_size}) ---")
    rng = np.random.default_rng(42)

    pb_mats = []
    pb_meta = []

    group_cols = ["sample_id", "group", "status", "dataset"]
    obs = adata.obs.copy()
    obs["_temp_idx"] = np.arange(adata.n_obs)

    for keys, sub in obs.groupby(group_cols, observed=True):
        sample_id, group, status, dataset = keys
        indices = sub["_temp_idx"].values

        if len(indices) < bin_size:
            continue

        rng.shuffle(indices)
        n_bins = len(indices) // bin_size

        for b in range(n_bins):
            chunk = indices[b * bin_size : (b + 1) * bin_size]

            summed_counts = X_raw[chunk].sum(axis=0)
            pb_mats.append(np.asarray(summed_counts).ravel())

            pb_meta.append({
                "pb_id": f"{sample_id}_bin{b+1}",
                "sample_id": sample_id,
                "group": group,
                "status": status,
                "dataset": dataset,
                "cell_count": len(chunk)
            })

    if len(pb_mats) == 0:
        print(f"--- Skipped: {celltype} (no bins generated) ---")
        return

    # 결과물 AnnData 생성
    pb_obs = pd.DataFrame(pb_meta).set_index("pb_id")
    pb_adata = ad.AnnData(X=np.vstack(pb_mats), obs=pb_obs, var=adata.var.copy())

    # 3) 정규화 (CPM-like -> log1p)
    print("--- Normalizing Pseudobulk Data ---")
    sc.pp.normalize_total(pb_adata, target_sum=1e6)
    sc.pp.log1p(pb_adata)

    # 0이 너무 많은 유전자 최종 2차 필터 (전체 pseudobulk의 5% 이상 발현)
    sc.pp.filter_genes(pb_adata, min_cells=int(pb_adata.n_obs * 0.05))

    # 4) 저장
    output_filename = f"{celltype}_80_pb_n{bin_size}.h5ad"
    save_path = os.path.join(out_dir, output_filename)
    pb_adata.write_h5ad(save_path)

    print(f"--- Done: {celltype} ---")
    print(f"Final Pseudobulk shape: {pb_adata.n_obs} samples x {pb_adata.n_vars} genes")
    print(f"Saved to: {save_path}")

def main():
    # ==========================================
    # 2. 데이터 로드
    # ==========================================
    print("--- Loading full dataset ---")
    adata_full = sc.read_h5ad(DATA_PATH)

    # NT/TC 그룹만 유지
    adata_full = adata_full[adata_full.obs["group"].isin(["NT", "TC"])].copy()

    # zero rate < 80% 유전자만 필터링
    print("--- Filtering to zero rate < 80% genes ---")
    gene_80 = pd.read_csv(GENE_80_PATH)
    keep_genes = gene_80["gene"].tolist()
    common_genes = [g for g in keep_genes if g in adata_full.var_names]
    adata_full = adata_full[:, common_genes].copy()
    print(f"Genes after 80% filter: {adata_full.n_vars}")

    # 모든 cell type 추출
    celltypes = adata_full.obs["predicted_celltype"].value_counts().index.tolist()
    print(f"Cell types to process: {celltypes}")

    # ==========================================
    # 3. Cell type별 Pseudobulk 생성
    # ==========================================
    for ct in celltypes:
        make_pseudobulk(adata_full, ct, OUT_DIR, MIN_CELLS_PER_BIN)

    print(f"\n{'='*50}")
    print("--- All cell types processed! ---")
    print(f"{'='*50}")

if __name__ == "__main__":
    main()
