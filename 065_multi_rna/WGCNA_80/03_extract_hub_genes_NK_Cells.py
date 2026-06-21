import os
import pandas as pd
import PyWGCNA

# ==========================================
# 1. 경로 및 설정
# ==========================================
BASE_DIR = "/data1/project/yeonu/065_multi_rna"
WGCNA_OBJ_PATH = os.path.join(BASE_DIR, "WGCNA_80/results_output/NK_Cells/NK_Cells_Network.p")
OUT_DIR = os.path.join(BASE_DIR, "WGCNA_80/results_output/NK_Cells/hub_genes")

os.makedirs(OUT_DIR, exist_ok=True)

def main():
    # ==========================================
    # 2. PyWGCNA 객체 로드
    # ==========================================
    print(f"--- Loading PyWGCNA object: {WGCNA_OBJ_PATH} ---")
    pywgcna = PyWGCNA.readWGCNA(WGCNA_OBJ_PATH)

    # ==========================================
    # 3. 유전자 정보 추출 (Module & kME)
    # ==========================================
    gene_info = pywgcna.datExpr.var.copy()
    gene_info = pd.concat([gene_info, pywgcna.signedKME], axis=1)

    # 전체 유전자 리스트 저장
    gene_info.to_csv(os.path.join(OUT_DIR, "all_genes_wgcna_info.csv"))

    # ==========================================
    # 4. 모듈별 Top Hub 유전자 추출
    # ==========================================
    modules = gene_info['moduleColors'].unique()
    print(f"Detected modules: {len(modules)}")

    hub_summary = []

    for module in modules:
        if module == 'grey': continue

        module_genes = gene_info[gene_info['moduleColors'] == module].copy()

        kme_col = f"kME{module}"
        if kme_col in module_genes.columns:
            module_genes = module_genes.sort_values(by=kme_col, ascending=False)
            module_genes.to_csv(os.path.join(OUT_DIR, f"genes_{module}.csv"))

            hub_summary.append({
                "Module": module,
                "Size": module_genes.shape[0],
                "Top_Hub": module_genes.index[0],
                "Top_kME": module_genes[kme_col].iloc[0]
            })

    pd.DataFrame(hub_summary).to_csv(os.path.join(OUT_DIR, "module_hub_summary.csv"), index=False)

    print(f"--- Extraction Done! ---")
    print(f"Hub gene lists saved in: {OUT_DIR}")

if __name__ == "__main__":
    main()
