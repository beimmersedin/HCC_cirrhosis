import os
import scanpy as sc
import PyWGCNA

# ==========================================
# 1. 경로 및 파일 설정
# ==========================================
BASE_DIR = "/data1/project/yeonu/065_multi_rna"
INPUT_PATH = os.path.join(BASE_DIR, "WGCNA_80/results_output/Endothelial_Cells_80_pb_n50.h5ad")
OUT_DIR = os.path.join(BASE_DIR, "WGCNA_80/results_output/Endothelial_Cells/")

os.makedirs(OUT_DIR, exist_ok=True)

def main():
    # ==========================================
    # 2. 데이터 로드
    # ==========================================
    print(f"--- Loading Pseudobulk: {INPUT_PATH} ---")
    adata = sc.read_h5ad(INPUT_PATH)

    # ==========================================
    # 3. PyWGCNA 객체 초기화
    # ==========================================
    pywgcna = PyWGCNA.WGCNA(
        name="Endothelial_Cells_Network",
        species="homo sapiens",
        anndata=adata,
        TPMcutoff=0,
        networkType="signed hybrid",
        TOMType="signed",
        minModuleSize=30,
        MEDissThres=0.2,
        outputPath=OUT_DIR,
        save=True,
        figureType="png"
    )

    # ==========================================
    # 4. 분석 파이프라인 실행 (모듈 탐색만)
    # ==========================================

    # [Step 1] Preprocess: outlier 샘플/유전자 제거
    pywgcna.preprocess(show=False)

    # [Step 2] 모듈 탐색 (soft threshold + adjacency + TOM + clustering + merge)
    pywgcna.findModules()

    # [Step 3] Module membership (kME) 계산
    pywgcna.CalculateSignedKME()

    # ==========================================
    # 5. 결과 저장
    # ==========================================
    pywgcna.saveWGCNA()

    print(f"--- Analysis Done! ---")
    print(f"Results saved in: {OUT_DIR}")

if __name__ == "__main__":
    main()
