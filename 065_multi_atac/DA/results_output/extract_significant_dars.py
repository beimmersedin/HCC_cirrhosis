"""
DA_results_*.csv에서 padj < 0.05 기준으로 Significant DARs를 추출하여 저장
출력 형식: Significant_DARs_Hepatocytes.csv와 동일 (peak=rowname, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
"""
import pandas as pd
import os

da_dir = "/data1/project/yeonu/065_multi_atac/DA/results_output/"
threshold = 0.05

cell_types = [
    "Hepatocytes", "T_Cells", "Mesenchymal", "Macrophages",
    "NK_Cells", "DCs", "Plasma_Cells", "B_Cells"
]

for ct in cell_types:
    infile = os.path.join(da_dir, "DA_results_{}.csv".format(ct))
    if not os.path.exists(infile):
        continue

    da = pd.read_csv(infile)
    sig = da[da["padj"] < threshold].copy()

    if len(sig) == 0:
        print("{:<20} → 0 DARs (skip)".format(ct))
        continue

    # 기존 Significant_DARs_Hepatocytes.csv 형식에 맞춤
    out = pd.DataFrame({
        "p_val": sig["pvalue"].values,
        "avg_log2FC": sig["log2FC"].values,
        "pct.1": sig["pct.1"].values,
        "pct.2": sig["pct.2"].values,
        "p_val_adj": sig["padj"].values,
    }, index=sig["peak"].values)

    outfile = os.path.join(da_dir, "Significant_DARs_{}.csv".format(ct))
    out.to_csv(outfile)
    n_up = (out["avg_log2FC"] > 0).sum()
    n_down = (out["avg_log2FC"] < 0).sum()
    print("{:<20} → {} DARs (Up: {}, Down: {}) → {}".format(ct, len(out), n_up, n_down, os.path.basename(outfile)))
