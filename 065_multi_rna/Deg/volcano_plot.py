import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

out_dir = "Deg/results_output"
# 출력 디렉토리가 없으면 생성
os.makedirs(out_dir, exist_ok=True)

csv_files = [f for f in os.listdir(out_dir) if f.startswith("MAST_DEGs_") and f.endswith(".csv")]

# 설정값
p_cutoff = 0.05
fc_cutoff = 0.5

# --- 1st pass: 전체 CSV에서 글로벌 축 범위 계산 ---
global_x_min, global_x_max = 0, 0
global_y_max = 0

for fname in csv_files:
    df = pd.read_csv(os.path.join(out_dir, fname))
    min_nonzero_p = df.loc[df["padj"] > 0, "padj"].min() if any(df["padj"] > 0) else 1e-300
    df["padj_fixed"] = df["padj"].replace(0, min_nonzero_p * 0.1)
    df["neg_log10_padj"] = -np.log10(df["padj_fixed"])

    global_x_min = min(global_x_min, df["log2FC"].min())
    global_x_max = max(global_x_max, df["log2FC"].max())
    global_y_max = max(global_y_max, df["neg_log10_padj"].max())

# 대칭 x축 + 여유 패딩
x_abs = max(abs(global_x_min), abs(global_x_max))
x_pad = x_abs * 0.1
xlim = (-(x_abs + x_pad), x_abs + x_pad)
ylim = (0, global_y_max * 1.05)

print(f"Fixed axis: xlim={xlim[0]:.2f}~{xlim[1]:.2f}, ylim={ylim[0]:.2f}~{ylim[1]:.2f}")

# --- 2nd pass: 고정 축으로 volcano plot 생성 ---
for fname in csv_files:
    print(f"Processing: {fname}")
    df = pd.read_csv(os.path.join(out_dir, fname))

    # -log10(padj) 계산 및 0 처리
    min_nonzero_p = df.loc[df["padj"] > 0, "padj"].min() if any(df["padj"] > 0) else 1e-300
    df["padj_fixed"] = df["padj"].replace(0, min_nonzero_p * 0.1)
    df["neg_log10_padj"] = -np.log10(df["padj_fixed"])

    # Significance 분류
    df["sig"] = "NS"
    df.loc[(df["log2FC"] >= fc_cutoff) & (df["padj"] < p_cutoff), "sig"] = "Up"
    df.loc[(df["log2FC"] <= -fc_cutoff) & (df["padj"] < p_cutoff), "sig"] = "Down"

    colors = {"Down": "#2166AC", "NS": "lightgrey", "Up": "#B2182B"}

    fig, ax = plt.subplots(figsize=(10, 8))

    # 산점도 그리기
    for sig, color in colors.items():
        sub = df[df["sig"] == sig]
        ax.scatter(sub["log2FC"], sub["neg_log10_padj"],
                   c=color, s=20, alpha=0.6, label=sig, edgecolors='none', zorder=3)

    # 가이드라인
    ax.axhline(-np.log10(p_cutoff), linestyle="--", color="black", linewidth=0.8, alpha=0.5, zorder=1)
    ax.axvline(fc_cutoff, linestyle="--", color="black", linewidth=0.8, alpha=0.5, zorder=1)
    ax.axvline(-fc_cutoff, linestyle="--", color="black", linewidth=0.8, alpha=0.5, zorder=1)

    # 텍스트 라벨링 (기본 annotate 사용)
    # Up/Down 각각에서 padj가 가장 낮은 상위 10개씩 추출
    top_up = df[df["sig"] == "Up"].nsmallest(10, "padj")
    top_down = df[df["sig"] == "Down"].nsmallest(10, "padj")
    top_genes = pd.concat([top_up, top_down])

    for _, row in top_genes.iterrows():
        ax.annotate(row["gene"],
                    xy=(row["log2FC"], row["neg_log10_padj"]),
                    xytext=(3, 3), # 점으로부터 약간 떨어진 위치에 표시
                    textcoords="offset points",
                    fontsize=8,
                    fontweight='bold',
                    ha="left",
                    va="bottom")

    # 스타일링
    label = fname.replace("MAST_DEGs_", "").replace(".csv", "")
    # ax.set_title(f"Volcano Plot: {label}", fontsize=14, fontweight='bold')
    ax.set_xlabel("log2(Fold Change)", fontsize=12)
    ax.set_ylabel("-log10(adj P-value)", fontsize=12)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.legend(loc="upper right", frameon=True)

    # 카운트 표시 (Up/Down 개수)
    n_up = (df["sig"] == "Up").sum()
    n_down = (df["sig"] == "Down").sum()
    ax.text(0.05, 0.95, f"Up: {n_up}\nDown: {n_down}",
            transform=ax.transAxes, fontsize=10, fontweight='bold',
            verticalalignment="top", bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    plt.tight_layout()
    out_name = f"volcano_{label}.png"
    fig.savefig(os.path.join(out_dir, out_name), dpi=300)
    plt.close()

print("Done! All volcano plots saved successfully.")
