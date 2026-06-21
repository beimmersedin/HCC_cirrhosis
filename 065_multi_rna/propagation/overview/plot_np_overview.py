import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# ==========================================
# 데이터 로드
# ==========================================
RESULT_DIR = "/data1/project/yeonu/065_multi_rna/propagation/results_output"
OUT_DIR = "/data1/project/yeonu/065_multi_rna/propagation/overview"

CELLTYPES = [
    "Hepatocytes", "Endothelial_Cells", "Mesenchymal",
    "T_Cells", "Macrophages", "NK_Cells", "Plasma_Cells"
]

all_data = {}
for ct in CELLTYPES:
    path = os.path.join(RESULT_DIR, f"{ct}_NP_Final_Results.csv")
    df = pd.read_csv(path)
    all_data[ct] = df

# ==========================================
# Figure 설정
# ==========================================
fig = plt.figure(figsize=(22, 18))
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1.3], hspace=0.35, wspace=0.3)

# ==========================================
# Panel A: Cell type별 Seed / Non-seed (p<0.05) 카운트
# ==========================================
ax_a = fig.add_subplot(gs[0, 0])

ct_labels = [ct.replace("_", "\n") for ct in CELLTYPES]
seed_counts = []
nonseed_sig_counts = []

for ct in CELLTYPES:
    df = all_data[ct]
    seed_counts.append((df["is_Seed"] == "Yes").sum())
    nonseed_sig_counts.append(((df["is_Seed"] == "No") & (df["P_value"] < 0.05)).sum())

x = np.arange(len(CELLTYPES))
width = 0.35

bars1 = ax_a.bar(x - width/2, seed_counts, width, color="#4DBBD5", label="Seed genes", edgecolor="white")
bars2 = ax_a.bar(x + width/2, nonseed_sig_counts, width, color="#E64B35", label="Non-seed (p < 0.05)", edgecolor="white")

for bar in bars1:
    ax_a.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 8,
              str(int(bar.get_height())), ha='center', va='bottom', fontsize=9, fontweight='bold')
for bar in bars2:
    ax_a.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 8,
              str(int(bar.get_height())), ha='center', va='bottom', fontsize=9, fontweight='bold')

ax_a.set_xticks(x)
ax_a.set_xticklabels(ct_labels, fontsize=10)
ax_a.set_ylabel("Number of genes", fontsize=12)
ax_a.set_title("A. Seed vs. Significant Non-seed Genes", fontsize=13, fontweight='bold', loc='left')
ax_a.legend(frameon=True, fontsize=10)
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)

# ==========================================
# Panel B: Cell type별 Top 1 Non-seed 유전자 NP Score
# ==========================================
ax_b = fig.add_subplot(gs[0, 1])

top1_genes = []
top1_scores = []
for ct in CELLTYPES:
    df = all_data[ct]
    nonseed = df[df["is_Seed"] == "No"].sort_values("Rank")
    top1_genes.append(nonseed.iloc[0]["Gene"])
    top1_scores.append(nonseed.iloc[0]["NP_Score"])

colors_b = ["#E64B35" if s == max(top1_scores) else "#3C5488" for s in top1_scores]
bars = ax_b.barh(ct_labels[::-1], top1_scores[::-1], color=colors_b[::-1], edgecolor="white", height=0.6)

for i, (gene, score) in enumerate(zip(top1_genes[::-1], top1_scores[::-1])):
    ax_b.text(score + max(top1_scores)*0.02, i, f"{gene}", va='center', fontsize=10, fontweight='bold')

ax_b.set_xlabel("NP Score", fontsize=12)
ax_b.set_title("B. Top-ranked Non-seed Gene per Cell Type", fontsize=13, fontweight='bold', loc='left')
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)

# ==========================================
# Panel C: Cross-cell type Heatmap (Top non-seed genes)
# ==========================================
ax_c = fig.add_subplot(gs[1, :])

# 각 cell type에서 non-seed top 10 유전자 수집
TOP_N = 10
top_genes_set = set()
for ct in CELLTYPES:
    df = all_data[ct]
    nonseed = df[df["is_Seed"] == "No"].nsmallest(TOP_N, "Rank")
    top_genes_set.update(nonseed["Gene"].tolist())

top_genes_list = sorted(top_genes_set)

# Rank matrix 생성
rank_matrix = pd.DataFrame(index=top_genes_list, columns=CELLTYPES, dtype=float)
for ct in CELLTYPES:
    df = all_data[ct]
    gene_rank = dict(zip(df["Gene"], df["Rank"]))
    for gene in top_genes_list:
        rank_matrix.loc[gene, ct] = gene_rank.get(gene, np.nan)

# NaN은 큰 값으로 (네트워크에 없는 경우)
rank_matrix = rank_matrix.fillna(rank_matrix.max().max() + 100)

# 평균 rank로 정렬
rank_matrix["mean_rank"] = rank_matrix[CELLTYPES].mean(axis=1)
rank_matrix = rank_matrix.sort_values("mean_rank")
rank_matrix = rank_matrix.drop("mean_rank", axis=1)

# -log10(rank) 변환 (높은 rank = 낮은 값, 시각적으로 직관적)
display_matrix = -np.log10(rank_matrix.values.astype(float))

ct_labels_h = [ct.replace("_", " ") for ct in CELLTYPES]

im = ax_c.imshow(display_matrix, aspect='auto', cmap='YlOrRd', interpolation='nearest')

ax_c.set_xticks(range(len(CELLTYPES)))
ax_c.set_xticklabels(ct_labels_h, fontsize=10, rotation=45, ha='right')
ax_c.set_yticks(range(len(rank_matrix)))
ax_c.set_yticklabels(rank_matrix.index.tolist(), fontsize=8)

# Rank 값 텍스트 표시
for i in range(display_matrix.shape[0]):
    for j in range(display_matrix.shape[1]):
        rank_val = int(rank_matrix.iloc[i, j])
        if rank_val < 16201:  # 실제 rank가 있는 경우만
            color = "white" if display_matrix[i, j] > np.percentile(display_matrix, 70) else "black"
            ax_c.text(j, i, str(rank_val), ha='center', va='center', fontsize=7, color=color)

ax_c.set_title("C. Cross-cell Type Rank Heatmap (Top Non-seed Genes)", fontsize=13, fontweight='bold', loc='left')

# colorbar
cbar = plt.colorbar(im, ax=ax_c, shrink=0.5, pad=0.02)
cbar.set_label("-log10(Rank)", fontsize=10)

plt.savefig(os.path.join(OUT_DIR, "network_propagation_overview.png"), dpi=300, bbox_inches='tight')
plt.close()
print("Done! Saved: network_propagation_overview.png")
