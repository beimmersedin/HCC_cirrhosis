"""
RWR Overview Figure (5 panels)
A: Pipeline schematic
B: Toy network (RWR concept)
C: Real subnetwork (Hepatocyte top genes)
D: Score distribution (Seed vs Non-seed)
E: Cross-cell type top candidates heatmap
"""

import os
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, ArrowStyle
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgba
from matplotlib.cm import ScalarMappable

# ==========================================
# Paths
# ==========================================
BASE_DIR = "/data1/project/yeonu/065_multi_rna/propagation"
RESULT_DIR = os.path.join(BASE_DIR, "results_output")
SEED_DIR = os.path.join(BASE_DIR, "seed")
NETWORK_FILE = os.path.join(BASE_DIR, "string_network.txt")
OUT_DIR = os.path.join(BASE_DIR, "overview")

CELLTYPES = [
    "Hepatocytes", "Endothelial_Cells", "Mesenchymal",
    "T_Cells", "Macrophages", "NK_Cells", "Plasma_Cells"
]

CT_LABELS = {
    "Hepatocytes": "Hepatocytes",
    "Endothelial_Cells": "Endothelial",
    "Mesenchymal": "Mesenchymal",
    "T_Cells": "T Cells",
    "Macrophages": "Macrophages",
    "NK_Cells": "NK Cells",
    "Plasma_Cells": "Plasma Cells"
}

# ==========================================
# Color scheme
# ==========================================
C_SEED = "#E64B35"
C_NONSEED_SIG = "#00A087"
C_OTHER = "#D9D9D9"
C_BLUE = "#3C5488"
C_ORANGE = "#F39B7F"
C_CYAN = "#4DBBD5"

# RWR propagation colormap (red → orange → yellow → light grey)
rwr_cmap = LinearSegmentedColormap.from_list(
    "rwr", ["#E64B35", "#F39B7F", "#FDBE85", "#FEE8C8", "#E0E0E0"]
)

# ==========================================
# Load all cell type results
# ==========================================
all_data = {}
for ct in CELLTYPES:
    path = os.path.join(RESULT_DIR, f"{ct}_NP_Final_Results.csv")
    all_data[ct] = pd.read_csv(path)


# ==========================================
# Figure setup
# ==========================================
fig = plt.figure(figsize=(24, 20))
gs_main = gridspec.GridSpec(2, 1, height_ratios=[1, 0.9], hspace=0.28)

# Top row: A (pipeline) | B (toy) | C (real)
gs_top = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs_main[0],
                                          width_ratios=[0.8, 1, 1.2], wspace=0.25)

# Bottom row: D (score dist) | E (heatmap)
gs_bot = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_main[1],
                                          width_ratios=[1, 1.3], wspace=0.3)

# ==========================================
# Panel A: Pipeline Schematic
# ==========================================
ax_a = fig.add_subplot(gs_top[0])
ax_a.set_xlim(0, 10)
ax_a.set_ylim(0, 10)
ax_a.axis('off')
ax_a.set_title("A", fontsize=18, fontweight='bold', loc='left', x=-0.05, y=1.02)

# Box drawing helper
def draw_pipeline_box(ax, x, y, w, h, color, text, fontsize=10):
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3",
                          facecolor=to_rgba(color, 0.15), edgecolor=color, linewidth=2)
    ax.add_patch(box)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center',
            fontsize=fontsize, fontweight='bold', color=color)

def draw_pipeline_arrow(ax, x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color='#666666', lw=2))

# Pipeline boxes
bw, bh = 8, 1.1
cx = 1.0

draw_pipeline_box(ax_a, cx, 8.3, bw, bh, "#4DBBD5",
                  "WGCNA Co-expression Modules", 11)
draw_pipeline_box(ax_a, cx, 6.6, bw, bh, "#E64B35",
                  "Cirrhosis Up-regulated DEGs\n(MAST: NT vs TC)", 10)
draw_pipeline_box(ax_a, cx, 4.6, bw, bh, "#8491B4",
                  "Seed Genes\n(Module \u2229 DEGs, 7 cell types)", 10)
draw_pipeline_box(ax_a, cx, 2.6, bw, bh, "#3C5488",
                  "STRING PPI Network\n(16,201 nodes, score \u2265 700)", 10)
draw_pipeline_box(ax_a, cx, 0.6, bw, bh, "#00A087",
                  "RWR (r = 0.1)\n\u2192 Candidate Gene Ranking", 11)

# Arrows between boxes
mid_x = cx + bw/2

# Simple vertical flow with ∩ between DEG and Seed
draw_pipeline_arrow(ax_a, mid_x, 8.3, mid_x, 7.75)
draw_pipeline_arrow(ax_a, mid_x, 6.6, mid_x, 6.05)

# ∩ symbol between DEG and Seed boxes
ax_a.text(mid_x, 5.85, "\u2229", ha='center', va='center', fontsize=18,
          fontweight='bold', color='#8491B4',
          bbox=dict(boxstyle='round,pad=0.15', facecolor='white', edgecolor='none'))

draw_pipeline_arrow(ax_a, mid_x, 4.6, mid_x, 3.75)
draw_pipeline_arrow(ax_a, mid_x, 2.6, mid_x, 1.75)


# ==========================================
# Panel B: Toy Network (RWR Concept)
# ==========================================
ax_b = fig.add_subplot(gs_top[1])
ax_b.set_title("B", fontsize=18, fontweight='bold', loc='left', x=-0.05, y=1.02)
ax_b.set_aspect('equal')
ax_b.axis('off')

# Build toy graph
np.random.seed(42)
toy = nx.Graph()

# Seed nodes (center cluster)
seeds = ['S1', 'S2', 'S3']
# 1st neighbors
n1 = ['N1a', 'N1b', 'N1c', 'N1d', 'N1e']
# 2nd neighbors
n2 = ['N2a', 'N2b', 'N2c', 'N2d', 'N2e', 'N2f']
# Far nodes
far = ['F1', 'F2', 'F3', 'F4']

# Edges: seed-seed
toy.add_edges_from([('S1','S2'), ('S2','S3'), ('S1','S3')])
# Edges: seed-N1
toy.add_edges_from([('S1','N1a'), ('S1','N1b'), ('S2','N1c'), ('S2','N1d'), ('S3','N1e')])
# Edges: N1-N1
toy.add_edges_from([('N1a','N1b'), ('N1c','N1d')])
# Edges: N1-N2
toy.add_edges_from([('N1a','N2a'), ('N1b','N2b'), ('N1b','N2c'),
                    ('N1c','N2d'), ('N1d','N2e'), ('N1e','N2f')])
# Edges: N2-N2
toy.add_edges_from([('N2a','N2b'), ('N2d','N2e')])
# Edges: N2-Far
toy.add_edges_from([('N2a','F1'), ('N2c','F2'), ('N2e','F3'), ('N2f','F4')])
# Edges: Far-Far
toy.add_edges_from([('F1','F2'), ('F3','F4')])

# Layout
pos = nx.spring_layout(toy, seed=123, k=1.8, iterations=100)

# Assign scores (simulate RWR)
score_map = {}
for n in seeds:
    score_map[n] = 1.0
for n in n1:
    score_map[n] = 0.65
for n in n2:
    score_map[n] = 0.3
for n in far:
    score_map[n] = 0.05

# Draw edges
nx.draw_networkx_edges(toy, pos, ax=ax_b, edge_color='#CCCCCC', width=1.5, alpha=0.6)

# Draw nodes by category
for node_list, label, size, edgecolor in [
    (far, 'Distant genes', 250, '#AAAAAA'),
    (n2, '2nd neighbors', 350, '#CCCCCC'),
    (n1, '1st neighbors', 450, '#CCCCCC'),
    (seeds, 'Seed genes', 600, 'white'),
]:
    scores = [score_map[n] for n in node_list]
    colors = [rwr_cmap(s) for s in scores]
    nx.draw_networkx_nodes(toy, pos, nodelist=node_list, node_color=colors,
                           node_size=size, edgecolors=edgecolor, linewidths=1.5,
                           ax=ax_b)

# Labels for seed nodes
seed_labels = {n: n for n in seeds}
nx.draw_networkx_labels(toy, pos, labels=seed_labels, font_size=9,
                        font_weight='bold', font_color='white', ax=ax_b)

# Draw "restart" arrow
s1_pos = pos['S1']
ax_b.annotate('', xy=(s1_pos[0]-0.08, s1_pos[1]+0.15),
              xytext=(s1_pos[0]-0.35, s1_pos[1]+0.45),
              arrowprops=dict(arrowstyle='->', color=C_SEED, lw=2.5,
                              connectionstyle='arc3,rad=-0.3'))
ax_b.text(s1_pos[0]-0.55, s1_pos[1]+0.55, 'Restart\n(r = 0.1)',
          fontsize=9, color=C_SEED, fontweight='bold', ha='center')

# Formula
ax_b.text(0.5, -0.08,
          r'$\mathbf{p}^{(t+1)} = (1-r) \cdot W \cdot \mathbf{p}^{(t)} + r \cdot \mathbf{p}_0$',
          transform=ax_b.transAxes, ha='center', va='center', fontsize=13,
          bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#cccccc'))

# Legend
legend_elements = [
    mpatches.Patch(facecolor=rwr_cmap(1.0), edgecolor='grey', label='Seed genes (high score)'),
    mpatches.Patch(facecolor=rwr_cmap(0.65), edgecolor='grey', label='1st neighbors'),
    mpatches.Patch(facecolor=rwr_cmap(0.3), edgecolor='grey', label='2nd neighbors'),
    mpatches.Patch(facecolor=rwr_cmap(0.05), edgecolor='grey', label='Distant genes (low score)'),
]
ax_b.legend(handles=legend_elements, loc='upper left', fontsize=8.5,
            frameon=True, framealpha=0.9, edgecolor='#cccccc')

ax_b.text(0.5, 1.04, "RWR Concept: Signal Propagation",
          transform=ax_b.transAxes, ha='center', fontsize=12, fontweight='bold')


# ==========================================
# Panel C: Real Subnetwork (Hepatocyte Top Genes)
# ==========================================
ax_c = fig.add_subplot(gs_top[2])
ax_c.set_title("C", fontsize=18, fontweight='bold', loc='left', x=-0.05, y=1.02)
ax_c.set_aspect('equal')
ax_c.axis('off')

# Load Hepatocyte results
hep_df = all_data["Hepatocytes"]
top_n = 40
top_genes = hep_df.head(top_n)
top_gene_set = set(top_genes["Gene"])

# Load seed
with open(os.path.join(SEED_DIR, "Hepatocytes_Module1_Cirrhosis_Up_434genes.txt")) as f:
    seed_genes = set(line.strip() for line in f if line.strip())

# Build subnetwork from STRING
G_sub = nx.Graph()
gene_scores = dict(zip(top_genes["Gene"], top_genes["NP_Score"]))
max_score = top_genes["NP_Score"].max()
min_score = top_genes["NP_Score"].min()

with open(NETWORK_FILE) as f:
    for line in f:
        g1, g2, w = line.strip().split('\t')
        if g1 in top_gene_set and g2 in top_gene_set:
            G_sub.add_edge(g1, g2, weight=float(w))

# Add isolated top genes
for g in top_gene_set:
    if g not in G_sub:
        G_sub.add_node(g)

# Layout
pos_c = nx.spring_layout(G_sub, seed=42, k=2.0, iterations=80)

# Normalize scores for coloring
norm = Normalize(vmin=min_score, vmax=max_score)

# Draw edges
nx.draw_networkx_edges(G_sub, pos_c, ax=ax_c, edge_color='#DDDDDD', width=0.8, alpha=0.5)

# Separate seed vs non-seed for drawing
seeds_in_top = [g for g in G_sub.nodes() if g in seed_genes]
nonseeds_in_top = [g for g in G_sub.nodes() if g not in seed_genes]

# Node colors based on NP score
for node_list, marker, edge_c in [(nonseeds_in_top, 'o', C_BLUE), (seeds_in_top, 'o', C_SEED)]:
    if not node_list:
        continue
    scores_norm = [norm(gene_scores.get(n, min_score)) for n in node_list]
    colors = [rwr_cmap(1.0 - s) for s in scores_norm]  # higher score = more red
    sizes = [200 + 600 * norm(gene_scores.get(n, min_score)) for n in node_list]
    nx.draw_networkx_nodes(G_sub, pos_c, nodelist=node_list, node_color=colors,
                           node_size=sizes, edgecolors=edge_c, linewidths=2.0,
                           ax=ax_c)

# Labels for top genes (non-seed top 10 + seed top 3)
top10_nonseed = hep_df[hep_df["is_Seed"] == "No"].head(10)["Gene"].tolist()
top3_seed = hep_df[hep_df["is_Seed"] == "Yes"].head(3)["Gene"].tolist()
label_genes = top10_nonseed + top3_seed
labels_c = {g: g for g in label_genes if g in G_sub}
for node, label in labels_c.items():
    x, y = pos_c[node]
    is_seed = node in seed_genes
    color = C_SEED if is_seed else '#333333'
    ax_c.annotate(label, xy=(x, y), xytext=(8, 8), textcoords='offset points',
                  fontsize=7.5, fontweight='bold', color=color,
                  bbox=dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.8, edgecolor='none'))

# Legend
legend_c = [
    plt.scatter([], [], c=C_SEED, s=100, edgecolors=C_SEED, linewidths=2, label='Seed gene'),
    plt.scatter([], [], c=C_ORANGE, s=100, edgecolors=C_BLUE, linewidths=2, label='Non-seed (high score)'),
    plt.scatter([], [], c='#E0E0E0', s=60, edgecolors=C_BLUE, linewidths=2, label='Non-seed (lower score)'),
]
ax_c.legend(handles=legend_c, loc='upper left', fontsize=8.5,
            frameon=True, framealpha=0.9, edgecolor='#cccccc')

ax_c.text(0.5, 1.04, f"Hepatocyte Top {top_n} Genes — Real PPI Subnetwork",
          transform=ax_c.transAxes, ha='center', fontsize=12, fontweight='bold')


# ==========================================
# Panel D: Score Distribution
# ==========================================
ax_d = fig.add_subplot(gs_bot[0])
ax_d.set_title("D", fontsize=18, fontweight='bold', loc='left', x=-0.05, y=1.02)

df_hep = all_data["Hepatocytes"].copy()
N = len(df_hep)
SCORE_CUTOFF = 1.0 / N

# Separate seed / non-seed
df_seed = df_hep[df_hep["is_Seed"] == "Yes"]
df_nonseed = df_hep[df_hep["is_Seed"] == "No"]

# -log10(P_value)
df_hep["log_p"] = -np.log10(df_hep["P_value"] + 1e-300)

# Histogram of NP scores
bins = np.linspace(0, df_hep["NP_Score"].quantile(0.999), 80)

ax_d.hist(df_nonseed["NP_Score"], bins=bins, color=C_OTHER, alpha=0.7,
          label=f'Non-seed (n={len(df_nonseed):,})', edgecolor='white', linewidth=0.3)
ax_d.hist(df_seed["NP_Score"], bins=bins, color=C_SEED, alpha=0.8,
          label=f'Seed (n={len(df_seed):,})', edgecolor='white', linewidth=0.3)

# 1/N cutoff
ax_d.axvline(SCORE_CUTOFF, color=C_BLUE, linestyle='--', linewidth=2, label=f'1/N = {SCORE_CUTOFF:.2e}')

ax_d.set_xlabel("NP Score", fontsize=12)
ax_d.set_ylabel("Number of Genes", fontsize=12)
ax_d.set_title("    Hepatocyte NP Score Distribution", fontsize=13, fontweight='bold', loc='left')
ax_d.legend(fontsize=10, frameon=True)
ax_d.spines['top'].set_visible(False)
ax_d.spines['right'].set_visible(False)
ax_d.set_xlim(left=0)

# Top gene annotation (after axis limits are set)
top1 = df_hep.iloc[0]
ylim_d = ax_d.get_ylim()
ax_d.annotate(f'{top1["Gene"]}\n(Non-seed, Rank 1)',
              xy=(top1["NP_Score"], 5), xytext=(top1["NP_Score"]*0.6, ylim_d[1]*0.7),
              arrowprops=dict(arrowstyle='->', color=C_BLUE, lw=1.5),
              fontsize=10, fontweight='bold', color=C_BLUE,
              bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=C_BLUE, alpha=0.9))


# ==========================================
# Panel E: Cross-cell Type Heatmap
# ==========================================
ax_e = fig.add_subplot(gs_bot[1])
ax_e.set_title("E", fontsize=18, fontweight='bold', loc='left', x=-0.05, y=1.02)

# Collect top 5 non-seed per cell type
TOP_K = 5
top_genes_per_ct = {}
for ct in CELLTYPES:
    df = all_data[ct]
    nonseed = df[df["is_Seed"] == "No"].nsmallest(TOP_K, "Rank")
    top_genes_per_ct[ct] = nonseed["Gene"].tolist()

# Union of all top genes
all_top_genes = []
seen = set()
for ct in CELLTYPES:
    for g in top_genes_per_ct[ct]:
        if g not in seen:
            all_top_genes.append(g)
            seen.add(g)

# Build rank matrix
rank_matrix = pd.DataFrame(index=all_top_genes, columns=CELLTYPES, dtype=float)
for ct in CELLTYPES:
    df = all_data[ct]
    gene_rank = dict(zip(df["Gene"], df["Rank"]))
    total = len(df)
    for gene in all_top_genes:
        r = gene_rank.get(gene, total)
        rank_matrix.loc[gene, ct] = r

# Percentile rank (lower = better)
pct_matrix = rank_matrix.copy()
for ct in CELLTYPES:
    total = len(all_data[ct])
    pct_matrix[ct] = rank_matrix[ct] / total * 100

# Sort by mean percentile
pct_matrix["mean_pct"] = pct_matrix[CELLTYPES].mean(axis=1)
pct_matrix = pct_matrix.sort_values("mean_pct")
pct_matrix = pct_matrix.drop("mean_pct", axis=1)
rank_matrix = rank_matrix.loc[pct_matrix.index]

# -log10(percentile/100) for display
display_vals = -np.log10(pct_matrix.values.astype(float) / 100 + 1e-5)

ct_display = [CT_LABELS[ct] for ct in CELLTYPES]

im = ax_e.imshow(display_vals, aspect='auto', cmap='YlOrRd', interpolation='nearest')

ax_e.set_xticks(range(len(CELLTYPES)))
ax_e.set_xticklabels(ct_display, fontsize=10, rotation=45, ha='right')
ax_e.set_yticks(range(len(pct_matrix)))
ax_e.set_yticklabels(pct_matrix.index.tolist(), fontsize=9)

# Rank text
for i in range(display_vals.shape[0]):
    for j in range(display_vals.shape[1]):
        rank_val = int(rank_matrix.iloc[i, j])
        color = "white" if display_vals[i, j] > np.percentile(display_vals, 60) else "black"
        ax_e.text(j, i, str(rank_val), ha='center', va='center', fontsize=7.5, color=color)

cbar = plt.colorbar(im, ax=ax_e, shrink=0.6, pad=0.02)
cbar.set_label("-log10(Percentile Rank)", fontsize=10)

# Mark seed status
for i, gene in enumerate(pct_matrix.index):
    for j, ct in enumerate(CELLTYPES):
        df = all_data[ct]
        row = df[df["Gene"] == gene]
        if len(row) > 0 and row.iloc[0]["is_Seed"] == "Yes":
            ax_e.plot(j, i, marker='*', color='white', markersize=6)

ax_e.set_title("    Top Non-seed Genes: Cross-cell Type Rank", fontsize=13, fontweight='bold', loc='left')

# Footnote for star
ax_e.text(1.0, -0.08, "* = seed gene in that cell type",
          transform=ax_e.transAxes, ha='right', fontsize=9, color='#666666', style='italic')


# ==========================================
# Save
# ==========================================
plt.savefig(os.path.join(OUT_DIR, "rwr_overview_figure.png"), dpi=300, bbox_inches='tight',
            facecolor='white')
plt.savefig(os.path.join(OUT_DIR, "rwr_overview_figure.pdf"), bbox_inches='tight',
            facecolor='white')
plt.close()
print("Done! Saved: rwr_overview_figure.png / .pdf")
