import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

fig, ax = plt.subplots(figsize=(18, 22))
ax.set_xlim(0, 18)
ax.set_ylim(0, 22)
ax.axis('off')

# ==========================================
# Color palette
# ==========================================
C_DB = "#3C5488"       # database blue
C_FILTER = "#4DBBD5"   # filter cyan
C_PROCESS = "#00A087"  # process green
C_SEED = "#E64B35"     # seed red
C_RWR = "#F39B7F"      # RWR orange
C_RESULT = "#8491B4"   # result purple
C_TEXT = "#2D2D2D"

def draw_box(ax, x, y, w, h, color, title, contents, title_size=13, content_size=10.5, alpha=0.15):
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3",
                          facecolor=color, edgecolor=color, alpha=alpha, linewidth=2)
    ax.add_patch(box)
    border = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3",
                             facecolor='none', edgecolor=color, linewidth=2.5)
    ax.add_patch(border)
    ax.text(x + w/2, y + h - 0.35, title, ha='center', va='top',
            fontsize=title_size, fontweight='bold', color=color)
    for i, line in enumerate(contents):
        ax.text(x + w/2, y + h - 0.85 - i*0.42, line, ha='center', va='top',
                fontsize=content_size, color=C_TEXT)

def draw_arrow(ax, x1, y1, x2, y2, color="#888888"):
    arrow = FancyArrowPatch((x1, y1), (x2, y2),
                             arrowstyle='->', mutation_scale=25,
                             linewidth=2.5, color=color, zorder=5)
    ax.add_patch(arrow)

def draw_step_label(ax, x, y, step_num):
    circle = plt.Circle((x, y), 0.35, color="#2D2D2D", zorder=10)
    ax.add_patch(circle)
    ax.text(x, y, str(step_num), ha='center', va='center',
            fontsize=14, fontweight='bold', color='white', zorder=11)

# ==========================================
# Title
# ==========================================
ax.text(9, 21.5, "Network Propagation Pipeline", ha='center', va='center',
        fontsize=22, fontweight='bold', color=C_TEXT)
ax.text(9, 21.0, "Random Walk with Restart (RWR) on STRING PPI Network",
        ha='center', va='center', fontsize=13, color='#666666')

# ==========================================
# Step 1: STRING DB
# ==========================================
y1 = 18.8
draw_step_label(ax, 1.5, y1 + 0.7, 1)
draw_box(ax, 2.2, y1 - 0.3, 5.5, 2.0, C_DB,
         "STRING Database v12.0",
         ["Protein-Protein Interaction (PPI)",
          "9606.protein.links.v12.0.txt",
          "9606.protein.info.v12.0.txt"])

draw_box(ax, 10.3, y1 - 0.3, 5.5, 2.0, C_DB,
         "ID Mapping (create_network.py)",
         ["ENSP protein ID  >  Gene Symbol",
          "Using preferred_name field",
          "e.g. ENSP00000269305 > TP53"])

draw_arrow(ax, 7.7, y1 + 0.7, 10.3, y1 + 0.7, C_DB)

# ==========================================
# Step 2: Filtering
# ==========================================
y2 = 16.3
draw_step_label(ax, 1.5, y2 + 0.7, 2)
draw_arrow(ax, 5, y1 - 0.3, 5, y2 + 1.7, "#888888")

draw_box(ax, 2.2, y2 - 0.3, 5.5, 2.0, C_FILTER,
         "High Confidence Filtering",
         ["Combined score >= 700 (high conf.)",
          "Remove low-confidence edges",
          "> string_network.txt"])

draw_box(ax, 10.3, y2 - 0.3, 5.5, 2.0, C_FILTER,
         "Filtered Network Statistics",
         ["Nodes (genes): 16,201",
          "Edges (interactions): 473,860",
          "Format: Gene1  Gene2  Score/1000"])

draw_arrow(ax, 7.7, y2 + 0.7, 10.3, y2 + 0.7, C_FILTER)

# ==========================================
# Step 3: Seed Genes
# ==========================================
y3 = 13.5
draw_step_label(ax, 1.5, y3 + 0.85, 3)
draw_arrow(ax, 5, y2 - 0.3, 5, y3 + 1.85, "#888888")

draw_box(ax, 2.2, y3 - 0.45, 5.5, 2.3, C_SEED,
         "Seed Gene Preparation",
         ["WGCNA Module (co-expression)",
          "  Cirrhosis Up DEGs (MAST)",
          "= Cell type-specific seed genes",
          "7 cell types x different seed sets"])

draw_box(ax, 10.3, y3 - 0.45, 5.5, 2.3, C_SEED,
         "Seed Genes per Cell Type",
         ["Endothelial: 511  |  Hepatocytes: 434",
          "Mesenchymal: 391  |  T Cells: 313",
          "Macrophages: 133  |  NK Cells: 116",
          "Plasma Cells: 18"])

draw_arrow(ax, 7.7, y3 + 0.7, 10.3, y3 + 0.7, C_SEED)

# ==========================================
# Step 4: RWR
# ==========================================
y4 = 10.2
draw_step_label(ax, 1.5, y4 + 1.1, 4)
draw_arrow(ax, 5, y3 - 0.45, 5, y4 + 2.1, "#888888")

draw_box(ax, 2.2, y4 - 0.7, 13.6, 2.8, C_RWR,
         "Random Walk with Restart (RWR)",
         [])

# RWR formula
ax.text(9, y4 + 1.35, r"$\mathbf{p}^{(t+1)} = (1 - r) \cdot W \cdot \mathbf{p}^{(t)} + r \cdot \mathbf{p}_0$",
        ha='center', va='center', fontsize=16, color=C_TEXT,
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#cccccc', linewidth=1.5))

# RWR details - left
ax.text(5.5, y4 + 0.5, "W : column-normalized adjacency matrix",
        ha='center', va='center', fontsize=10.5, color=C_TEXT)
ax.text(5.5, y4 + 0.1, "r = 0.1 : restart probability",
        ha='center', va='center', fontsize=10.5, color=C_TEXT)
ax.text(5.5, y4 - 0.3, "Convergence: L1 norm < 1e-6",
        ha='center', va='center', fontsize=10.5, color=C_TEXT)

# RWR details - right
ax.text(12.5, y4 + 0.5, "p0 : initial probability (seed genes = 1/N)",
        ha='center', va='center', fontsize=10.5, color=C_TEXT)
ax.text(12.5, y4 + 0.1, "constantWeight = True (uniform seed weight)",
        ha='center', va='center', fontsize=10.5, color=C_TEXT)
ax.text(12.5, y4 - 0.3, "normalize = True (sum of p0 = 1)",
        ha='center', va='center', fontsize=10.5, color=C_TEXT)

# ==========================================
# Step 5: Convergence
# ==========================================
y5 = 7.2
draw_step_label(ax, 1.5, y5 + 0.85, 5)
draw_arrow(ax, 9, y4 - 0.7, 9, y5 + 1.85, "#888888")

draw_box(ax, 2.2, y5 - 0.45, 5.5, 2.3, C_PROCESS,
         "Iterative Propagation",
         ["Propagate probability to neighbors",
          "10% restart to seed at each step",
          "Iterate until convergence",
          "> Final probability = NP Score"])

draw_box(ax, 10.3, y5 - 0.45, 5.5, 2.3, C_PROCESS,
         "Interpretation",
         ["High NP Score =",
          "  Close to seeds in the network",
          "  Signal converges from multiple seeds",
          "> Disease-associated candidate gene"])

draw_arrow(ax, 7.7, y5 + 0.7, 10.3, y5 + 0.7, C_PROCESS)

# ==========================================
# Step 6: Output
# ==========================================
y6 = 4.3
draw_step_label(ax, 1.5, y6 + 0.85, 6)
draw_arrow(ax, 5, y5 - 0.45, 5, y6 + 1.85, "#888888")

draw_box(ax, 2.2, y6 - 0.45, 5.5, 2.3, C_RESULT,
         "Post-processing",
         ["Rank by NP Score",
          "Label Seed vs Non-seed",
          "Calculate P-value",
          "> {CellType}_NP_Final_Results.csv"])

draw_box(ax, 10.3, y6 - 0.45, 5.5, 2.3, C_RESULT,
         "Key Output",
         ["Non-seed + significant genes =",
          "  Novel candidates from propagation",
          "e.g. TP53 (Rank 1 in 4/7 cell types)",
          "     CTNNB1, EGFR, EP300, etc."])

draw_arrow(ax, 7.7, y6 + 0.7, 10.3, y6 + 0.7, C_RESULT)

# ==========================================
# Bottom: Reference
# ==========================================
ax.text(9, 3.2, "References: Kohler et al., AJHG (2008); Szklarczyk et al., NAR (2023); Mohsen et al., Genome Biology (2021)",
        ha='center', va='center', fontsize=9.5, color='#999999', style='italic')

plt.savefig("/data1/project/yeonu/065_multi_rna/propagation/overview/network_propagation_pipeline.png",
            dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("Done! Saved: network_propagation_pipeline.png")
