#!/usr/bin/env python
# coding: utf-8
"""
TF Motif Analysis — Clean Overview Diagram for PPT
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

fig, ax = plt.subplots(figsize=(18, 11))
ax.set_xlim(0, 18)
ax.set_ylim(0, 11)
ax.axis('off')
fig.patch.set_facecolor('white')

# ── Colors ──
C_BG1   = '#EEF5FF'   # Step1 bg
C_BD1   = '#4A90D9'
C_BG2   = '#FFF5EE'   # Step2 bg
C_BD2   = '#E8834A'
C_BG3   = '#F0FFF0'   # Step3 bg
C_BD3   = '#5BAA5B'
C_INPUT = '#F7F7F7'
C_NUM   = '#C0392B'    # emphasis numbers
C_TXT   = '#2C3E50'
C_SUB   = '#7F8C8D'
C_ARROW = '#34495E'

def rounded_box(ax, x, y, w, h, fc, ec, lw=2, alpha=1.0):
    box = FancyBboxPatch((x, y), w, h,
                         boxstyle="round,pad=0.2,rounding_size=0.3",
                         facecolor=fc, edgecolor=ec, linewidth=lw, alpha=alpha, zorder=1)
    ax.add_patch(box)

def small_box(ax, x, y, w, h, fc, ec, lw=1.2):
    box = FancyBboxPatch((x, y), w, h,
                         boxstyle="round,pad=0.1,rounding_size=0.15",
                         facecolor=fc, edgecolor=ec, linewidth=lw, zorder=2)
    ax.add_patch(box)

def arrow_down(ax, x, y1, y2):
    ax.annotate('', xy=(x, y2), xytext=(x, y1),
                arrowprops=dict(arrowstyle='->', color=C_ARROW, lw=2.5), zorder=3)

def arrow_right(ax, x1, y, x2):
    ax.annotate('', xy=(x2, y), xytext=(x1, y),
                arrowprops=dict(arrowstyle='->', color=C_ARROW, lw=2.5), zorder=3)

# ════════════════════════════════════════════════════════════════════
# Title
# ════════════════════════════════════════════════════════════════════
ax.text(9, 10.5, 'TF Motif Analysis Pipeline', ha='center', fontsize=22,
        fontweight='bold', color=C_TXT)
ax.text(9, 10.05, 'TC Hepatocytes  |  Cirrhosis vs Hepatitis  |  JASPAR2022',
        ha='center', fontsize=11, color=C_SUB)

# ════════════════════════════════════════════════════════════════════
# Input: DAR peaks (top-left)
# ════════════════════════════════════════════════════════════════════
small_box(ax, 0.3, 8.4, 3.0, 1.2, '#FFFFFF', '#95A5A6', lw=1.5)
ax.text(1.8, 9.2, 'Input', ha='center', fontsize=10, color=C_SUB, fontweight='bold')
ax.text(1.8, 8.75, '689 Significant DARs', ha='center', fontsize=12,
        fontweight='bold', color=C_TXT)

# Input: JASPAR (top-right area for step1)
small_box(ax, 4.5, 8.4, 3.0, 1.2, '#FFFFFF', '#95A5A6', lw=1.5)
ax.text(6.0, 9.2, 'Motif DB', ha='center', fontsize=10, color=C_SUB, fontweight='bold')
ax.text(6.0, 8.75, '692 TF Motifs', ha='center', fontsize=12,
        fontweight='bold', color=C_TXT)
ax.text(6.0, 8.45, '(JASPAR2022, Human)', ha='center', fontsize=8.5, color=C_SUB)

# Arrows from inputs to Step1
arrow_down(ax, 1.8, 8.4, 7.75)
arrow_down(ax, 6.0, 8.4, 7.75)

# ════════════════════════════════════════════════════════════════════
# STEP 1 — Motif Enrichment
# ════════════════════════════════════════════════════════════════════
rounded_box(ax, 0.2, 5.6, 7.4, 2.15, C_BG1, C_BD1)

# Step label
ax.text(0.6, 7.35, 'Step 1', fontsize=13, fontweight='bold', color=C_BD1)
ax.text(2.3, 7.35, 'Motif Enrichment', fontsize=13, fontweight='bold', color=C_TXT)

# Content
ax.text(0.6, 6.85, 'FindMotifs()  —  Hypergeometric test', fontsize=10, color=C_TXT)
ax.text(0.6, 6.4, 'Motif frequency in DARs  vs  all peaks (background)',
        fontsize=9.5, color=C_SUB)

# Example box
small_box(ax, 0.5, 5.75, 6.8, 0.55, 'white', C_BD1, lw=1)
ax.text(0.7, 6.0, 'e.g. KLF15:  DAR ', fontsize=9, color=C_TXT, family='monospace')
ax.text(3.2, 6.0, '57.2%', fontsize=9, color=C_NUM, fontweight='bold', family='monospace')
ax.text(4.2, 6.0, ' vs  Background ', fontsize=9, color=C_TXT, family='monospace')
ax.text(6.15, 6.0, '24.7%', fontsize=9, color=C_NUM, fontweight='bold', family='monospace')
ax.text(7.0, 6.0, '  (2.3x)', fontsize=9, color=C_NUM, fontweight='bold', family='monospace')

# ── Step 1 output arrow ──
arrow_right(ax, 7.6, 6.7, 8.3)

# Step 1 output
small_box(ax, 8.3, 6.1, 3.8, 1.3, 'white', C_BD1, lw=1.5)
ax.text(10.2, 7.05, '291 Enriched Motifs', ha='center', fontsize=12,
        fontweight='bold', color=C_BD1)
ax.text(10.2, 6.6, '(FDR < 0.05)', ha='center', fontsize=10, color=C_SUB)
ax.text(10.2, 6.25, '692 > 291  (42% passed)', ha='center', fontsize=9, color=C_NUM)

# ════════════════════════════════════════════════════════════════════
# Arrow Step1 → Step2
# ════════════════════════════════════════════════════════════════════
arrow_down(ax, 3.9, 5.6, 4.95)
arrow_down(ax, 10.2, 6.1, 4.95)

# ════════════════════════════════════════════════════════════════════
# STEP 2 — Peak-TF Binding
# ════════════════════════════════════════════════════════════════════
rounded_box(ax, 0.2, 2.8, 7.4, 2.15, C_BG2, C_BD2)

ax.text(0.6, 4.55, 'Step 2', fontsize=13, fontweight='bold', color=C_BD2)
ax.text(2.3, 4.55, 'Peak-TF Binding Map', fontsize=13, fontweight='bold', color=C_TXT)

ax.text(0.6, 4.05, 'motifmatchr  —  PWM sliding window scan', fontsize=10, color=C_TXT)
ax.text(0.6, 3.6, 'Scan 291 motif PWMs across 689 DAR sequences (p < 5e-05)',
        fontsize=9.5, color=C_SUB)

# Illustration of scanning
small_box(ax, 0.5, 2.95, 6.8, 0.55, 'white', C_BD2, lw=1)
ax.text(0.7, 3.2, 'Peak (~1,000bp):', fontsize=9, color=C_TXT)
ax.text(2.65, 3.2, '|████████████████████|', fontsize=9, color='#BDC3C7', family='monospace')
ax.text(5.0, 3.2, '← [TF motif ~10bp] scan →', fontsize=8.5, color=C_BD2, family='monospace')

# Step 2 output
arrow_right(ax, 7.6, 3.9, 8.3)

small_box(ax, 8.3, 3.3, 3.8, 1.3, 'white', C_BD2, lw=1.5)
ax.text(10.2, 4.25, '35,430 Binding Pairs', ha='center', fontsize=12,
        fontweight='bold', color=C_BD2)
ax.text(10.2, 3.8, '689 peaks × 277 TFs', ha='center', fontsize=10, color=C_SUB)
ax.text(10.2, 3.45, '+ PWM affinity score', ha='center', fontsize=9, color=C_SUB)

# ════════════════════════════════════════════════════════════════════
# Arrow Step2 → Step3
# ════════════════════════════════════════════════════════════════════
arrow_down(ax, 3.9, 2.8, 2.15)
arrow_down(ax, 10.2, 3.3, 2.15)

# ════════════════════════════════════════════════════════════════════
# STEP 3 — chromVAR
# ════════════════════════════════════════════════════════════════════
rounded_box(ax, 0.2, 0.2, 7.4, 1.95, C_BG3, C_BD3)

ax.text(0.6, 1.75, 'Step 3', fontsize=13, fontweight='bold', color=C_BD3)
ax.text(2.3, 1.75, 'Cell-level TF Activity', fontsize=13, fontweight='bold', color=C_TXT)

ax.text(0.6, 1.3, 'chromVAR  —  Per-cell TF activity deviation z-scores',
        fontsize=10, color=C_TXT)
ax.text(0.6, 0.85, 'Quantify how accessible each TF binding site is in each cell',
        fontsize=9.5, color=C_SUB)
ax.text(0.6, 0.45, '> Downstream: GRN construction, Causal Evidence integration',
        fontsize=9.5, color=C_BD3, fontweight='bold')

# Step 3 output
arrow_right(ax, 7.6, 1.2, 8.3)

small_box(ax, 8.3, 0.55, 3.8, 1.3, 'white', C_BD3, lw=1.5)
ax.text(10.2, 1.5, 'TF Activity Matrix', ha='center', fontsize=12,
        fontweight='bold', color=C_BD3)
ax.text(10.2, 1.05, '692 motifs × 2,734 cells', ha='center', fontsize=10, color=C_SUB)
ax.text(10.2, 0.7, 'z-score per cell per TF', ha='center', fontsize=9, color=C_SUB)

# ════════════════════════════════════════════════════════════════════
# Right side: Funnel summary
# ════════════════════════════════════════════════════════════════════
# Funnel visualization
fx = 14.5
rounded_box(ax, 12.8, 0.2, 4.8, 9.35, '#FAFAFA', '#D5D8DC', lw=1.5)
ax.text(fx, 9.15, 'Filtering Funnel', ha='center', fontsize=13,
        fontweight='bold', color=C_TXT)

# Funnel bars
bar_data = [
    (8.3, 3.8, '692', 'Total Motifs (JASPAR2022)', '#D6EAF8', C_BD1),
    (7.1, 2.6, '291', 'Enriched (FDR < 0.05)', C_BG1, C_BD1),
    (5.6, 2.0, '277', 'TFs with binding', C_BG2, C_BD2),
    (4.2, 1.6, '79', 'Tier 1 Master Regulators', C_BG3, C_BD3),
]

for y, bw, num, label, fc, ec in bar_data:
    bx = fx - bw/2
    small_box(ax, bx, y, bw, 0.7, fc, ec, lw=1.5)
    ax.text(fx, y + 0.45, num, ha='center', fontsize=16, fontweight='bold', color=ec)
    ax.text(fx, y + 0.1, label, ha='center', fontsize=8.5, color=C_SUB)

# Arrows between funnel bars
for i in range(len(bar_data)-1):
    y_from = bar_data[i][0]
    y_to = bar_data[i+1][0] + 0.7
    arrow_down(ax, fx, y_from, y_to)

# Percentage labels on arrows
ax.text(fx + 2.0, 7.85, '42%', fontsize=9, color=C_NUM, fontweight='bold', ha='center')
ax.text(fx + 2.0, 7.55, 'pass', fontsize=8, color=C_SUB, ha='center')

ax.text(fx + 1.6, 6.5, '95%', fontsize=9, color=C_NUM, fontweight='bold', ha='center')
ax.text(fx + 1.6, 6.2, 'mapped', fontsize=8, color=C_SUB, ha='center')

ax.text(fx + 1.2, 5.1, '29%', fontsize=9, color=C_NUM, fontweight='bold', ha='center')
ax.text(fx + 1.2, 4.8, 'causal', fontsize=8, color=C_SUB, ha='center')

# ════════════════════════════════════════════════════════════════════
# Save
# ════════════════════════════════════════════════════════════════════
out_path = '/data1/project/yeonu/065_multi_atac/TFmotif/results_output/TF_pipeline_overview.png'
plt.savefig(out_path, dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {out_path}')
