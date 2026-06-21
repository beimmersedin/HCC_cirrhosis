#!/usr/bin/env python
# coding: utf-8
"""
ATAC Analysis Pipeline Overview Figure — 7-Step Version (TF Motif 분리)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 11,
    'figure.facecolor': 'white',
})

# ─── 색상 ────────────────────────────────────────────────────────────────────
C = {
    's1': '#E74C3C',   # DA
    's2': '#E67E22',   # Motif Enrichment
    's3': '#D35400',   # Peak-TF Mapping
    's4': '#F39C12',   # chromVAR
    's5': '#F1C40F',   # Linkage
    's6': '#2ECC71',   # GRN
    's7': '#3498DB',   # Tier
    'hl': '#8E44AD',   # HNF4A
    'ar': '#95A5A6',
    'dk': '#2C3E50',
    'gr': '#7F8C8D',
    'lg': '#BDC3C7',
    'ev': ['#E74C3C', '#E67E22', '#F39C12', '#27AE60', '#3498DB'],
}

fig = plt.figure(figsize=(18, 28))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 18)
ax.set_ylim(0, 28)
ax.axis('off')

# ─── Helpers ─────────────────────────────────────────────────────────────────
def rbox(x, y, w, h, color, alpha=1.0, lw=0, ec='none', zorder=2):
    b = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15",
                       facecolor=color, edgecolor=ec, linewidth=lw,
                       alpha=alpha, zorder=zorder)
    ax.add_patch(b)

def arr(x1, y1, x2, y2, color='#95A5A6', lw=2.5, style='-'):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                                linestyle=style))

# ═══════════════════════════════════════════════════════════════════════════
# TITLE
# ═══════════════════════════════════════════════════════════════════════════
ax.text(9, 27.5, 'ATAC-seq Regulatory Analysis Pipeline',
        fontsize=24, fontweight='bold', ha='center', color=C['dk'])
ax.text(9, 26.95, 'HCC Multiome  |  TC Hepatocytes  |  Cirrhosis vs Hepatitis',
        fontsize=13, ha='center', color=C['gr'])
ax.plot([1, 17], [26.7, 26.7], color=C['lg'], lw=1)

# ═══════════════════════════════════════════════════════════════════════════
# LEFT: 7-Step Pipeline
# ═══════════════════════════════════════════════════════════════════════════
bw, bh = 7.0, 1.85
gap = 0.65  # step 간 간격

steps = [
    {'y': 24.5, 'c': C['s1'], 'n': '01',
     'title': 'Differential Accessibility',
     'method': 'Signac FindMarkers  |  Logistic Regression',
     'detail': 'latent.var = nCount_ATAC  |  padj < 0.05',
     'result': '689 DARs', 'sub': '98.6% ↑ in Cirrhosis',
     'input': 'ATAC peaks × TC Hepatocytes'},
    {'y': 22.0, 'c': C['s2'], 'n': '02',
     'title': 'TF Motif Enrichment',
     'method': 'FindMotifs()  |  Hypergeometric test',
     'detail': 'JASPAR2022 (692 human motifs)  |  FDR < 0.05',
     'result': '291 enriched motifs', 'sub': 'KLF/SP family dominant',
     'input': '687 up-DARs vs background peaks'},
    {'y': 19.5, 'c': C['s3'], 'n': '03',
     'title': 'Peak-TF Binding Map',
     'method': 'motifmatchr  |  PWM sliding window scan',
     'detail': 'p < 5e-05 threshold  |  291 enriched motifs only',
     'result': '35,430 peak-TF pairs', 'sub': '689 peaks × 277 TFs',
     'input': 'DAR sequences + enriched motifs'},
    {'y': 17.0, 'c': C['s4'], 'n': '04',
     'title': 'chromVAR TF Activity',
     'method': 'RunChromVAR()  |  deviation z-score',
     'detail': 'Accessibility bias-corrected per cell per motif',
     'result': '692 × 2,734 matrix', 'sub': 'Cell-level TF activity scores',
     'input': 'Full ATAC matrix + motif-peak map'},
    {'y': 14.5, 'c': C['s5'], 'n': '05',
     'title': 'Peak-to-Gene Linkage',
     'method': 'Signac LinkPeaks  |  Pearson cor, ≤ 500 kb',
     'detail': 'ATAC accessibility ↔ RNA expression per cell',
     'result': '291 peaks → 316 genes', 'sub': '23,571 TF→Peak→Gene triplets',
     'input': 'ATAC + RNA (same cells) + gene annotation'},
    {'y': 12.0, 'c': C['s6'], 'n': '06',
     'title': 'GRN Construction',
     'method': 'Spearman cor  (chromVAR ↔ RNA expression)',
     'detail': 'Cell-level TF activity vs target gene expression',
     'result': '473 edges (86 TFs)', 'sub': 'Spearman ρ > 0.1',
     'input': 'chromVAR scores + RNA counts'},
    {'y': 9.5, 'c': C['s7'], 'n': '07',
     'title': '5-Layer Causal Integration',
     'method': 'Geometric mean of 5 normalized evidence scores',
     'detail': 'DA sig + DA FC + PWM + Linkage + GRN cor',
     'result': '419 Tier 1 triplets', 'sub': '79 Master Regulators',
     'input': 'All upstream results merged'},
]

for i, s in enumerate(steps):
    x0 = 0.5
    y0 = s['y']

    # Box
    rbox(x0, y0, bw, bh, s['c'], alpha=0.07, lw=2, ec=s['c'])

    # Number circle
    circ = plt.Circle((x0 + 0.5, y0 + bh - 0.38), 0.3, color=s['c'], zorder=3)
    ax.add_patch(circ)
    ax.text(x0 + 0.5, y0 + bh - 0.38, s['n'], fontsize=10, fontweight='bold',
            ha='center', va='center', color='white', zorder=4)

    # Title
    ax.text(x0 + 1.05, y0 + bh - 0.38, s['title'],
            fontsize=13, fontweight='bold', color=s['c'], va='center')

    # Method + detail
    ax.text(x0 + 0.3, y0 + bh - 0.85, s['method'],
            fontsize=8.5, color='#666', va='center', style='italic')
    ax.text(x0 + 0.3, y0 + bh - 1.2, s['detail'],
            fontsize=8, color='#999', va='center')

    # Result badge
    rbox(x0 + 0.2, y0 + 0.06, 3.0, 0.4, s['c'], alpha=0.15)
    ax.text(x0 + 1.7, y0 + 0.26, s['result'], fontsize=10, fontweight='bold',
            ha='center', va='center', color=s['c'])

    # Sub
    ax.text(x0 + 5.0, y0 + 0.26, s['sub'], fontsize=8.5,
            ha='center', va='center', color='#888')

    # Arrow to next
    if i < len(steps) - 1:
        ax.annotate('', xy=(x0 + bw/2, s['y'] - 0.3),
                    xytext=(x0 + bw/2, y0),
                    arrowprops=dict(arrowstyle='->', color=C['ar'], lw=2))

# ─── TF Motif bracket (Steps 2-3-4 grouping) ────────────────────────────
bracket_x = 7.7
bracket_top = 23.5
bracket_bot = 17.2
ax.plot([bracket_x, bracket_x + 0.3, bracket_x + 0.3, bracket_x],
        [bracket_top, bracket_top, bracket_bot, bracket_bot],
        color=C['s3'], lw=1.5, alpha=0.5)
ax.text(bracket_x + 0.5, (bracket_top + bracket_bot) / 2,
        'TF Motif\nAnalysis', fontsize=9, fontweight='bold',
        color=C['s3'], va='center', alpha=0.6, rotation=0,
        ha='left', linespacing=1.4)


# ═══════════════════════════════════════════════════════════════════════════
# RIGHT TOP: Causal Evidence Chain
# ═══════════════════════════════════════════════════════════════════════════
rx = 9.2

ax.text(rx + 4, 26.3, 'Causal Evidence Chain', fontsize=17, fontweight='bold',
        ha='center', color=C['dk'])
ax.text(rx + 4, 25.85, 'TF → DAR → Gene  (3-node regulatory path)',
        fontsize=10, ha='center', color='#999')

# 3-node diagram
ny = 24.95
nr = 0.48
nodes = [
    (rx + 1.0, ny, 'TF',   '79 TFs',    C['s1']),
    (rx + 4.0, ny, 'DAR',  '102 DARs',  C['s2']),
    (rx + 7.0, ny, 'Gene', '109 genes', C['s7']),
]
for (nx, ny_, label, sub, col) in nodes:
    circ = plt.Circle((nx, ny_), nr, color=col, zorder=3, alpha=0.9)
    ax.add_patch(circ)
    ax.text(nx, ny_, label, fontsize=12, fontweight='bold',
            ha='center', va='center', color='white', zorder=4)
    ax.text(nx, ny_ - 0.7, sub, fontsize=8.5, ha='center', color='#777')

arr(rx + 1.0 + nr + 0.1, ny, rx + 4.0 - nr - 0.1, ny, C['ar'], 2.5)
arr(rx + 4.0 + nr + 0.1, ny, rx + 7.0 - nr - 0.1, ny, C['ar'], 2.5)
ax.text(rx + 2.5, ny + 0.32, 'binds', fontsize=8, ha='center',
        color='#AAA', style='italic')
ax.text(rx + 5.5, ny + 0.32, 'regulates', fontsize=8, ha='center',
        color='#AAA', style='italic')


# ─── 5 Evidence Layers ──────────────────────────────────────────────────
ey0 = 21.2
eh = 2.65
rbox(rx, ey0, 8, eh, '#F5F5F5', alpha=0.7, lw=1.5, ec='#DDD')
ax.text(rx + 4, ey0 + eh - 0.4, '5 Evidence Layers',
        fontsize=14, fontweight='bold', ha='center', color=C['dk'])

evs = [
    ('DA Significance',     '−log₁₀(padj) / 30',    C['ev'][0], 'Step 01'),
    ('DA Effect Size',      '|log₂FC| / 5',          C['ev'][1], 'Step 01'),
    ('Motif PWM Score',     'PWM score / 20',         C['ev'][2], 'Step 03'),
    ('Peak-Gene Linkage',   '|z-score| / 10',         C['ev'][3], 'Step 05'),
    ('TF-Gene Correlation', '|Spearman ρ| / 0.4',     C['ev'][4], 'Step 06'),
]
for j, (name, formula, col, src) in enumerate(evs):
    ey = ey0 + eh - 0.95 - j * 0.32
    ax.plot(rx + 0.4, ey, 'o', color=col, markersize=7, zorder=3)
    ax.text(rx + 0.75, ey, name, fontsize=9, fontweight='bold',
            va='center', color=col)
    ax.text(rx + 4.3, ey, formula, fontsize=8.5, va='center',
            color='#888', family='monospace')
    ax.text(rx + 7.2, ey, src, fontsize=7, va='center',
            color='#BBB', ha='center')

ax.text(rx + 4, ey0 + 0.22,
        'Composite = ( ev₁ × ev₂ × ev₃ × ev₄ × ev₅ )^(1/5)',
        fontsize=8.5, ha='center', color='#666', family='monospace', style='italic')


# ─── Tier Classification ────────────────────────────────────────────────
ty0 = 19.0
th = 1.8
rbox(rx, ty0, 8, th, '#F5F5F5', alpha=0.7, lw=1.5, ec='#DDD')
ax.text(rx + 4, ty0 + th - 0.4, 'Evidence Tier Classification',
        fontsize=14, fontweight='bold', ha='center', color=C['dk'])

tiers_data = [
    ('Tier 1', 'All 5 evidence layers > 0',  '419',    '1.8%',  C['s1']),
    ('Tier 2', 'PWM + Linkage + DA > 0',     '17,200', '73.0%', '#E8A838'),
    ('Tier 3', 'DA + Linkage > 0',           '5,952',  '25.2%', C['lg']),
]
for j, (tier, criteria, count, pct, col) in enumerate(tiers_data):
    ty = ty0 + th - 0.9 - j * 0.38
    rbox(rx + 0.2, ty - 0.14, 0.95, 0.3, col, alpha=0.9)
    ax.text(rx + 0.67, ty, tier, fontsize=8, fontweight='bold',
            ha='center', va='center', color='white')
    ax.text(rx + 1.4, ty, criteria, fontsize=8.5, va='center', color='#555')
    ax.text(rx + 6.2, ty, count, fontsize=9.5, fontweight='bold',
            va='center', ha='right', color=col)
    ax.text(rx + 7.4, ty, f'({pct})', fontsize=8, va='center',
            ha='right', color='#AAA')


# ─── Dashed arrows (pipeline → evidence) ────────────────────────────────
arr(7.5, 25.0, rx, ey0 + eh - 0.75, C['s1'], 1.0, '--')
arr(7.5, 20.2, rx, ey0 + eh - 1.55, C['s3'], 1.0, '--')
arr(7.5, 15.2, rx, ey0 + eh - 1.9, C['s5'], 1.0, '--')
arr(7.5, 12.7, rx, ey0 + 0.5, C['s6'], 1.0, '--')


# ═══════════════════════════════════════════════════════════════════════════
# RIGHT MIDDLE: Top 10 Master Regulators Table
# ═══════════════════════════════════════════════════════════════════════════
tb_y0 = 12.0
tb_h = 6.5
rbox(rx, tb_y0, 8, tb_h, C['hl'], alpha=0.04, lw=2, ec=C['hl'])

ax.text(rx + 4, tb_y0 + tb_h - 0.45, 'Top 10 Tier 1 Master Regulators',
        fontsize=14, fontweight='bold', ha='center', color=C['hl'])

# Header
hdr_y = tb_y0 + tb_h - 1.05
headers = [('Rank', 0.5), ('TF', 1.4), ('Targets', 4.0), ('DARs', 5.3), ('Score', 6.6)]
for txt, dx in headers:
    ha = 'center' if txt != 'TF' else 'left'
    ax.text(rx + dx, hdr_y, txt, fontsize=9, fontweight='bold', color='#999',
            va='center', ha=ha)
ax.plot([rx + 0.3, rx + 7.7], [hdr_y - 0.2, hdr_y - 0.2], color='#DDD', lw=0.8)

top10 = [
    ('1',  'HNF4A',         '55', '43', '0.301'),
    ('2',  'NR1H2::RXRA',   '31', '18', '0.248'),
    ('3',  'RXRG',          '28', '23', '0.229'),
    ('4',  'NR2F1',         '22', '19', '0.260'),
    ('5',  'RXRB',          '22', '16', '0.221'),
    ('6',  'HNF4G',         '21', '17', '0.321'),
    ('7',  'CTCF',          '19', '17', '0.306'),
    ('8',  'HNF1B',         '16', '18', '0.277'),
    ('9',  'TEAD4',         '13', '10', '0.255'),
    ('10', 'HNF1A',         '11', '14', '0.251'),
]

for j, (rank, tf, tgt, dars, score) in enumerate(top10):
    ry = hdr_y - 0.6 - j * 0.45
    is_top = rank == '1'
    tf_color = C['hl'] if is_top else C['dk']
    tf_weight = 'bold' if is_top else 'normal'
    num_color = C['hl'] if is_top else '#666'

    if is_top:
        rbox(rx + 0.15, ry - 0.18, 7.6, 0.4, C['hl'], alpha=0.08)

    ax.text(rx + 0.5, ry, rank, fontsize=9, ha='center', va='center', color='#BBB')
    ax.text(rx + 1.4, ry, tf, fontsize=10, fontweight=tf_weight,
            va='center', color=tf_color)
    ax.text(rx + 4.0, ry, tgt, fontsize=10, fontweight=tf_weight,
            ha='center', va='center', color=num_color)
    ax.text(rx + 5.3, ry, dars, fontsize=10, ha='center', va='center', color=num_color)
    ax.text(rx + 6.6, ry, score, fontsize=9, ha='center', va='center', color=num_color)


# ═══════════════════════════════════════════════════════════════════════════
# BOTTOM: HNF4A Key Finding
# ═══════════════════════════════════════════════════════════════════════════
hy0 = 0.5
hh = 8.8
rbox(0.5, hy0, 17, hh, C['hl'], alpha=0.05, lw=2.5, ec=C['hl'])

ax.text(9, hy0 + hh - 0.55, 'Key Finding',
        fontsize=20, fontweight='bold', ha='center', color=C['hl'])
ax.plot([3, 15], [hy0 + hh - 0.95, hy0 + hh - 0.95], color=C['hl'],
        lw=0.8, alpha=0.4)

# HNF4A title
ax.text(1.5, hy0 + hh - 1.8, 'HNF4A', fontsize=40, fontweight='bold',
        color=C['hl'], va='center')
ax.text(1.5, hy0 + hh - 2.6, 'Hepatocyte Nuclear Factor 4 Alpha',
        fontsize=12, color='#999', va='center', style='italic')
ax.text(1.5, hy0 + hh - 3.2, '#1 Master Regulator of chromatin remodeling',
        fontsize=11, color='#777', va='center')
ax.text(1.5, hy0 + hh - 3.7, 'in Cirrhosis-associated Hepatocytes',
        fontsize=11, color='#777', va='center')

# Stat boxes
stat_items = [
    ('55', 'target genes'),
    ('43', 'mediating DARs'),
    ('0.301', 'causal score'),
]
for k, (val, label) in enumerate(stat_items):
    sx = 1.3 + k * 3.2
    sy = hy0 + 0.6
    rbox(sx, sy, 2.7, 2.3, C['hl'], alpha=0.08, lw=1.5, ec=C['hl'])
    ax.text(sx + 1.35, sy + 1.5, val, fontsize=30, fontweight='bold',
            ha='center', va='center', color=C['hl'])
    ax.text(sx + 1.35, sy + 0.5, label, fontsize=11,
            ha='center', va='center', color='#888')

# Notable target genes
ntx = 11.0
ax.text(ntx, hy0 + hh - 1.8, 'Notable HNF4A Target Genes',
        fontsize=14, fontweight='bold', color=C['hl'], va='center')

targets = [
    ('CYP3A4',  'Drug metabolism (cytochrome P450)'),
    ('AGT',     'Angiotensinogen — blood pressure'),
    ('GPT',     'ALT enzyme — liver function marker'),
    ('SEC14L2', 'Lipid transport & cholesterol'),
    ('TF',      'Transferrin — iron transport'),
    ('POR',     'P450 oxidoreductase'),
    ('SLC27A5', 'Bile acid metabolism'),
]

for k, (gene, func) in enumerate(targets):
    ty = hy0 + hh - 2.5 - k * 0.55
    rbox(ntx, ty - 0.18, 1.6, 0.38, C['hl'], alpha=0.12)
    ax.text(ntx + 0.8, ty, gene, fontsize=10, fontweight='bold',
            ha='center', va='center', color=C['hl'])
    ax.text(ntx + 1.9, ty, func, fontsize=9, va='center', color='#888')


# ═══════════════════════════════════════════════════════════════════════════
# SAVE
# ═══════════════════════════════════════════════════════════════════════════
save_dir = '/data1/project/yeonu/065_multi_atac/analysis/figures'
for fmt in ['png', 'pdf']:
    path = f'{save_dir}/RE_pipeline_overview.{fmt}'
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f'Saved: {path}')

plt.close()
print('Done!')
