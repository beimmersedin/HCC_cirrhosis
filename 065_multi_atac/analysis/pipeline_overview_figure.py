#!/usr/bin/env python
# coding: utf-8
"""
ATAC Analysis Pipeline Overview Figure — PPT용 (v2)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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
    'step1':     '#E74C3C',
    'step2':     '#E67E22',
    'step3':     '#F1C40F',
    'step4':     '#2ECC71',
    'step5':     '#3498DB',
    'highlight': '#8E44AD',
    'arrow':     '#95A5A6',
    'dark':      '#2C3E50',
    'grey':      '#7F8C8D',
    'lightgrey': '#BDC3C7',
    'ev': ['#E74C3C', '#E67E22', '#F39C12', '#27AE60', '#3498DB'],
}

fig = plt.figure(figsize=(18, 22))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 18)
ax.set_ylim(0, 22)
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
ax.text(9, 21.5, 'ATAC-seq Regulatory Analysis Pipeline',
        fontsize=24, fontweight='bold', ha='center', color=C['dark'])
ax.text(9, 20.95, 'HCC Multiome  |  TC Hepatocytes  |  Cirrhosis vs Hepatitis',
        fontsize=13, ha='center', color=C['grey'])

# thin separator line
ax.plot([1, 17], [20.7, 20.7], color=C['lightgrey'], lw=1)

# ═══════════════════════════════════════════════════════════════════════════
# LEFT: Pipeline Steps (x: 0.5 ~ 7.5)
# ═══════════════════════════════════════════════════════════════════════════
bw, bh = 7.0, 2.0
steps = [
    {'y': 18.2, 'c': C['step1'], 'n': '01',
     'title': 'Differential Accessibility',
     'method': 'Signac FindMarkers  |  Logistic Regression',
     'detail': 'latent.var = nCount_ATAC  |  padj < 0.05',
     'result': '689 DARs', 'sub': '98.6% ↑ in Cirrhosis'},
    {'y': 15.6, 'c': C['step2'], 'n': '02',
     'title': 'TF Motif Enrichment',
     'method': 'JASPAR2022  |  motifmatchr  |  chromVAR',
     'detail': '692 human motifs scanned  |  FDR < 0.05',
     'result': '291 motifs, 277 TFs', 'sub': '35,430 peak-TF pairs'},
    {'y': 13.0, 'c': C['step3'], 'n': '03',
     'title': 'Peak-to-Gene Linkage',
     'method': 'Signac LinkPeaks  |  500 kb window',
     'detail': 'ATAC accessibility ↔ RNA expression correlation',
     'result': '291 peaks → 316 genes', 'sub': '23,571 TF→Peak→Gene triplets'},
    {'y': 10.4, 'c': C['step4'], 'n': '04',
     'title': 'GRN Construction',
     'method': 'Spearman cor  (chromVAR TF activity ↔ RNA)',
     'detail': 'Cell-level TF activity vs gene expression',
     'result': '473 edges (86 TFs)', 'sub': 'Spearman ρ > 0.1'},
    {'y': 7.8, 'c': C['step5'], 'n': '05',
     'title': '5-Layer Causal Integration',
     'method': 'Geometric mean of 5 normalized evidence scores',
     'detail': 'DA sig + DA FC + PWM + Linkage + GRN cor',
     'result': '419 Tier 1 triplets', 'sub': '79 Master Regulators'},
]

for i, s in enumerate(steps):
    x0 = 0.5
    y0 = s['y']

    rbox(x0, y0, bw, bh, s['c'], alpha=0.07, lw=2, ec=s['c'])

    # number circle
    circ = plt.Circle((x0 + 0.55, y0 + bh - 0.42), 0.32, color=s['c'], zorder=3)
    ax.add_patch(circ)
    ax.text(x0 + 0.55, y0 + bh - 0.42, s['n'], fontsize=11, fontweight='bold',
            ha='center', va='center', color='white', zorder=4)

    # title
    ax.text(x0 + 1.15, y0 + bh - 0.42, s['title'],
            fontsize=14, fontweight='bold', color=s['c'], va='center')

    # method + detail
    ax.text(x0 + 0.35, y0 + bh - 0.95, s['method'],
            fontsize=9, color='#666', va='center', style='italic')
    ax.text(x0 + 0.35, y0 + bh - 1.35, s['detail'],
            fontsize=8.5, color='#999', va='center')

    # result badge
    rbox(x0 + 0.25, y0 + 0.08, 3.3, 0.45, s['c'], alpha=0.15)
    ax.text(x0 + 1.9, y0 + 0.3, s['result'], fontsize=11, fontweight='bold',
            ha='center', va='center', color=s['c'])

    # sub
    ax.text(x0 + 5.2, y0 + 0.3, s['sub'], fontsize=9,
            ha='center', va='center', color='#888')

    # arrow down
    if i < len(steps) - 1:
        ax.annotate('', xy=(x0 + bw/2, s['y'] - 0.35),
                    xytext=(x0 + bw/2, y0),
                    arrowprops=dict(arrowstyle='->', color=C['arrow'], lw=2))


# ═══════════════════════════════════════════════════════════════════════════
# RIGHT TOP: Causal Evidence Chain (x: 8.5 ~ 17)
# ═══════════════════════════════════════════════════════════════════════════
rx = 8.8

# --- Section: 3-node diagram ---
ax.text(rx + 4, 20.3, 'Causal Evidence Chain', fontsize=17, fontweight='bold',
        ha='center', color=C['dark'])
ax.text(rx + 4, 19.85, 'TF → DAR → Gene  (3-node regulatory path)',
        fontsize=10, ha='center', color='#999')

ny = 18.85
nr = 0.5
positions = [
    (rx + 1.0, ny, 'TF',   '79 TFs',    C['step1']),
    (rx + 4.0, ny, 'DAR',  '102 DARs',  C['step2']),
    (rx + 7.0, ny, 'Gene', '109 genes', C['step5']),
]
for (nx, ny_, label, sub, col) in positions:
    circ = plt.Circle((nx, ny_), nr, color=col, zorder=3, alpha=0.9)
    ax.add_patch(circ)
    ax.text(nx, ny_, label, fontsize=13, fontweight='bold',
            ha='center', va='center', color='white', zorder=4)
    ax.text(nx, ny_ - 0.75, sub, fontsize=9, ha='center', color='#777')

# arrows
arr(rx + 1.0 + nr + 0.1, ny, rx + 4.0 - nr - 0.1, ny, C['arrow'], 2.5)
arr(rx + 4.0 + nr + 0.1, ny, rx + 7.0 - nr - 0.1, ny, C['arrow'], 2.5)
ax.text(rx + 2.5, ny + 0.35, 'binds', fontsize=8, ha='center',
        color='#AAA', style='italic')
ax.text(rx + 5.5, ny + 0.35, 'regulates', fontsize=8, ha='center',
        color='#AAA', style='italic')

# --- Section: 5 Evidence Layers ---
ey0 = 15.3
eh = 2.7
rbox(rx, ey0, 8, eh, '#F5F5F5', alpha=0.7, lw=1.5, ec='#DDD')
ax.text(rx + 4, ey0 + eh - 0.4, '5 Evidence Layers',
        fontsize=14, fontweight='bold', ha='center', color=C['dark'])

evs = [
    ('DA Significance',     '−log₁₀(padj) / 30',    C['ev'][0]),
    ('DA Effect Size',      '|log₂FC| / 5',          C['ev'][1]),
    ('Motif PWM Score',     'PWM score / 20',         C['ev'][2]),
    ('Peak-Gene Linkage',   '|z-score| / 10',         C['ev'][3]),
    ('TF-Gene Correlation', '|Spearman ρ| / 0.4',     C['ev'][4]),
]
for j, (name, formula, col) in enumerate(evs):
    ey = ey0 + eh - 1.0 - j * 0.3
    ax.plot(rx + 0.4, ey, 'o', color=col, markersize=8, zorder=3)
    ax.text(rx + 0.8, ey, name, fontsize=9, fontweight='bold',
            va='center', color=col)
    ax.text(rx + 4.5, ey, formula, fontsize=8.5, va='center',
            color='#888', family='monospace')

# composite formula - separate line below
ax.text(rx + 4, ey0 + 0.3,
        'Composite = ( ev₁ × ev₂ × ev₃ × ev₄ × ev₅ )^(1/5)',
        fontsize=8.5, ha='center', color='#666', family='monospace', style='italic')

# --- Section: Tier Classification ---
ty0 = 13.0
th = 1.9
rbox(rx, ty0, 8, th, '#F5F5F5', alpha=0.7, lw=1.5, ec='#DDD')
ax.text(rx + 4, ty0 + th - 0.4, 'Evidence Tier Classification',
        fontsize=14, fontweight='bold', ha='center', color=C['dark'])

tiers_data = [
    ('Tier 1', 'All 5 evidence layers > 0',  '419',    '1.8%',  C['step1']),
    ('Tier 2', 'PWM + Linkage + DA > 0',     '17,200', '73.0%', '#E8A838'),
    ('Tier 3', 'DA + Linkage > 0',           '5,952',  '25.2%', C['lightgrey']),
]
for j, (tier, criteria, count, pct, col) in enumerate(tiers_data):
    ty = ty0 + th - 0.95 - j * 0.42
    rbox(rx + 0.2, ty - 0.16, 1.0, 0.34, col, alpha=0.9)
    ax.text(rx + 0.7, ty, tier, fontsize=8.5, fontweight='bold',
            ha='center', va='center', color='white')
    ax.text(rx + 1.5, ty, criteria, fontsize=9, va='center', color='#555')
    ax.text(rx + 6.3, ty, count, fontsize=10, fontweight='bold',
            va='center', ha='right', color=col)
    ax.text(rx + 7.5, ty, f'({pct})', fontsize=8.5, va='center',
            ha='right', color='#AAA')

# --- Connecting dashed arrows (pipeline → evidence) ---
arr(7.5, 18.7, rx, ey0 + eh - 0.8, C['step1'], 1.2, '--')
arr(7.5, 16.3, rx, ey0 + eh - 1.55, C['step2'], 1.2, '--')
arr(7.5, 13.7, rx, ey0 + eh - 2.0, C['step3'], 1.2, '--')
arr(7.5, 11.1, rx, ey0 + 0.55, C['step4'], 1.2, '--')


# ═══════════════════════════════════════════════════════════════════════════
# RIGHT BOTTOM: Top Master Regulators Table (x: 8.8 ~ 16.8)
# ═══════════════════════════════════════════════════════════════════════════
tb_y0 = 7.8
tb_h = 4.6
rbox(rx, tb_y0, 8, tb_h, C['highlight'], alpha=0.04, lw=2, ec=C['highlight'])

ax.text(rx + 4, tb_y0 + tb_h - 0.45, 'Top 10 Tier 1 Master Regulators',
        fontsize=14, fontweight='bold', ha='center', color=C['highlight'])

# Table header
hdr_y = tb_y0 + tb_h - 1.0
ax.text(rx + 0.5, hdr_y, 'Rank', fontsize=9, fontweight='bold', color='#999', va='center')
ax.text(rx + 1.4, hdr_y, 'TF', fontsize=9, fontweight='bold', color='#999', va='center')
ax.text(rx + 4.2, hdr_y, 'Targets', fontsize=9, fontweight='bold', color='#999',
        va='center', ha='center')
ax.text(rx + 5.5, hdr_y, 'DARs', fontsize=9, fontweight='bold', color='#999',
        va='center', ha='center')
ax.text(rx + 6.8, hdr_y, 'Score', fontsize=9, fontweight='bold', color='#999',
        va='center', ha='center')

ax.plot([rx + 0.3, rx + 7.7], [hdr_y - 0.18, hdr_y - 0.18],
        color='#DDD', lw=0.8)

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
    ry = hdr_y - 0.55 - j * 0.34
    is_top = rank == '1'
    tf_color = C['highlight'] if is_top else C['dark']
    tf_weight = 'bold' if is_top else 'normal'
    num_color = C['highlight'] if is_top else '#666'

    if is_top:
        rbox(rx + 0.2, ry - 0.14, 7.5, 0.32, C['highlight'], alpha=0.08)

    ax.text(rx + 0.7, ry, rank, fontsize=9, ha='center', va='center', color='#BBB')
    ax.text(rx + 1.4, ry, tf, fontsize=10, fontweight=tf_weight,
            va='center', color=tf_color)
    ax.text(rx + 4.2, ry, tgt, fontsize=10, fontweight=tf_weight,
            ha='center', va='center', color=num_color)
    ax.text(rx + 5.5, ry, dars, fontsize=10, ha='center', va='center', color=num_color)
    ax.text(rx + 6.8, ry, score, fontsize=9, ha='center', va='center', color=num_color)


# ═══════════════════════════════════════════════════════════════════════════
# BOTTOM: HNF4A Key Finding
# ═══════════════════════════════════════════════════════════════════════════
hy0 = 0.5
hh = 6.8
rbox(0.5, hy0, 17, hh, C['highlight'], alpha=0.05, lw=2.5, ec=C['highlight'])

ax.text(9, hy0 + hh - 0.55, 'Key Finding',
        fontsize=18, fontweight='bold', ha='center', color=C['highlight'])
ax.plot([3, 15], [hy0 + hh - 0.95, hy0 + hh - 0.95], color=C['highlight'],
        lw=0.8, alpha=0.4)

# Left: HNF4A title
ax.text(1.5, hy0 + hh - 1.7, 'HNF4A', fontsize=36, fontweight='bold',
        color=C['highlight'], va='center')
ax.text(1.5, hy0 + hh - 2.4, 'Hepatocyte Nuclear Factor 4 Alpha',
        fontsize=11, color='#999', va='center', style='italic')
ax.text(1.5, hy0 + hh - 2.9, '#1 Master Regulator of chromatin remodeling',
        fontsize=10, color='#777', va='center')
ax.text(1.5, hy0 + hh - 3.3, 'in Cirrhosis-associated Hepatocytes',
        fontsize=10, color='#777', va='center')

# Stat boxes
stat_items = [
    ('55', 'target genes'),
    ('43', 'mediating DARs'),
    ('0.301', 'causal score'),
]
for k, (val, label) in enumerate(stat_items):
    sx = 1.3 + k * 3.0
    sy = hy0 + 0.5
    rbox(sx, sy, 2.5, 2.0, C['highlight'], alpha=0.08, lw=1.5, ec=C['highlight'])
    ax.text(sx + 1.25, sy + 1.3, val, fontsize=26, fontweight='bold',
            ha='center', va='center', color=C['highlight'])
    ax.text(sx + 1.25, sy + 0.45, label, fontsize=10,
            ha='center', va='center', color='#888')

# Right: Notable target genes
ntx = 10.5
ax.text(ntx, hy0 + hh - 1.7, 'Notable HNF4A Target Genes',
        fontsize=13, fontweight='bold', color=C['highlight'], va='center')

targets = [
    ('CYP3A4',  'Drug metabolism (cytochrome P450)'),
    ('AGT',     'Angiotensinogen — blood pressure regulation'),
    ('GPT',     'ALT enzyme — liver function marker'),
    ('SEC14L2', 'Lipid transport & cholesterol metabolism'),
    ('TF',      'Transferrin — iron transport'),
    ('POR',     'P450 oxidoreductase'),
    ('SLC27A5', 'Bile acid metabolism'),
]

for k, (gene, func) in enumerate(targets):
    ty = hy0 + hh - 2.3 - k * 0.5
    # gene badge
    rbox(ntx, ty - 0.15, 1.7, 0.35, C['highlight'], alpha=0.12)
    ax.text(ntx + 0.85, ty, gene, fontsize=10, fontweight='bold',
            ha='center', va='center', color=C['highlight'])
    ax.text(ntx + 2.0, ty, func, fontsize=9, va='center', color='#888')


# ═══════════════════════════════════════════════════════════════════════════
# SAVE
# ═══════════════════════════════════════════════════════════════════════════
save_dir = '/data1/project/yeonu/065_multi_atac/analysis/figures'
for fmt in ['png', 'pdf']:
    path = f'{save_dir}/pipeline_overview.{fmt}'
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f'Saved: {path}')

plt.close()
print('Done!')
