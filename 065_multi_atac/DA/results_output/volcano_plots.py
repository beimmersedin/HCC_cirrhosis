#!/usr/bin/env python
# coding: utf-8
"""
Volcano Plots for DA Results (TC: Cirrhosis vs Hepatitis)
Each cell type gets its own volcano plot.
"""

import os
import glob
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings('ignore')

# ── Settings ──
RESULT_DIR = '/data1/project/yeonu/065_multi_atac/DA/results_output'
SAVE_DIR = os.path.join(RESULT_DIR, 'volcano_figures')
os.makedirs(SAVE_DIR, exist_ok=True)

FC_THRESH = 0.5       # |log2FC| threshold
PADJ_THRESH = 0.05    # adjusted p-value threshold
TOP_N_LABEL = 10      # number of top peaks to label

plt.rcParams.update({
    'figure.dpi': 120,
    'font.size': 10,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
})

# ── Find all DA result CSVs ──
csv_files = sorted(glob.glob(os.path.join(RESULT_DIR, 'DA_results_*.csv')))
print(f'Found {len(csv_files)} DA result files.\n')

for csv_path in csv_files:
    # Parse cell type name from filename
    fname = os.path.basename(csv_path)
    cell_type = fname.replace('DA_results_', '').replace('.csv', '')

    # Load data
    df = pd.read_csv(csv_path, index_col=0)
    df = df.dropna(subset=['padj', 'log2FC'])

    # -log10(adjusted p-value), cap at 300 to avoid inf
    df['-log10_padj'] = -np.log10(df['padj'].clip(lower=1e-300))

    # Classify peaks
    df['category'] = 'NS'  # not significant
    df.loc[(df['padj'] < PADJ_THRESH) & (df['log2FC'] > FC_THRESH), 'category'] = 'Up_Cirrhosis'
    df.loc[(df['padj'] < PADJ_THRESH) & (df['log2FC'] < -FC_THRESH), 'category'] = 'Up_Hepatitis'
    df.loc[(df['padj'] < PADJ_THRESH) &
           (df['log2FC'].abs() <= FC_THRESH), 'category'] = 'Sig_lowFC'

    n_up = (df['category'] == 'Up_Cirrhosis').sum()
    n_down = (df['category'] == 'Up_Hepatitis').sum()
    n_sig_low = (df['category'] == 'Sig_lowFC').sum()
    n_total = len(df)

    # ── Plot ──
    fig, ax = plt.subplots(figsize=(8, 6))

    colors = {
        'NS': '#BBBBBB',
        'Sig_lowFC': '#FFC107',
        'Up_Cirrhosis': '#E53935',
        'Up_Hepatitis': '#1E88E5',
    }

    # Plot each category (NS first for z-order)
    for cat in ['NS', 'Sig_lowFC', 'Up_Hepatitis', 'Up_Cirrhosis']:
        sub = df[df['category'] == cat]
        if len(sub) == 0:
            continue
        ax.scatter(sub['log2FC'], sub['-log10_padj'],
                   s=6, alpha=0.5, color=colors[cat], edgecolors='none',
                   label=cat, rasterized=True)

    # Threshold lines
    ax.axhline(y=-np.log10(PADJ_THRESH), color='grey', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.axvline(x=FC_THRESH, color='grey', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.axvline(x=-FC_THRESH, color='grey', linestyle='--', linewidth=0.8, alpha=0.6)

    # Label top significant peaks
    sig_df = df[df['padj'] < PADJ_THRESH].nsmallest(TOP_N_LABEL, 'padj')
    for peak, row in sig_df.iterrows():
        # Shorten peak label: chr8-8227872-8228698 -> chr8:8.2M
        parts = peak.split('-')
        if len(parts) == 3:
            chrom = parts[0]
            start_mb = int(parts[1]) / 1e6
            label = f'{chrom}:{start_mb:.1f}M'
        else:
            label = peak[:20]
        ax.annotate(label, (row['log2FC'], row['-log10_padj']),
                    fontsize=6.5, alpha=0.85,
                    xytext=(5, 5), textcoords='offset points',
                    arrowprops=dict(arrowstyle='-', color='grey', lw=0.5))

    # Labels & title
    ax.set_xlabel('log2FC (Cirrhosis vs Hepatitis)')
    ax.set_ylabel('-log10(adjusted p-value)')
    ax.set_title(f'Volcano Plot — {cell_type}\n'
                 f'(n={n_total:,} peaks | Up Cirrhosis: {n_up} | Up Hepatitis: {n_down})',
                 fontweight='bold')

    # Legend
    # legend_elements = [
    #     Line2D([0], [0], marker='o', color='w', markerfacecolor=colors['Up_Cirrhosis'],
    #            markersize=7, label=f'Up in Cirrhosis ({n_up})'),
    #     Line2D([0], [0], marker='o', color='w', markerfacecolor=colors['Up_Hepatitis'],
    #            markersize=7, label=f'Up in Hepatitis ({n_down})'),
    #     Line2D([0], [0], marker='o', color='w', markerfacecolor=colors['Sig_lowFC'],
    #            markersize=7, label=f'Sig, |FC|≤{FC_THRESH} ({n_sig_low})'),
    #     Line2D([0], [0], marker='o', color='w', markerfacecolor=colors['NS'],
    #            markersize=7, label=f'NS ({n_total - n_up - n_down - n_sig_low})'),
    # ]
    # ax.legend(handles=legend_elements, fontsize=8, loc='upper right',
    #           framealpha=0.9, edgecolor='lightgrey')

    # plt.tight_layout()

    # Save
    out_path = os.path.join(SAVE_DIR, f'volcano_{cell_type}.png')
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'  [{cell_type}] Up={n_up}, Down={n_down} → {out_path}')

print(f'\nAll volcano plots saved to: {SAVE_DIR}')
