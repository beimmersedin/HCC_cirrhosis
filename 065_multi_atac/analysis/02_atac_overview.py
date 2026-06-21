#!/usr/bin/env python
# coding: utf-8
# # ATAC-seq Overview — HCC Multiome
# Comprehensive ATAC QC & feature exploration

import os
import matplotlib
matplotlib.use('Agg')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings('ignore')

# Style
sns.set_style('whitegrid')
plt.rcParams.update({
    'figure.dpi': 120,
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 11,
})

# Paths
BASE = '/data1/project2/yeonu/cirrhosis'
MULTI = f'{BASE}/multiome'
SAVE_DIR = '/data1/project/yeonu/065_multi_atac/analysis/figures_atac'
os.makedirs(SAVE_DIR, exist_ok=True)

_fig_counter = [0]
def savefig(name=None):
    _fig_counter[0] += 1
    fname = name or f'fig_{_fig_counter[0]:02d}'
    path = f'{SAVE_DIR}/{fname}.png'
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {path}')

# Colors
group_color_map = {'NT': '#4C72B0', 'TC': '#DD8452', 'PL': '#55A868'}
status_color = {'Cirrhosis': '#C44E52', 'Hepatitis': '#8172B3'}
group_order = ['NT', 'PL', 'TC']

print('Libraries loaded.')

# ============================================================
# Data Loading
# ============================================================
rna_meta = pd.read_csv(f'{MULTI}/rna_cell_metadata.csv', index_col=0)
umap_df  = pd.read_csv(f'{MULTI}/umap_celltype.csv', index_col=0)
peak_meta = pd.read_csv(f'{MULTI}/atac_peak_metadata.csv', index_col=0)

df = rna_meta.join(umap_df[['UMAP_1', 'UMAP_2']], how='left')
# celltype columns already exist in rna_meta; grab from umap_df only if missing
if 'celltype_assign' not in df.columns:
    df = df.join(umap_df[['celltype_assign', 'celltype_assign_fine']], how='left')

celltypes = sorted(df['celltype_assign'].dropna().unique())
n_ct = len(celltypes)
ct_cmap = plt.cm.get_cmap('tab10', max(n_ct, 10))
ct_colors = {ct: ct_cmap(i) for i, ct in enumerate(celltypes)}

print(f'Cells loaded: {len(df):,}')
print(f'Peaks loaded: {len(peak_meta):,}')

# ============================================================
# 1. ATAC QC Summary Table — Group × Status
# ============================================================
print('\n' + '=' * 70)
print('1. ATAC QC SUMMARY TABLE (Group × Status)')
print('=' * 70)

qc_cols = ['nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'nucleosome_signal']

# (A) By Group (NT → PL → TC)
qc_by_group = df.groupby('group')[qc_cols].agg(['median', 'mean', 'std']).round(2)
qc_by_group = qc_by_group.reindex(group_order)
print('\n[By Group]')
print(qc_by_group.to_string())

# (B) By Group × Status
qc_by_gs = df.groupby(['group', 'status'])[qc_cols].agg(['count', 'median']).round(2)
qc_by_gs = qc_by_gs.reindex(group_order, level=0)
print('\n[By Group × Status]')
print(qc_by_gs.to_string())

# (C) Cell count table
ct_table = pd.crosstab(df['status'], df['group'])
ct_table = ct_table.reindex(columns=group_order)
print('\n[Cell Count: Status × Group]')
print(ct_table.to_string())

# Summary figure
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
qc_labels = ['ATAC Fragment Count', 'Accessible Peaks', 'TSS Enrichment', 'Nucleosome Signal']

for ax, col, label in zip(axes.flat, qc_cols, qc_labels):
    data_plot = []
    for g in group_order:
        data_plot.append(df.loc[df['group'] == g, col].values)
    bp = ax.boxplot(data_plot, labels=group_order, patch_artist=True,
                    flierprops=dict(markersize=1, alpha=0.3),
                    medianprops=dict(color='black', linewidth=1.5))
    for patch, g in zip(bp['boxes'], group_order):
        patch.set_facecolor(group_color_map[g])
        patch.set_alpha(0.7)
    ax.set_title(label, fontsize=13, fontweight='bold')
    ax.set_xlabel('Group (NT → PL → TC)')
    # Add median text
    for i, g in enumerate(group_order):
        med = df.loc[df['group'] == g, col].median()
        ax.text(i + 1, med, f' {med:.1f}', va='bottom', ha='left', fontsize=9, color='red')

plt.suptitle('1. ATAC QC Summary by Group (NT → PL → TC)', fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
savefig('01_atac_qc_summary_by_group')

# ============================================================
# 2. Peak Count Distribution by Group
# ============================================================
print('\n' + '=' * 70)
print('2. PEAK COUNT DISTRIBUTION BY GROUP')
print('=' * 70)

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# (A) Histogram
for g in group_order:
    sub = df[df['group'] == g]
    axes[0].hist(sub['nFeature_ATAC'], bins=80, alpha=0.5,
                 color=group_color_map[g], label=f'{g} (n={len(sub):,})', density=True)
axes[0].set_xlabel('nFeature_ATAC (Accessible Peaks per Cell)')
axes[0].set_ylabel('Density')
axes[0].set_title('(A) Peak Count Distribution')
axes[0].legend(fontsize=9)

# (B) Violin
sns.violinplot(data=df, x='group', y='nFeature_ATAC', ax=axes[1],
               order=group_order, palette=group_color_map, inner='box', cut=0)
axes[1].set_title('(B) Peak Count Violin')
axes[1].set_xlabel('Group')

# (C) Per sample boxplot
sample_order = sorted(df['orig.ident'].unique())
sample_palette = {s: group_color_map.get(df.loc[df['orig.ident']==s, 'group'].iloc[0], '#999')
                  for s in sample_order}
sns.boxplot(data=df, x='orig.ident', y='nFeature_ATAC', ax=axes[2],
            order=sample_order, palette=sample_palette, fliersize=0.3)
axes[2].set_xticklabels(axes[2].get_xticklabels(), rotation=45, ha='right', fontsize=8)
axes[2].set_title('(C) Peak Count per Sample')
axes[2].set_xlabel('')

plt.suptitle('2. ATAC Peak Count Distribution', fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
savefig('02_peak_count_distribution_by_group')

# ============================================================
# 3. Peak Count Distribution by Cell Type
# ============================================================
print('\n' + '=' * 70)
print('3. PEAK COUNT DISTRIBUTION BY CELL TYPE')
print('=' * 70)

ct_median_peaks = df.groupby('celltype_assign')['nFeature_ATAC'].median().sort_values()
print('\n[Median Accessible Peaks per Cell Type]')
print(ct_median_peaks.to_string())

fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# (A) Boxplot by cell type
ct_order_peaks = ct_median_peaks.index.tolist()
sns.boxplot(data=df, y='celltype_assign', x='nFeature_ATAC', ax=axes[0],
            order=ct_order_peaks, palette=[ct_colors.get(ct, '#999') for ct in ct_order_peaks],
            fliersize=0.3, linewidth=0.7)
axes[0].set_xlabel('nFeature_ATAC (Accessible Peaks)')
axes[0].set_ylabel('')
axes[0].set_title('(A) Accessible Peaks by Cell Type')
for i, ct in enumerate(ct_order_peaks):
    med = ct_median_peaks[ct]
    axes[0].text(med + 50, i, f'{med:.0f}', va='center', fontsize=8, color='red')

# (B) Violin by cell type, split by group
# Simplified: stacked bar of median peaks
ct_group_median = df.groupby(['celltype_assign', 'group'])['nFeature_ATAC'].median().unstack()
ct_group_median = ct_group_median.reindex(columns=group_order)
ct_group_median = ct_group_median.loc[ct_order_peaks]
ct_group_median.plot.barh(ax=axes[1],
    color=[group_color_map[g] for g in group_order], width=0.7)
axes[1].set_xlabel('Median nFeature_ATAC')
axes[1].set_ylabel('')
axes[1].set_title('(B) Median Peaks by Cell Type & Group')
axes[1].legend(title='Group', fontsize=9)

plt.suptitle('3. ATAC Peak Count by Cell Type', fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
savefig('03_peak_count_by_celltype')

# ============================================================
# 4. TSS Enrichment by Cell Type
# ============================================================
print('\n' + '=' * 70)
print('4. TSS ENRICHMENT BY CELL TYPE')
print('=' * 70)

ct_median_tss = df.groupby('celltype_assign')['TSS.enrichment'].median().sort_values()
print('\n[Median TSS Enrichment per Cell Type]')
print(ct_median_tss.to_string())

fig, axes = plt.subplots(1, 3, figsize=(20, 6))

# (A) Boxplot by cell type
ct_order_tss = ct_median_tss.index.tolist()
sns.boxplot(data=df, y='celltype_assign', x='TSS.enrichment', ax=axes[0],
            order=ct_order_tss, palette=[ct_colors.get(ct, '#999') for ct in ct_order_tss],
            fliersize=0.3, linewidth=0.7)
axes[0].set_xlabel('TSS Enrichment')
axes[0].set_ylabel('')
axes[0].set_title('(A) TSS Enrichment by Cell Type')
axes[0].axvline(x=2, color='red', linestyle='--', alpha=0.5, label='QC threshold (2)')
axes[0].legend(fontsize=8)

# (B) By group
sns.boxplot(data=df, x='group', y='TSS.enrichment', ax=axes[1],
            order=group_order, palette=group_color_map, fliersize=0.3)
axes[1].set_title('(B) TSS Enrichment by Group')
axes[1].axhline(y=2, color='red', linestyle='--', alpha=0.5)

# (C) By sample
sns.boxplot(data=df, x='orig.ident', y='TSS.enrichment', ax=axes[2],
            order=sample_order, palette=sample_palette, fliersize=0.3)
axes[2].set_xticklabels(axes[2].get_xticklabels(), rotation=45, ha='right', fontsize=8)
axes[2].set_title('(C) TSS Enrichment by Sample')
axes[2].axhline(y=2, color='red', linestyle='--', alpha=0.5)
axes[2].set_xlabel('')

plt.suptitle('4. TSS Enrichment (higher = better ATAC quality)', fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
savefig('04_tss_enrichment_by_celltype')

# ============================================================
# 5. Nucleosome Signal by Cell Type
# ============================================================
print('\n' + '=' * 70)
print('5. NUCLEOSOME SIGNAL BY CELL TYPE')
print('=' * 70)

ct_median_nuc = df.groupby('celltype_assign')['nucleosome_signal'].median().sort_values()
print('\n[Median Nucleosome Signal per Cell Type]')
print(ct_median_nuc.to_string())

fig, axes = plt.subplots(1, 3, figsize=(20, 6))

# (A) Boxplot by cell type
ct_order_nuc = ct_median_nuc.index.tolist()
sns.boxplot(data=df, y='celltype_assign', x='nucleosome_signal', ax=axes[0],
            order=ct_order_nuc, palette=[ct_colors.get(ct, '#999') for ct in ct_order_nuc],
            fliersize=0.3, linewidth=0.7)
axes[0].set_xlabel('Nucleosome Signal')
axes[0].set_ylabel('')
axes[0].set_title('(A) Nucleosome Signal by Cell Type')
axes[0].axvline(x=2, color='red', linestyle='--', alpha=0.5, label='QC threshold (2)')
axes[0].legend(fontsize=8)

# (B) By group
sns.boxplot(data=df, x='group', y='nucleosome_signal', ax=axes[1],
            order=group_order, palette=group_color_map, fliersize=0.3)
axes[1].set_title('(B) Nucleosome Signal by Group')
axes[1].axhline(y=2, color='red', linestyle='--', alpha=0.5)

# (C) By sample
sns.boxplot(data=df, x='orig.ident', y='nucleosome_signal', ax=axes[2],
            order=sample_order, palette=sample_palette, fliersize=0.3)
axes[2].set_xticklabels(axes[2].get_xticklabels(), rotation=45, ha='right', fontsize=8)
axes[2].set_title('(C) Nucleosome Signal by Sample')
axes[2].axhline(y=2, color='red', linestyle='--', alpha=0.5)
axes[2].set_xlabel('')

plt.suptitle('5. Nucleosome Signal (lower = better ATAC quality)', fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
savefig('05_nucleosome_signal_by_celltype')

# ============================================================
# 6. ATAC Fragment Distribution by Group
# ============================================================
print('\n' + '=' * 70)
print('6. ATAC FRAGMENT DISTRIBUTION BY GROUP')
print('=' * 70)

fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# Row 1: nCount_ATAC
# (A) Histogram
for g in group_order:
    sub = df[df['group'] == g]
    axes[0, 0].hist(sub['nCount_ATAC'], bins=80, alpha=0.5,
                    color=group_color_map[g], label=g, density=True)
axes[0, 0].set_xlabel('nCount_ATAC')
axes[0, 0].set_ylabel('Density')
axes[0, 0].set_title('(A) Fragment Count Distribution')
axes[0, 0].legend()

# (B) Log-scale histogram
for g in group_order:
    sub = df[df['group'] == g]
    axes[0, 1].hist(np.log10(sub['nCount_ATAC'] + 1), bins=80, alpha=0.5,
                    color=group_color_map[g], label=g, density=True)
axes[0, 1].set_xlabel('log10(nCount_ATAC + 1)')
axes[0, 1].set_ylabel('Density')
axes[0, 1].set_title('(B) Fragment Count (log10)')
axes[0, 1].legend()

# (C) ECDF
for g in group_order:
    sub = df[df['group'] == g]
    sorted_vals = np.sort(sub['nCount_ATAC'].values)
    cdf = np.arange(1, len(sorted_vals) + 1) / len(sorted_vals)
    axes[0, 2].plot(sorted_vals, cdf, color=group_color_map[g], label=g, linewidth=1.5)
axes[0, 2].set_xlabel('nCount_ATAC')
axes[0, 2].set_ylabel('Cumulative Fraction')
axes[0, 2].set_title('(C) ECDF of Fragment Count')
axes[0, 2].legend()

# Row 2: nFeature_ATAC
# (D) Histogram
for g in group_order:
    sub = df[df['group'] == g]
    axes[1, 0].hist(sub['nFeature_ATAC'], bins=80, alpha=0.5,
                    color=group_color_map[g], label=g, density=True)
axes[1, 0].set_xlabel('nFeature_ATAC')
axes[1, 0].set_ylabel('Density')
axes[1, 0].set_title('(D) Peak Count Distribution')
axes[1, 0].legend()

# (E) Log-scale
for g in group_order:
    sub = df[df['group'] == g]
    axes[1, 1].hist(np.log10(sub['nFeature_ATAC'] + 1), bins=80, alpha=0.5,
                    color=group_color_map[g], label=g, density=True)
axes[1, 1].set_xlabel('log10(nFeature_ATAC + 1)')
axes[1, 1].set_ylabel('Density')
axes[1, 1].set_title('(E) Peak Count (log10)')
axes[1, 1].legend()

# (F) nCount vs nFeature scatter
for g in group_order:
    sub = df[df['group'] == g]
    axes[1, 2].scatter(sub['nCount_ATAC'], sub['nFeature_ATAC'],
                       s=0.5, alpha=0.15, color=group_color_map[g], label=g)
axes[1, 2].set_xlabel('nCount_ATAC (Fragments)')
axes[1, 2].set_ylabel('nFeature_ATAC (Peaks)')
axes[1, 2].set_title('(F) Fragments vs Peaks')
axes[1, 2].legend(markerscale=5)

plt.suptitle('6. ATAC Fragment Distribution by Group', fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
savefig('06_atac_fragment_distribution')

# ============================================================
# 7. RNA vs ATAC Correlation by Cell Type
# ============================================================
print('\n' + '=' * 70)
print('7. RNA vs ATAC CORRELATION BY CELL TYPE')
print('=' * 70)

fig, axes = plt.subplots(2, 3, figsize=(18, 11))

# (A) nCount_RNA vs nCount_ATAC colored by cell type
for ct in celltypes:
    sub = df[df['celltype_assign'] == ct]
    axes[0, 0].scatter(sub['nCount_RNA'], sub['nCount_ATAC'],
                       s=0.5, alpha=0.15, color=ct_colors[ct], label=ct)
axes[0, 0].set_xlabel('nCount_RNA')
axes[0, 0].set_ylabel('nCount_ATAC')
axes[0, 0].set_title('(A) RNA vs ATAC Counts (by Cell Type)')
axes[0, 0].legend(markerscale=5, fontsize=7, loc='upper right')

# (B) nFeature_RNA vs nFeature_ATAC colored by cell type
for ct in celltypes:
    sub = df[df['celltype_assign'] == ct]
    axes[0, 1].scatter(sub['nFeature_RNA'], sub['nFeature_ATAC'],
                       s=0.5, alpha=0.15, color=ct_colors[ct], label=ct)
axes[0, 1].set_xlabel('nFeature_RNA (Genes)')
axes[0, 1].set_ylabel('nFeature_ATAC (Peaks)')
axes[0, 1].set_title('(B) RNA vs ATAC Features (by Cell Type)')
axes[0, 1].legend(markerscale=5, fontsize=7, loc='upper right')

# (C) RNA/ATAC ratio by group (NT → PL → TC)
df['rna_atac_ratio'] = df['nCount_RNA'] / (df['nCount_ATAC'] + 1)
sns.boxplot(data=df, x='group', y='rna_atac_ratio', ax=axes[0, 2],
            order=group_order, palette=group_color_map, fliersize=0.3)
axes[0, 2].set_ylim(0, df['rna_atac_ratio'].quantile(0.95))
axes[0, 2].set_title('(C) RNA/ATAC Ratio by Group')
axes[0, 2].set_xlabel('Group')

# (D) RNA/ATAC ratio by cell type
ct_order_ratio = df.groupby('celltype_assign')['rna_atac_ratio'].median().sort_values().index
sns.boxplot(data=df, y='celltype_assign', x='rna_atac_ratio', ax=axes[1, 0],
            order=ct_order_ratio, palette=[ct_colors.get(ct, '#999') for ct in ct_order_ratio],
            fliersize=0.3, linewidth=0.7)
axes[1, 0].set_xlim(0, df['rna_atac_ratio'].quantile(0.95))
axes[1, 0].set_xlabel('RNA / ATAC Count Ratio')
axes[1, 0].set_ylabel('')
axes[1, 0].set_title('(D) RNA/ATAC Ratio by Cell Type')

# (E) Correlation: median RNA vs median ATAC per cell type
ct_medians = df.groupby('celltype_assign').agg(
    med_rna=('nCount_RNA', 'median'),
    med_atac=('nCount_ATAC', 'median'),
    n_cells=('nCount_RNA', 'size')
)
for ct in ct_medians.index:
    axes[1, 1].scatter(ct_medians.loc[ct, 'med_rna'], ct_medians.loc[ct, 'med_atac'],
                       s=ct_medians.loc[ct, 'n_cells'] / 30, color=ct_colors[ct],
                       edgecolors='black', linewidth=0.5)
    axes[1, 1].annotate(ct, (ct_medians.loc[ct, 'med_rna'], ct_medians.loc[ct, 'med_atac']),
                        fontsize=8, ha='left', va='bottom')
axes[1, 1].set_xlabel('Median nCount_RNA')
axes[1, 1].set_ylabel('Median nCount_ATAC')
axes[1, 1].set_title('(E) Cell Type: Median RNA vs ATAC\n(size = n_cells)')

# (F) Heatmap: correlation between QC metrics per cell type
corr_cols = ['nCount_RNA', 'nFeature_RNA', 'nCount_ATAC', 'nFeature_ATAC',
             'TSS.enrichment', 'nucleosome_signal', 'percent.mt']
available_corr = [c for c in corr_cols if c in df.columns]
corr_mat = df[available_corr].corr()
sns.heatmap(corr_mat, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
            ax=axes[1, 2], square=True, linewidths=0.5,
            cbar_kws={'shrink': 0.8}, annot_kws={'fontsize': 8})
axes[1, 2].set_title('(F) QC Metric Correlations')
axes[1, 2].tick_params(axis='both', labelsize=8)

plt.suptitle('7. RNA vs ATAC Correlation', fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
savefig('07_rna_vs_atac_correlation')

# ============================================================
# 8. Peak Chromosome Distribution
# ============================================================
print('\n' + '=' * 70)
print('8. PEAK CHROMOSOME DISTRIBUTION')
print('=' * 70)

# Parse peak coordinates
peaks = peak_meta.index.to_series()
peak_info = peaks.str.extract(r'^(chr[\dXYM]+)-(\d+)-(\d+)$')
peak_info.columns = ['chr', 'start', 'end']
peak_info['start'] = peak_info['start'].astype(float)
peak_info['end'] = peak_info['end'].astype(float)
peak_info['width'] = peak_info['end'] - peak_info['start']

# Chromosome order
chr_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
chr_order = [c for c in chr_order if c in peak_info['chr'].values]

chr_counts = peak_info['chr'].value_counts().reindex(chr_order).fillna(0).astype(int)
chr_widths = peak_info.groupby('chr')['width'].median().reindex(chr_order)

print('\n[Peak Count per Chromosome]')
print(chr_counts.to_string())
print(f'\nTotal peaks parsed: {len(peak_info):,}')
print(f'Median peak width: {peak_info["width"].median():.0f} bp')

fig, axes = plt.subplots(2, 2, figsize=(16, 10))

# (A) Peak count per chromosome
axes[0, 0].bar(range(len(chr_order)), chr_counts.values, color='steelblue', edgecolor='white')
axes[0, 0].set_xticks(range(len(chr_order)))
axes[0, 0].set_xticklabels(chr_order, rotation=45, ha='right', fontsize=8)
axes[0, 0].set_ylabel('Number of Peaks')
axes[0, 0].set_title('(A) Peak Count per Chromosome')

# (B) Peak width distribution
axes[0, 1].hist(peak_info['width'].dropna(), bins=100, color='coral', edgecolor='white')
axes[0, 1].set_xlabel('Peak Width (bp)')
axes[0, 1].set_ylabel('Number of Peaks')
axes[0, 1].set_title(f'(B) Peak Width Distribution (median={peak_info["width"].median():.0f} bp)')
axes[0, 1].axvline(x=peak_info['width'].median(), color='red', linestyle='--', alpha=0.7)

# (C) Median peak width per chromosome
axes[1, 0].bar(range(len(chr_order)), chr_widths.values, color='mediumpurple', edgecolor='white')
axes[1, 0].set_xticks(range(len(chr_order)))
axes[1, 0].set_xticklabels(chr_order, rotation=45, ha='right', fontsize=8)
axes[1, 0].set_ylabel('Median Peak Width (bp)')
axes[1, 0].set_title('(C) Median Peak Width per Chromosome')

# (D) Cumulative peak coverage
chr_total_bp = peak_info.groupby('chr')['width'].sum().reindex(chr_order).fillna(0) / 1e6
axes[1, 1].bar(range(len(chr_order)), chr_total_bp.values, color='seagreen', edgecolor='white')
axes[1, 1].set_xticks(range(len(chr_order)))
axes[1, 1].set_xticklabels(chr_order, rotation=45, ha='right', fontsize=8)
axes[1, 1].set_ylabel('Total Peak Coverage (Mb)')
axes[1, 1].set_title('(D) Total Peak Coverage per Chromosome')

plt.suptitle('8. Peak Chromosome Distribution', fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
savefig('08_peak_chromosome_distribution')

# ============================================================
# BONUS: ATAC UMAP colored by QC metrics
# ============================================================
print('\n' + '=' * 70)
print('BONUS: UMAP colored by ATAC QC metrics')
print('=' * 70)

umap_valid = df.dropna(subset=['UMAP_1', 'UMAP_2'])

fig, axes = plt.subplots(2, 3, figsize=(20, 12))
umap_metrics = [
    ('nCount_ATAC', 'Fragments per Cell', 'YlOrRd'),
    ('nFeature_ATAC', 'Accessible Peaks per Cell', 'YlOrRd'),
    ('TSS.enrichment', 'TSS Enrichment', 'RdYlGn'),
    ('nucleosome_signal', 'Nucleosome Signal', 'RdYlGn_r'),
    ('nCount_RNA', 'RNA UMI Count', 'YlGnBu'),
    ('percent.mt', 'Mitochondrial %', 'RdYlBu_r'),
]

for ax, (col, label, cmap) in zip(axes.flat, umap_metrics):
    vmin = umap_valid[col].quantile(0.02)
    vmax = umap_valid[col].quantile(0.98)
    sc = ax.scatter(umap_valid['UMAP_1'], umap_valid['UMAP_2'],
                    c=umap_valid[col], cmap=cmap, s=0.5, alpha=0.4,
                    vmin=vmin, vmax=vmax)
    plt.colorbar(sc, ax=ax, shrink=0.7, label=col)
    ax.set_title(label, fontsize=12, fontweight='bold')
    ax.set_xlabel('UMAP_1')
    ax.set_ylabel('UMAP_2')

plt.suptitle('UMAP colored by ATAC & RNA QC Metrics', fontsize=15, fontweight='bold', y=1.01)
plt.tight_layout()
savefig('09_umap_atac_qc_metrics')

# ============================================================
# BONUS 2: ATAC QC by Cell Type — Heatmap Summary
# ============================================================
print('\n' + '=' * 70)
print('BONUS 2: ATAC QC Heatmap by Cell Type')
print('=' * 70)

heatmap_data = df.groupby('celltype_assign')[qc_cols].median()
heatmap_data.columns = ['Fragment Count', 'Accessible Peaks', 'TSS Enrichment', 'Nucleosome Signal']

fig, ax = plt.subplots(figsize=(10, 6))
# Z-score normalization for visualization
heatmap_z = (heatmap_data - heatmap_data.mean()) / heatmap_data.std()
sns.heatmap(heatmap_z, annot=heatmap_data.round(1).values, fmt='', cmap='RdBu_r',
            center=0, linewidths=0.5, ax=ax,
            cbar_kws={'label': 'Z-score (color) / Raw median (text)'})
ax.set_title('ATAC QC Metrics by Cell Type\n(color=z-score, text=raw median)',
             fontsize=14, fontweight='bold')
ax.set_ylabel('')

plt.tight_layout()
savefig('10_atac_qc_heatmap_by_celltype')

# ============================================================
# Done
# ============================================================
print('\n' + '=' * 70)
print(f'All figures saved to: {SAVE_DIR}')
print(f'Total figures: {_fig_counter[0]}')
print('=' * 70)
