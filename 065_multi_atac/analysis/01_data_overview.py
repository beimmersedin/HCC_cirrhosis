#!/usr/bin/env python
# coding: utf-8

# # HCC Multi-omics Data Overview
#
# **Project**: Hepatocellular Carcinoma (HCC) Multi-omics Analysis
# **Data Platforms**: 10x Genomics Multiome (scRNA-seq + scATAC-seq) & NanoString GeoMx (Spatial Transcriptomics)
#
# ---
#
# ## Contents
# 1. **Data Basics** — 플랫폼, 모달리티, 샘플 구성
# 2. **Quality Control (QC)** — QC 지표, Doublet, Batch Effect
# 3. **Cell Annotation & Visualization** — UMAP, Marker, Cell Type 분포
# 4. **Multi-omics Integration** — RNA-ATAC 연결
# 5. **GeoMx Spatial Overview** — Spatial Transcriptomics 요약
# 6. **Future Plan** — 향후 분석 계획

# In[1]:

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
GEOMX = f'{BASE}/geomx'

SAVE_DIR = '/data1/project/yeonu/065_multi_atac/analysis/figures'
os.makedirs(SAVE_DIR, exist_ok=True)

_fig_counter = [0]
def savefig(name=None):
    _fig_counter[0] += 1
    fname = name or f'fig_{_fig_counter[0]:02d}'
    path = f'{SAVE_DIR}/{fname}.png'
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {path}')

print('Libraries loaded.')


# In[2]:


# ── Load Multiome metadata ──
rna_meta = pd.read_csv(f'{MULTI}/rna_cell_metadata.csv', index_col=0)
umap_df  = pd.read_csv(f'{MULTI}/umap_celltype.csv', index_col=0)

# Merge UMAP coords into rna_meta
df = rna_meta.join(umap_df[['UMAP_1', 'UMAP_2']], how='left')

# ── Load GeoMx metadata ──
geo_meta = pd.read_csv(f'{GEOMX}/segment_meta_merged_with_groups.csv')
geo_expr = pd.read_csv(f'{GEOMX}/expression_matrix_q3.csv')

print(f'Multiome cells : {len(df):,}')
print(f'GeoMx segments : {len(geo_meta):,}')
print(f'GeoMx genes    : {len(geo_expr):,}')


# ---
# ## 1. Data Basics — 데이터 기본 구성
# 
# ### 1.1 Platform & Modalities
# 
# | Platform | Technology | Modalities | Resolution |
# |----------|-----------|------------|------------|
# | **10x Multiome** | scRNA-seq + scATAC-seq | Transcriptome + Chromatin Accessibility | Single-cell |
# | **NanoString GeoMx** | Whole Transcriptome Atlas (WTA) | Spatial Transcriptomics | ROI-level (CD31+ / CD31- segments) |

# In[3]:


# ── 1.2 Sample Summary (Multiome) ──
print('=' * 60)
print('MULTIOME — Sample Summary')
print('=' * 60)

sample_summary = df.groupby(['orig.ident', 'group', 'status']).size().reset_index(name='n_cells')
sample_summary = sample_summary.sort_values('orig.ident')
print(f"\nTotal cells: {len(df):,}")
print(f"Total samples: {df['orig.ident'].nunique()}")
print(f"Groups: {sorted(df['group'].unique())}")
print(f"Status: {sorted(df['status'].unique())}")
print()
print(sample_summary.to_string())


# In[4]:


# ── 1.3 Cell counts by group & status ──
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

# (A) Cells per sample
sample_counts = df['orig.ident'].value_counts().sort_index()
colors_sample = []
group_color_map = {'NT': '#4C72B0', 'TC': '#DD8452', 'PL': '#55A868'}
for s in sample_counts.index:
    g = df.loc[df['orig.ident'] == s, 'group'].iloc[0]
    colors_sample.append(group_color_map.get(g, '#999999'))

axes[0].bar(range(len(sample_counts)), sample_counts.values, color=colors_sample)
axes[0].set_xticks(range(len(sample_counts)))
axes[0].set_xticklabels(sample_counts.index, rotation=45, ha='right', fontsize=8)
axes[0].set_ylabel('Number of Cells')
axes[0].set_title('(A) Cells per Sample')
legend_elements = [Line2D([0],[0], marker='s', color='w', markerfacecolor=c, markersize=10, label=g)
                   for g, c in group_color_map.items()]
axes[0].legend(handles=legend_elements, title='Group', fontsize=9)

# (B) Cells per group
group_counts = df['group'].value_counts()
axes[1].bar(group_counts.index, group_counts.values,
            color=[group_color_map.get(g, '#999') for g in group_counts.index])
axes[1].set_ylabel('Number of Cells')
axes[1].set_title('(B) Cells per Group')
for i, (g, v) in enumerate(zip(group_counts.index, group_counts.values)):
    axes[1].text(i, v + 100, f'{v:,}', ha='center', fontsize=9)

# (C) Cells per status
status_color = {'Cirrhosis': '#C44E52', 'Hepatitis': '#8172B3'}
status_counts = df['status'].value_counts()
axes[2].bar(status_counts.index, status_counts.values,
            color=[status_color.get(s, '#999') for s in status_counts.index])
axes[2].set_ylabel('Number of Cells')
axes[2].set_title('(C) Cells per Status')
for i, (s, v) in enumerate(zip(status_counts.index, status_counts.values)):
    axes[2].text(i, v + 100, f'{v:,}', ha='center', fontsize=9)

plt.tight_layout()
savefig('01_cell_counts_by_group_status')


# ---
# ## 2. Quality Control (QC)
# 
# ### 2.1 RNA QC Metrics
# - **nCount_RNA**: Total UMI counts per cell
# - **nFeature_RNA**: Number of detected genes per cell
# - **percent.mt**: Mitochondrial gene percentage (high = dying cells)

# In[5]:


# ── 2.1 RNA QC Violin Plots ──
qc_cols_rna = ['nCount_RNA', 'nFeature_RNA', 'percent.mt']
qc_labels = ['Total UMI Counts', 'Detected Genes', 'Mitochondrial %']

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

for ax, col, label in zip(axes, qc_cols_rna, qc_labels):
    sns.violinplot(data=df, x='group', y=col, ax=ax,
                   palette=group_color_map, inner='box', cut=0, scale='width')
    ax.set_title(label)
    ax.set_xlabel('Group')
    ax.set_ylabel(col)

plt.suptitle('RNA QC Metrics by Group', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
savefig('02_rna_qc_violin_by_group')


# In[6]:


# ── 2.1b RNA QC per Sample ──
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

sample_order = sorted(df['orig.ident'].unique())
sample_palette = {s: group_color_map.get(df.loc[df['orig.ident']==s, 'group'].iloc[0], '#999')
                  for s in sample_order}

for ax, col, label in zip(axes, qc_cols_rna, qc_labels):
    sns.boxplot(data=df, x='orig.ident', y=col, ax=ax,
                order=sample_order, palette=sample_palette,
                fliersize=0.5, linewidth=0.7)
    ax.set_title(label)
    ax.set_xlabel('')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=7)

plt.suptitle('RNA QC Metrics by Sample', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
savefig('03_rna_qc_boxplot_by_sample')


# In[7]:


# ── 2.1c RNA QC Scatter: nCount vs nFeature, colored by percent.mt ──
fig, ax = plt.subplots(figsize=(8, 6))
sc = ax.scatter(df['nCount_RNA'], df['nFeature_RNA'],
                c=df['percent.mt'], cmap='RdYlBu_r', s=1, alpha=0.3,
                vmin=0, vmax=df['percent.mt'].quantile(0.95))
plt.colorbar(sc, ax=ax, label='Mitochondrial %')
ax.set_xlabel('nCount_RNA (Total UMIs)')
ax.set_ylabel('nFeature_RNA (Detected Genes)')
ax.set_title('RNA QC: UMI Counts vs Gene Counts')
plt.tight_layout()
savefig('04_rna_qc_scatter')


# ### 2.2 ATAC QC Metrics
# - **nCount_ATAC**: Total ATAC fragments per cell
# - **nFeature_ATAC**: Number of accessible peaks per cell
# - **TSS.enrichment**: Transcription Start Site enrichment (higher = better quality)
# - **nucleosome_signal**: Nucleosome banding pattern (lower = better)

# In[8]:


# ── 2.2 ATAC QC Violin Plots ──
qc_cols_atac = ['nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'nucleosome_signal']
qc_labels_atac = ['ATAC Fragment Count', 'Accessible Peaks', 'TSS Enrichment', 'Nucleosome Signal']

fig, axes = plt.subplots(1, 4, figsize=(20, 5))

for ax, col, label in zip(axes, qc_cols_atac, qc_labels_atac):
    sns.violinplot(data=df, x='group', y=col, ax=ax,
                   palette=group_color_map, inner='box', cut=0, scale='width')
    ax.set_title(label)
    ax.set_xlabel('Group')
    ax.set_ylabel(col)

plt.suptitle('ATAC QC Metrics by Group', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
savefig('05_atac_qc_violin_by_group')


# In[9]:


# ── 2.2b ATAC QC Scatter: TSS enrichment vs nucleosome signal ──
fig, ax = plt.subplots(figsize=(8, 6))
for g in sorted(df['group'].unique()):
    sub = df[df['group'] == g]
    ax.scatter(sub['nucleosome_signal'], sub['TSS.enrichment'],
              s=1, alpha=0.2, label=g, color=group_color_map.get(g, '#999'))
ax.set_xlabel('Nucleosome Signal (lower = better)')
ax.set_ylabel('TSS Enrichment (higher = better)')
ax.set_title('ATAC QC: TSS Enrichment vs Nucleosome Signal')
ax.legend(markerscale=5)
plt.tight_layout()
savefig('06_atac_qc_scatter')


# ### 2.3 Doublet Detection
# 
# Doublet detection was performed using **scDblFinder**.  
# Cells classified as `singlet` were retained for downstream analysis.

# In[10]:


# ── 2.3 Doublet Detection Summary ──
if 'scDblFinder.class' in df.columns:
    dbl_counts = df['scDblFinder.class'].value_counts()
    print('Doublet Detection (scDblFinder):')
    print(dbl_counts)
    print(f"\nDoublet rate: {dbl_counts.get('doublet', 0) / len(df) * 100:.2f}%")
    print(f"All cells in this dataset are: {df['scDblFinder.class'].unique()}")
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    # (A) Doublet score distribution
    axes[0].hist(df['scDblFinder.score'], bins=50, color='steelblue', edgecolor='white')
    axes[0].set_xlabel('scDblFinder Score')
    axes[0].set_ylabel('Number of Cells')
    axes[0].set_title('(A) Doublet Score Distribution')
    
    # (B) Doublet class counts
    dbl_counts.plot.bar(ax=axes[1], color=['#55A868', '#C44E52'])
    axes[1].set_title('(B) Doublet Classification')
    axes[1].set_ylabel('Number of Cells')
    axes[1].tick_params(axis='x', rotation=0)
    for i, v in enumerate(dbl_counts.values):
        axes[1].text(i, v + 100, f'{v:,}', ha='center', fontsize=10)
    
    plt.tight_layout()
    savefig('07_doublet_detection')
else:
    print('scDblFinder column not found.')


# ### 2.4 Batch Effect Assessment
# 
# Multiple samples from different patients were integrated.  
# Below: UMAP colored by sample to assess batch mixing.

# In[11]:


# ── 2.4 Batch Effect: UMAP by Sample ──
umap_valid = df.dropna(subset=['UMAP_1', 'UMAP_2'])

samples = sorted(umap_valid['orig.ident'].unique())
n_samples = len(samples)
cmap_samples = plt.cm.get_cmap('tab20', n_samples)
sample_cmap = {s: cmap_samples(i) for i, s in enumerate(samples)}

fig, axes = plt.subplots(1, 2, figsize=(18, 7))

# (A) Colored by sample
for s in samples:
    sub = umap_valid[umap_valid['orig.ident'] == s]
    axes[0].scatter(sub['UMAP_1'], sub['UMAP_2'], s=0.5, alpha=0.3,
                    color=sample_cmap[s], label=s)
axes[0].set_title('UMAP — Colored by Sample (Batch Check)', fontsize=13)
axes[0].set_xlabel('UMAP_1')
axes[0].set_ylabel('UMAP_2')
axes[0].legend(markerscale=8, fontsize=7, loc='upper right', ncol=2,
               bbox_to_anchor=(1.25, 1.0))

# (B) Colored by group
for g in sorted(umap_valid['group'].unique()):
    sub = umap_valid[umap_valid['group'] == g]
    axes[1].scatter(sub['UMAP_1'], sub['UMAP_2'], s=0.5, alpha=0.3,
                    color=group_color_map.get(g, '#999'), label=g)
axes[1].set_title('UMAP — Colored by Group', fontsize=13)
axes[1].set_xlabel('UMAP_1')
axes[1].set_ylabel('UMAP_2')
axes[1].legend(markerscale=8, fontsize=10)

plt.tight_layout()
savefig('08_umap_batch_effect')


# In[12]:


# ── 2.4b QC Summary Table ──
qc_summary = df.groupby('orig.ident').agg(
    n_cells=('nCount_RNA', 'size'),
    median_nCount_RNA=('nCount_RNA', 'median'),
    median_nFeature_RNA=('nFeature_RNA', 'median'),
    median_percent_mt=('percent.mt', 'median'),
    median_nCount_ATAC=('nCount_ATAC', 'median'),
    median_nFeature_ATAC=('nFeature_ATAC', 'median'),
    median_TSS_enrichment=('TSS.enrichment', 'median'),
    median_nucleosome=('nucleosome_signal', 'median'),
).round(2)

print('QC Summary per Sample (median values):')
print(qc_summary.to_string())


# ---
# ## 3. Cell Annotation & Visualization
# 
# ### 3.1 UMAP — Major Cell Types
# 
# Cell type annotation was performed based on known marker genes.  
# - **Major types (celltype_assign)**: 9 categories  
# - **Fine types (celltype_assign_fine)**: ~20 subtypes

# In[13]:


# ── 3.1 UMAP — Major Cell Types ──
umap_valid = df.dropna(subset=['UMAP_1', 'UMAP_2'])

celltypes = sorted(umap_valid['celltype_assign'].dropna().unique())
n_ct = len(celltypes)
ct_cmap = plt.cm.get_cmap('tab10', max(n_ct, 10))
ct_colors = {ct: ct_cmap(i) for i, ct in enumerate(celltypes)}

fig, ax = plt.subplots(figsize=(10, 8))

# Shuffle for better visibility
plot_df = umap_valid.sample(frac=1, random_state=42)

for ct in celltypes:
    sub = plot_df[plot_df['celltype_assign'] == ct]
    ax.scatter(sub['UMAP_1'], sub['UMAP_2'], s=1, alpha=0.4,
              color=ct_colors[ct], label=f'{ct} ({len(sub):,})')

ax.set_title('UMAP — Major Cell Types', fontsize=14, fontweight='bold')
ax.set_xlabel('UMAP_1')
ax.set_ylabel('UMAP_2')
ax.legend(markerscale=8, fontsize=9, loc='upper right',
          bbox_to_anchor=(1.35, 1.0), title='Cell Type (n)')

plt.tight_layout()
savefig('09_umap_major_celltypes')


# In[14]:


# ── 3.1b UMAP — Fine Cell Types ──
fine_types = sorted(umap_valid['celltype_assign_fine'].dropna().unique())
n_fine = len(fine_types)
fine_cmap = plt.cm.get_cmap('tab20', max(n_fine, 20))
fine_colors = {ft: fine_cmap(i) for i, ft in enumerate(fine_types)}

fig, ax = plt.subplots(figsize=(10, 8))

for ft in fine_types:
    sub = plot_df[plot_df['celltype_assign_fine'] == ft]
    ax.scatter(sub['UMAP_1'], sub['UMAP_2'], s=1, alpha=0.4,
              color=fine_colors[ft], label=f'{ft} ({len(sub):,})')

ax.set_title('UMAP — Fine Cell Types', fontsize=14, fontweight='bold')
ax.set_xlabel('UMAP_1')
ax.set_ylabel('UMAP_2')
ax.legend(markerscale=8, fontsize=7, loc='upper right',
          bbox_to_anchor=(1.4, 1.0), title='Fine Cell Type (n)', ncol=1)

plt.tight_layout()
savefig('10_umap_fine_celltypes')


# In[15]:


# ── 3.1c UMAP — Split by Group ──
groups = sorted(umap_valid['group'].unique())
fig, axes = plt.subplots(1, len(groups), figsize=(6 * len(groups), 6))
if len(groups) == 1:
    axes = [axes]

for ax, g in zip(axes, groups):
    # Background: all cells in grey
    ax.scatter(umap_valid['UMAP_1'], umap_valid['UMAP_2'],
              s=0.3, alpha=0.05, color='lightgrey')
    # Foreground: current group
    sub = umap_valid[umap_valid['group'] == g]
    for ct in celltypes:
        ct_sub = sub[sub['celltype_assign'] == ct]
        ax.scatter(ct_sub['UMAP_1'], ct_sub['UMAP_2'], s=1, alpha=0.5,
                  color=ct_colors[ct], label=ct if g == groups[0] else '')
    ax.set_title(f'{g} (n={len(sub):,})', fontsize=13, fontweight='bold')
    ax.set_xlabel('UMAP_1')
    ax.set_ylabel('UMAP_2')

# Shared legend (from first panel)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='center right', bbox_to_anchor=(1.12, 0.5),
           markerscale=8, fontsize=9, title='Cell Type', title_fontsize=10)

plt.suptitle('UMAP Split by Group', fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
savefig('11_umap_split_by_group')


# In[16]:


# ── 3.1d UMAP — Split by Status ──
statuses = ['Hepatitis', 'Cirrhosis']
fig, axes = plt.subplots(1, len(statuses), figsize=(6 * len(statuses), 6))
if len(statuses) == 1:
    axes = [axes]

for ax, s in zip(axes, statuses):
    ax.scatter(umap_valid['UMAP_1'], umap_valid['UMAP_2'],
              s=0.3, alpha=0.05, color='lightgrey')
    sub = umap_valid[umap_valid['status'] == s]
    for ct in celltypes:
        ct_sub = sub[sub['celltype_assign'] == ct]
        ax.scatter(ct_sub['UMAP_1'], ct_sub['UMAP_2'], s=1, alpha=0.5,
                  color=ct_colors[ct], label=ct if s == statuses[0] else '')
    ax.set_title(f'{s} (n={len(sub):,})', fontsize=13, fontweight='bold')
    ax.set_xlabel('UMAP_1')
    ax.set_ylabel('UMAP_2')

# Shared legend (from first panel)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='center right', bbox_to_anchor=(1.12, 0.5),
           markerscale=8, fontsize=9, title='Cell Type', title_fontsize=10)

plt.suptitle('UMAP Split by Status', fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
savefig('12_umap_split_by_status')


# ### 3.2 Cell Type Proportions

# In[17]:


# ── 3.2 Cell Type Proportions ──
fig, axes = plt.subplots(1, 3, figsize=(20, 6))

# (A) Overall cell type distribution
ct_counts = df['celltype_assign'].value_counts().sort_values(ascending=True)
ct_counts.plot.barh(ax=axes[0], color=[ct_colors.get(ct, '#999') for ct in ct_counts.index])
axes[0].set_xlabel('Number of Cells')
axes[0].set_title('(A) Cell Type Distribution')
for i, (ct, v) in enumerate(zip(ct_counts.index, ct_counts.values)):
    axes[0].text(v + 50, i, f'{v:,}', va='center', fontsize=8)

# (B) Proportions by Group (stacked bar)
ct_group = pd.crosstab(df['group'], df['celltype_assign'], normalize='index') * 100
ct_group[celltypes].plot.bar(stacked=True, ax=axes[1],
                              color=[ct_colors.get(ct, '#999') for ct in celltypes],
                              width=0.7)
axes[1].set_ylabel('Proportion (%)')
axes[1].set_title('(B) Cell Type Proportion by Group')
axes[1].legend(fontsize=7, bbox_to_anchor=(1.02, 1.0), loc='upper left')
axes[1].tick_params(axis='x', rotation=0)

# (C) Proportions by Status
ct_status = pd.crosstab(df['status'], df['celltype_assign'], normalize='index') * 100
ct_status[celltypes].plot.bar(stacked=True, ax=axes[2],
                               color=[ct_colors.get(ct, '#999') for ct in celltypes],
                               width=0.7)
axes[2].set_ylabel('Proportion (%)')
axes[2].set_title('(C) Cell Type Proportion by Status')
axes[2].legend(fontsize=7, bbox_to_anchor=(1.02, 1.0), loc='upper left')
axes[2].tick_params(axis='x', rotation=0)

plt.tight_layout()
savefig('13_celltype_proportions')


# In[18]:


# ── 3.2b Cell Type Proportions per Sample ──
ct_sample = pd.crosstab(df['orig.ident'], df['celltype_assign'], normalize='index') * 100
sample_order = ['NC1_NT', 'NC1_PL', 'NC1_TC', 'NC2_TC', 'C1_TC']
sample_order = [s for s in sample_order if s in ct_sample.index]
ct_sample = ct_sample.loc[sample_order]

fig, ax = plt.subplots(figsize=(14, 6))
ct_sample[celltypes].plot.bar(stacked=True, ax=ax,
                               color=[ct_colors.get(ct, '#999') for ct in celltypes],
                               width=0.8)
ax.set_ylabel('Proportion (%)')
ax.set_title('Cell Type Proportion per Sample', fontsize=14, fontweight='bold')
ax.legend(fontsize=8, bbox_to_anchor=(1.02, 1.0), loc='upper left', title='Cell Type')
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)

plt.tight_layout()
savefig('14_celltype_proportion_per_sample')


# In[19]:


# ── 3.2c Fine Cell Type Proportions by Group ──
fine_group = pd.crosstab(df['celltype_assign_fine'], df['group'])
fine_group['total'] = fine_group.sum(axis=1)
fine_group = fine_group.sort_values('total', ascending=True)

fig, ax = plt.subplots(figsize=(10, 8))
fine_group.drop(columns='total').plot.barh(stacked=True, ax=ax,
    color=[group_color_map.get(g, '#999') for g in fine_group.drop(columns='total').columns],
    width=0.8)
ax.set_xlabel('Number of Cells')
ax.set_title('Fine Cell Type Distribution by Group', fontsize=14, fontweight='bold')
ax.legend(title='Group')

plt.tight_layout()
savefig('15_fine_celltype_by_group')


# ---
# ## 4. Multi-omics Integration — RNA + ATAC
# 
# Since the 10x Multiome platform captures RNA and ATAC from the **same cell**,  
# both modalities share the same cell barcodes.

# In[20]:


# ── 4.1 RNA vs ATAC per Cell ──
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# (A) nCount_RNA vs nCount_ATAC
for g in sorted(df['group'].unique()):
    sub = df[df['group'] == g]
    axes[0].scatter(sub['nCount_RNA'], sub['nCount_ATAC'],
                   s=1, alpha=0.2, color=group_color_map.get(g, '#999'), label=g)
axes[0].set_xlabel('nCount_RNA')
axes[0].set_ylabel('nCount_ATAC')
axes[0].set_title('(A) RNA vs ATAC Counts')
axes[0].legend(markerscale=5)

# (B) nFeature_RNA vs nFeature_ATAC
for g in sorted(df['group'].unique()):
    sub = df[df['group'] == g]
    axes[1].scatter(sub['nFeature_RNA'], sub['nFeature_ATAC'],
                   s=1, alpha=0.2, color=group_color_map.get(g, '#999'), label=g)
axes[1].set_xlabel('nFeature_RNA (Genes)')
axes[1].set_ylabel('nFeature_ATAC (Peaks)')
axes[1].set_title('(B) RNA Features vs ATAC Features')
axes[1].legend(markerscale=5)

# (C) RNA/ATAC ratio by cell type
df_temp = df.copy()
df_temp['rna_atac_ratio'] = df_temp['nCount_RNA'] / (df_temp['nCount_ATAC'] + 1)
ct_order = df_temp.groupby('celltype_assign')['rna_atac_ratio'].median().sort_values().index
sns.boxplot(data=df_temp, y='celltype_assign', x='rna_atac_ratio', ax=axes[2],
            order=ct_order, palette=[ct_colors.get(ct, '#999') for ct in ct_order],
            fliersize=0.3, linewidth=0.7)
axes[2].set_xlabel('RNA / ATAC Count Ratio')
axes[2].set_ylabel('')
axes[2].set_title('(C) RNA/ATAC Ratio by Cell Type')
axes[2].set_xlim(0, df_temp['rna_atac_ratio'].quantile(0.95))

plt.suptitle('Multi-omics: RNA-ATAC Joint QC', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
savefig('16_rna_atac_joint_qc')


# In[21]:


# ── 4.2 Multiome Data Summary Table ──
rna_gene = pd.read_csv(f'{MULTI}/rna_gene_metadata.csv', index_col=0)

print('=' * 60)
print('MULTIOME DATA SUMMARY')
print('=' * 60)
print(f"Total cells:              {len(df):,}")
print(f"Total RNA genes:          {len(rna_gene):,}")
print(f"  - Highly variable:      {rna_gene['highly_variable'].sum():,}")
print(f"Total ATAC peaks:         76,827")
print(f"Samples:                  {df['orig.ident'].nunique()}")
print(f"Groups:                   {', '.join(sorted(df['group'].unique()))}")
print(f"Status:                   {', '.join(sorted(df['status'].unique()))}")
print(f"Major cell types:         {df['celltype_assign'].nunique()}")
print(f"Fine cell types:          {df['celltype_assign_fine'].nunique()}")
print(f"Median RNA UMIs/cell:     {df['nCount_RNA'].median():.0f}")
print(f"Median RNA genes/cell:    {df['nFeature_RNA'].median():.0f}")
print(f"Median ATAC frags/cell:   {df['nCount_ATAC'].median():.0f}")
print(f"Median ATAC peaks/cell:   {df['nFeature_ATAC'].median():.0f}")


# ---
# ## 5. GeoMx Spatial Transcriptomics Overview
# 
# NanoString GeoMx DSP — Whole Transcriptome Atlas (WTA)  
# Q3-normalized expression data with CD31+ / CD31- segment distinction.

# In[22]:


# ── 5.1 GeoMx Sample Summary ──
print('=' * 60)
print('GeoMx DATA SUMMARY')
print('=' * 60)
print(f"Total segments:       {len(geo_meta)}")
print(f"Total genes (Q3):     {len(geo_expr):,}")
print(f"Slides:               {geo_meta['SlideName'].nunique()} ({', '.join(geo_meta['SlideName'].unique())})")
print(f"ROIs:                 {geo_meta['ROILabel'].nunique()}")
print()

# Group breakdown
print('Group breakdown:')
for col in ['Group_NTTC', 'Group_CD31', 'Group_Cirr', 'Group_detailed']:
    print(f"\n  {col}:")
    print(geo_meta[col].value_counts().to_string().replace('\n', '\n    '))


# In[23]:


# ── 5.2 GeoMx Segment Composition ──
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# (A) NT vs TC
geo_meta['Group_NTTC'].value_counts().plot.bar(
    ax=axes[0], color=['#4C72B0', '#DD8452'])
axes[0].set_title('NT vs TC')
axes[0].set_ylabel('Number of Segments')
axes[0].tick_params(axis='x', rotation=0)

# (B) CD31+ vs CD31-
geo_meta['Group_CD31'].value_counts().plot.bar(
    ax=axes[1], color=['#55A868', '#C44E52'])
axes[1].set_title('CD31+ vs CD31-')
axes[1].set_ylabel('Number of Segments')
axes[1].tick_params(axis='x', rotation=0)

# (C) Detailed groups
geo_meta['Group_detailed'].value_counts().sort_index().plot.bar(
    ax=axes[2], color='steelblue')
axes[2].set_title('Detailed Groups')
axes[2].set_ylabel('Number of Segments')
axes[2].set_xticklabels(axes[2].get_xticklabels(), rotation=45, ha='right', fontsize=8)

plt.suptitle('GeoMx — Segment Composition', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
savefig('17_geomx_segment_composition')


# In[24]:


# ── 5.3 GeoMx QC Metrics ──
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# (A) AOI Nuclei Count
sns.boxplot(data=geo_meta, x='Group_NTTC', y='AOINucleiCount', ax=axes[0],
            palette=['#4C72B0', '#DD8452'])
axes[0].set_title('AOI Nuclei Count')

# (B) Sequencing Saturation
sns.boxplot(data=geo_meta, x='Group_NTTC', y='SequencingSaturation', ax=axes[1],
            palette=['#4C72B0', '#DD8452'])
axes[1].set_title('Sequencing Saturation')

# (C) AOI Surface Area by CD31
sns.boxplot(data=geo_meta, x='Group_CD31', y='AOISurfaceArea', ax=axes[2],
            palette=['#55A868', '#C44E52'])
axes[2].set_title('AOI Surface Area by CD31')

plt.suptitle('GeoMx — QC Metrics', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
savefig('18_geomx_qc_metrics')

print('\n' + '=' * 60)
print(f'All figures saved to: {SAVE_DIR}')
print('=' * 60)


# ---
# ## 6. Future Plan — 향후 분석 계획
# 
# ### 6.1 Differential Expression Analysis (DEG)
# | Comparison | Method | Description |
# |-----------|--------|-------------|
# | TC vs NT (Multiome) | Wilcoxon / MAST / pseudobulk DESeq2 | 세포 타입별 종양 vs 정상 DEG |
# | CD31+ vs CD31- (GeoMx) | limma-trend | 내피세포 구획 vs 비내피세포 구획 DEG |
# | Cirrhosis+ vs Cirrhosis- | limma / MAST | 간경변 유무에 따른 DEG |
# 
# ### 6.2 Chromatin Accessibility (ATAC-seq)
# - **Peak Calling & Differential Accessibility**: 세포 타입별, 조건별로 오픈 크로마틴 영역 차이 분석
# - **Motif Enrichment**: 전사인자 결합 모티프 분석 (chromVAR)
# - **Peak-to-Gene Linkage**: ATAC peak과 RNA 발현의 상관관계를 통해 *cis*-regulatory element 추론
# 
# ### 6.3 Gene Regulatory Network (GRN)
# - RNA + ATAC 데이터를 통합하여 전사인자 → 타겟 유전자 조절 네트워크 구축
# - Tools: SCENIC+ / FigR / Signac
# 
# ### 6.4 Multi-modal Integration
# - **WNN (Weighted Nearest Neighbor)**: RNA + ATAC joint embedding
# - **Spatial ↔ Single-cell Deconvolution**: GeoMx ROI의 세포 구성을 Multiome reference로 추정
# - **Cross-platform Validation**: Multiome DEG 결과를 GeoMx spatial data에서 검증
# 
# ### 6.5 Trajectory & Cell Communication
# - **Pseudotime Analysis**: Monocle3 — 종양 미세환경 내 세포 분화 궤적
# - **Cell-Cell Communication**: CellChat / NicheNet — 리간드-수용체 상호작용 분석
# 
# ### 6.6 Clinical Relevance
# - HCC subtype 분류와 예후 연관 분석
# - 간경변(Cirrhosis) 배경에 따른 종양 미세환경 차이 규명

# ---
# ## Summary Table
# 
# | Item | Detail |
# |------|--------|
# | **Data Type** | 10x Multiome (scRNA + scATAC) + NanoString GeoMx (Spatial WTA) |
# | **Disease** | Hepatocellular Carcinoma (HCC) |
# | **Multiome Cells** | 28,859 cells |
# | **Multiome Genes** | 36,601 (RNA) / 76,827 peaks (ATAC) |
# | **Multiome Groups** | NT (Non-Tumor), TC (Tumor Core), PL (Peripheral Liver) |
# | **Multiome Status** | Cirrhosis, Hepatitis |
# | **Cell Types** | 9 major / ~20 fine subtypes |
# | **Doublet Detection** | scDblFinder — singlets retained |
# | **GeoMx Segments** | 88 segments (CD31+ / CD31-) |
# | **GeoMx Genes** | 16,818 (Q3-normalized WTA) |
# | **GeoMx Groups** | NT/TC × CD31+/CD31- × Cirrhosis+/- |

# 

# 
