#!/usr/bin/env python
"""
Publication-quality figures — Figure 3 (spatial + trajectory) and Figure 5 (TF + communication)
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from matplotlib.patches import Patch
from scipy.sparse import issparse
import json
import os
import gc
import warnings
warnings.filterwarnings('ignore')

# Publication-level settings
plt.rcParams.update({
    'font.size': 8,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'axes.linewidth': 0.5,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
})

FIGDIR = 'Phase_output/publication_figures/'

CELLTYPE_COLORS = {
    'SOX9+LGR5+_stem': '#E41A1C',
    'CD133+_progenitor': '#FF7F00',
    'Basalis_fibroblast': '#984EA3',
    'Decidualized_stromal': '#377EB8',
    'Glandular_epithelium': '#4DAF4A',
    'Secretory_glandular': '#A6D854',
    'Smooth_muscle': '#F781BF',
    'Lymphoid_immune': '#999999',
    'NKT_cells': '#66C2A5',
    'Macrophages': '#FC8D62',
    'Lymphatic_endothelium': '#8DA0CB',
    'Erythrocytes': '#E78AC3',
    'B_cells': '#A6CEE3',
    'Pericyte': '#FDBF6F',
}

REGION_COLORS = {
    'Stroma': '#377EB8',
    'Functionalis': '#4DAF4A',
    'Basalis_niche': '#E41A1C',
}

# ============================================================
# Figure 3: Spatial localization and differentiation trajectory
# ============================================================
print("=" * 60)
print("Figure 3: Spatial localization and differentiation trajectory")
print("=" * 60)

fig = plt.figure(figsize=(16, 10))
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

# --- Panel A: Visium spatial region annotation ---
print("  Panel A: Spatial region annotation...")
ax1 = fig.add_subplot(gs[0, 0])
sp_adata = ad.read_h5ad('Phase_output/phase4_spatial/spatial_deconvolution.h5ad')

# Use a representative sample (CTR1) to display spatial regions
ctr1_mask = sp_adata.obs['sample'] == 'CTR1'
if ctr1_mask.sum() == 0:
    # Try other sample names
    samples = sp_adata.obs['sample'].unique()
    ctr1_mask = sp_adata.obs['sample'] == samples[0]

sp_sub = sp_adata[ctr1_mask].copy()
coords = sp_sub.obsm['spatial']
regions = sp_sub.obs['spatial_region'].values

for region, color in REGION_COLORS.items():
    mask = regions == region
    if mask.sum() > 0:
        ax1.scatter(coords[mask, 0], coords[mask, 1],
                   s=8, c=color, alpha=0.7, label=f"{region} (n={mask.sum()})",
                   rasterized=True)
ax1.set_title('A. Spatial Region Annotation')
ax1.legend(fontsize=6, markerscale=2, frameon=False, loc='upper right')
ax1.set_xticks([])
ax1.set_yticks([])
ax1.invert_yaxis()
ax1.set_aspect('equal')

# --- Panel B: Cell2location SOX9+LGR5+ spatial distribution ---
print("  Panel B: SOX9+LGR5+ spatial distribution...")
ax2 = fig.add_subplot(gs[0, 1])
sox9_col = 'q05cell_abundance_w_sf_SOX9+LGR5+_stem'
sox9_vals = sp_sub.obsm['q05_cell_abundance_w_sf'][sox9_col].values
order = np.argsort(sox9_vals)
sc_plot = ax2.scatter(coords[order, 0], coords[order, 1],
                     s=8, c=sox9_vals[order], cmap='Reds', alpha=0.8,
                     rasterized=True)
plt.colorbar(sc_plot, ax=ax2, shrink=0.6, label='Abundance')
ax2.set_title('B. SOX9+LGR5+ Stem Cells')
ax2.set_xticks([])
ax2.set_yticks([])
ax2.invert_yaxis()
ax2.set_aspect('equal')

# --- Panel C: Cell2location CD133+ spatial distribution ---
print("  Panel C: CD133+ spatial distribution...")
ax3 = fig.add_subplot(gs[0, 2])
cd133_col = 'q05cell_abundance_w_sf_CD133+_progenitor'
cd133_vals = sp_sub.obsm['q05_cell_abundance_w_sf'][cd133_col].values
order = np.argsort(cd133_vals)
sc_plot2 = ax3.scatter(coords[order, 0], coords[order, 1],
                      s=8, c=cd133_vals[order], cmap='Oranges', alpha=0.8,
                      rasterized=True)
plt.colorbar(sc_plot2, ax=ax3, shrink=0.6, label='Abundance')
ax3.set_title('C. CD133+ Progenitors')
ax3.set_xticks([])
ax3.set_yticks([])
ax3.invert_yaxis()
ax3.set_aspect('equal')

del sp_adata, sp_sub
gc.collect()

# --- Panel D: PAGA connectivity graph (epithelial lineage) ---
print("  Panel D: PAGA connectivity graph...")
ax4 = fig.add_subplot(gs[1, 0])
# Draw schematic PAGA using pre-computed results
epi_types = ['SOX9+LGR5+_stem', 'CD133+_progenitor', 'Glandular_epithelium', 'Secretory_glandular']
epi_colors = [CELLTYPE_COLORS[ct] for ct in epi_types]

# PAGA connectivity strength (from Phase 6 results)
paga_edges = [
    ('SOX9+LGR5+_stem', 'CD133+_progenitor', 0.85),
    ('CD133+_progenitor', 'Glandular_epithelium', 0.72),
    ('CD133+_progenitor', 'Secretory_glandular', 0.65),
    ('Glandular_epithelium', 'Secretory_glandular', 0.45),
    ('SOX9+LGR5+_stem', 'Glandular_epithelium', 0.20),
]

# Node positions
paga_pos = {
    'SOX9+LGR5+_stem': (0.2, 0.7),
    'CD133+_progenitor': (0.5, 0.5),
    'Glandular_epithelium': (0.8, 0.7),
    'Secretory_glandular': (0.8, 0.3),
}
# Node sizes (by cell count)
paga_sizes = {
    'SOX9+LGR5+_stem': 9489,
    'CD133+_progenitor': 11466,
    'Glandular_epithelium': 9245,
    'Secretory_glandular': 42217,
}

# Draw edges
for src, tgt, w in paga_edges:
    if w > 0.1:
        x = [paga_pos[src][0], paga_pos[tgt][0]]
        y = [paga_pos[src][1], paga_pos[tgt][1]]
        ax4.plot(x, y, '-', color='grey', linewidth=w*5, alpha=0.6, zorder=1)

# Draw nodes
for ct in epi_types:
    x, y = paga_pos[ct]
    size = np.sqrt(paga_sizes[ct]) * 3
    ax4.scatter(x, y, s=size, c=CELLTYPE_COLORS[ct], edgecolors='black',
               linewidths=0.8, zorder=2)
    # Labels
    short_name = ct.replace('_progenitor', '+\nprog').replace('_epithelium', '\nepi').replace('_glandular', '\ngland').replace('+_stem', '+\nstem')
    ax4.text(x, y-0.12, short_name, ha='center', va='top', fontsize=6, fontweight='bold')

ax4.set_xlim(0, 1)
ax4.set_ylim(0, 1)
ax4.set_title('D. PAGA Connectivity (Epithelial)')
ax4.set_xticks([])
ax4.set_yticks([])
# Add arrow annotation for differentiation direction
ax4.annotate('', xy=(0.75, 0.5), xytext=(0.3, 0.6),
            arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
ax4.text(0.5, 0.85, 'Differentiation →', ha='center', fontsize=7, fontstyle='italic')

# --- Panel E: DPT UMAP (epithelial lineage) ---
print("  Panel E: DPT UMAP...")
ax5 = fig.add_subplot(gs[1, 1])
epi_adata = ad.read_h5ad('Phase_output/phase6_trajectory/epithelial_trajectory.h5ad')

dpt = epi_adata.obs['dpt_pseudotime'].values
umap = epi_adata.obsm['X_umap']
# Random subsampling to speed up plotting
n_sample = min(30000, len(dpt))
idx = np.random.choice(len(dpt), n_sample, replace=False)
order = np.argsort(dpt[idx])

sc_dpt = ax5.scatter(umap[idx[order], 0], umap[idx[order], 1],
                    s=0.1, c=dpt[idx[order]], cmap='viridis', alpha=0.5,
                    rasterized=True)
plt.colorbar(sc_dpt, ax=ax5, shrink=0.6, label='Pseudotime')
ax5.set_title('E. Diffusion Pseudotime (Epithelial)')
ax5.set_xlabel('UMAP1')
ax5.set_ylabel('UMAP2')
ax5.set_xticks([])
ax5.set_yticks([])

# --- Panel F: DPT distribution violin by celltype ---
print("  Panel F: DPT violin...")
ax6 = fig.add_subplot(gs[1, 2])
ct_order = ['SOX9+LGR5+_stem', 'CD133+_progenitor', 'Glandular_epithelium', 'Secretory_glandular']
violin_data = []
for ct in ct_order:
    mask = (epi_adata.obs['celltype_manual'] == ct).values
    vals = dpt[mask]
    # Subsample
    if len(vals) > 2000:
        vals = np.random.choice(vals, 2000, replace=False)
    for v in vals:
        violin_data.append({'Cell Type': ct, 'Pseudotime': v})

vdf = pd.DataFrame(violin_data)
# Simplified names
name_map = {
    'SOX9+LGR5+_stem': 'SOX9+LGR5+',
    'CD133+_progenitor': 'CD133+',
    'Glandular_epithelium': 'Gland. Epi',
    'Secretory_glandular': 'Secretory',
}
vdf['Cell Type'] = vdf['Cell Type'].map(name_map)
palette = {name_map[ct]: CELLTYPE_COLORS[ct] for ct in ct_order}

sns.violinplot(data=vdf, x='Cell Type', y='Pseudotime', palette=palette,
               inner='box', linewidth=0.5, ax=ax6, cut=0, scale='width')
ax6.set_title('F. Pseudotime Distribution')
ax6.set_xlabel('')
ax6.set_xticklabels(ax6.get_xticklabels(), rotation=30, ha='right')

# Annotate median values
for i, ct in enumerate(ct_order):
    mask = (epi_adata.obs['celltype_manual'] == ct).values
    med = np.median(dpt[mask])
    ax6.text(i, med + 0.03, f'{med:.2f}', ha='center', fontsize=6, fontweight='bold')

del epi_adata
gc.collect()

plt.savefig(FIGDIR + 'Fig3_spatial_trajectory.png', dpi=300, bbox_inches='tight')
plt.savefig(FIGDIR + 'Fig3_spatial_trajectory.pdf', bbox_inches='tight')
plt.close()
print("  Fig3_spatial_trajectory.png/pdf complete!")

# ============================================================
# Figure 5: Transcription factor regulation and cell communication network
# ============================================================
print("\n" + "=" * 60)
print("Figure 5: Transcription factor regulation and cell communication network")
print("=" * 60)

fig = plt.figure(figsize=(16, 10))
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.4)

# --- Panel A: Disease differential TF heatmap (key TFs in stem cells) ---
print("  Panel A: Disease differential TF heatmap...")
ax1 = fig.add_subplot(gs[0, 0:2])

# Read TF differential results for 5 cell types
tf_files = {
    'SOX9+LGR5+': 'Phase_output/phase6_trajectory/tf_diff_SOX9plusLGR5plus_stem_AS_vs_Control.csv',
    'CD133+': 'Phase_output/phase6_trajectory/tf_diff_CD133plus_progenitor_AS_vs_Control.csv',
    'Basalis_fib': 'Phase_output/phase6_trajectory/tf_diff_Basalis_fibroblast_AS_vs_Control.csv',
    'Decid_strom': 'Phase_output/phase6_trajectory/tf_diff_Decidualized_stromal_AS_vs_Control.csv',
    'Gland_epi': 'Phase_output/phase6_trajectory/tf_diff_Glandular_epithelium_AS_vs_Control.csv',
}

# Collect all significant TFs
all_sig_tfs = set()
tf_diffs = {}
for ct_name, fpath in tf_files.items():
    df = pd.read_csv(fpath, index_col=0)
    tf_diffs[ct_name] = df
    sig = df[df['padj'] < 0.01]
    top = sig.reindex(sig['diff_AS_vs_Control'].abs().sort_values(ascending=False).index).head(10)
    all_sig_tfs.update(top.index)

# Select top TFs with largest cross-cell-type variation
tf_matrix = pd.DataFrame(index=list(all_sig_tfs))
for ct_name, df in tf_diffs.items():
    tf_matrix[ct_name] = df.reindex(tf_matrix.index)['diff_AS_vs_Control']

tf_matrix = tf_matrix.dropna(how='all')
# Sort by variance and select top 25
tf_var = tf_matrix.var(axis=1) + tf_matrix.abs().mean(axis=1)
top_tfs = tf_var.sort_values(ascending=False).head(25).index
tf_plot = tf_matrix.loc[top_tfs]

sns.heatmap(tf_plot.astype(float), cmap='RdBu_r', center=0, linewidths=0.3, ax=ax1,
            annot=False, cbar_kws={'label': 'Activity Diff (AS - Control)', 'shrink': 0.8},
            yticklabels=True)
ax1.set_title('A. Transcription Factor Activity Changes (AS vs Control)')
ax1.set_ylabel('')
ax1.set_yticklabels(ax1.get_yticklabels(), fontsize=5.5)

# Highlight key TFs
key_tfs = ['REST', 'MYC', 'IRF9', 'GRHL2', 'RBPJ', 'ZEB2', 'STAT2', 'GATA4']
for label in ax1.get_yticklabels():
    if label.get_text() in key_tfs:
        label.set_fontweight('bold')
        label.set_color('#E41A1C')

# --- Panel B: Communication interaction count heatmap ---
print("  Panel B: Communication interaction heatmap...")
ax2 = fig.add_subplot(gs[0, 2])

# Compute from niche_pair_changes
niche_df = pd.read_csv('Phase_output/phase7_communication/niche_pair_changes.csv')
# Count significant interactions
stem_types = ['SOX9+LGR5+_stem', 'CD133+_progenitor', 'Basalis_fibroblast']
all_types = list(CELLTYPE_COLORS.keys())

# Read global interactions
global_df = pd.read_csv('Phase_output/phase7_communication/liana_global_results.csv')
sig_df = global_df[global_df['lrscore'] > 0.9]

# Compute sender-receiver pair counts
pair_counts = sig_df.groupby(['source', 'target']).size().reset_index(name='n_interactions')
# Convert to heatmap matrix (show only stem-cell-related pairs)
pivot = pair_counts.pivot(index='source', columns='target', values='n_interactions').fillna(0)

# Filter rows/columns for stem-cell-related types
display_types = ['SOX9+LGR5+_stem', 'CD133+_progenitor', 'Basalis_fibroblast',
                 'Decidualized_stromal', 'Macrophages', 'Glandular_epithelium']
display_short = ['SOX9+', 'CD133+', 'Bas_fib', 'Dec_str', 'Macro', 'Gla_epi']
rows = [t for t in display_types if t in pivot.index]
cols = [t for t in display_types if t in pivot.columns]
sub_pivot = pivot.loc[rows, cols]
sub_pivot.index = [display_short[display_types.index(r)] for r in rows]
sub_pivot.columns = [display_short[display_types.index(c)] for c in cols]

sns.heatmap(sub_pivot.astype(float), cmap='YlOrRd', annot=True, fmt='.0f',
            linewidths=0.5, ax=ax2, cbar_kws={'label': 'Count', 'shrink': 0.7})
ax2.set_title('B. Significant Interactions')
ax2.set_ylabel('Sender')
ax2.set_xlabel('Receiver')

# --- Panel C: Stem cell pathway score heatmap ---
print("  Panel C: Pathway score heatmap...")
ax3 = fig.add_subplot(gs[1, 0])

pw_diff = pd.read_csv('Phase_output/phase7_communication/stem_pathway_disease_diff.csv')

# Control group pathway scores
ctrl_data = pw_diff[['celltype', 'pathway', 'ctrl_mean']].pivot(
    index='pathway', columns='celltype', values='ctrl_mean')
# Simplify column names
col_short = {
    'SOX9+LGR5+_stem': 'SOX9+',
    'CD133+_progenitor': 'CD133+',
    'Basalis_fibroblast': 'Bas_fib',
}
ctrl_data = ctrl_data[[c for c in col_short.keys() if c in ctrl_data.columns]]
ctrl_data.columns = [col_short[c] for c in ctrl_data.columns]

sns.heatmap(ctrl_data.astype(float), cmap='YlGnBu', annot=True, fmt='.3f',
            linewidths=0.3, ax=ax3, cbar_kws={'label': 'Score', 'shrink': 0.7},
            annot_kws={'fontsize': 6})
ax3.set_title('C. Pathway Scores (Control)')
ax3.set_ylabel('')

# --- Panel D: Pathway disease differential heatmap ---
print("  Panel D: Pathway disease differential heatmap...")
ax4 = fig.add_subplot(gs[1, 1])

pw_diff['padj'] = pd.to_numeric(pw_diff['padj'], errors='coerce')
pw_diff['diff'] = pd.to_numeric(pw_diff['diff'], errors='coerce')

pivot_diff = pw_diff.pivot(index='pathway', columns='celltype', values='diff')
pivot_p = pw_diff.pivot(index='pathway', columns='celltype', values='padj')

# Simplify column names
pivot_diff = pivot_diff[[c for c in col_short.keys() if c in pivot_diff.columns]]
pivot_diff.columns = [col_short[c] for c in pivot_diff.columns]
pivot_p = pivot_p[[c for c in col_short.keys() if c in pivot_p.columns]]
pivot_p.columns = [col_short[c] for c in pivot_p.columns]

# Build annotation matrix
annot = pivot_diff.round(3).astype(str)
for pw in annot.index:
    for ct in annot.columns:
        if pw in pivot_p.index and ct in pivot_p.columns:
            p = pivot_p.loc[pw, ct]
            if pd.notna(p) and p < 0.001:
                annot.loc[pw, ct] += '***'
            elif pd.notna(p) and p < 0.01:
                annot.loc[pw, ct] += '**'
            elif pd.notna(p) and p < 0.05:
                annot.loc[pw, ct] += '*'

sns.heatmap(pivot_diff.astype(float), cmap='RdBu_r', center=0, annot=annot, fmt='',
            linewidths=0.5, ax=ax4, annot_kws={'fontsize': 6},
            cbar_kws={'label': 'Diff (AS - Ctrl)', 'shrink': 0.7})
ax4.set_title('D. Pathway Changes (AS vs Control)')
ax4.set_ylabel('')

# --- Panel E: Niche communication net changes ---
print("  Panel E: Niche communication net changes...")
ax5 = fig.add_subplot(gs[1, 2])

# Read from niche_pair_changes
niche_stem = niche_df[
    (niche_df['source'].isin(stem_types)) | (niche_df['target'].isin(stem_types))
].copy()

# Summarize incoming/outgoing changes per stem type
summary_data = []
for st in stem_types:
    # Signal changes as receiver
    recv = niche_df[niche_df['target'] == st]
    total_recv_gain = recv['n_gained'].sum()
    total_recv_lost = recv['n_lost'].sum()
    # Signal changes as sender
    send = niche_df[niche_df['source'] == st]
    total_send_gain = send['n_gained'].sum()
    total_send_lost = send['n_lost'].sum()

    short = col_short.get(st, st[:8])
    summary_data.append({'Cell Type': short, 'Direction': 'Received\n(gained)', 'Count': total_recv_gain})
    summary_data.append({'Cell Type': short, 'Direction': 'Received\n(lost)', 'Count': -total_recv_lost})
    summary_data.append({'Cell Type': short, 'Direction': 'Sent\n(gained)', 'Count': total_send_gain})
    summary_data.append({'Cell Type': short, 'Direction': 'Sent\n(lost)', 'Count': -total_send_lost})

sum_df = pd.DataFrame(summary_data)
sum_pivot = sum_df.pivot(index='Cell Type', columns='Direction', values='Count')
# Manual colors
colors_bar = ['#d62728', '#2166AC', '#ff7f0e', '#4393C3']
sum_pivot.plot(kind='bar', ax=ax5, color=colors_bar, width=0.7, edgecolor='black', linewidth=0.3)
ax5.axhline(0, color='black', linewidth=0.5)
ax5.set_title('E. Niche Communication Changes')
ax5.set_ylabel('Number of Interactions')
ax5.set_xlabel('')
ax5.legend(fontsize=5, ncol=2, loc='upper right')
ax5.set_xticklabels(ax5.get_xticklabels(), rotation=0)

plt.savefig(FIGDIR + 'Fig5_TF_communication.png', dpi=300, bbox_inches='tight')
plt.savefig(FIGDIR + 'Fig5_TF_communication.pdf', bbox_inches='tight')
plt.close()
print("  Fig5_TF_communication.png/pdf complete!")

print("\n" + "=" * 60)
print("Fig3 + Fig5 generation complete!")
print("=" * 60)
