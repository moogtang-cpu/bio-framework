#!/usr/bin/env python
"""
Post-Analysis Step 4: Publication-quality figures — Figure 1-6
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

# Nature-style color palette
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

DATASET_COLORS = {
    'GSE215968_AS_vs_WOI': '#1f77b4',
    'GSE215968_CD133': '#ff7f0e',
    'GSE111976': '#2ca02c',
    'E-MTAB-10287': '#d62728',
    'GSE260658': '#9467bd',
}

print("=" * 60)
print("Generating publication-quality figures")
print("=" * 60)

# Load data
print("\nLoading data...")
adata = ad.read_h5ad('Phase_output/phase3_integration/integrated_harmony.h5ad')
print(f"  Main atlas: {adata.shape}")

# ============================================================
# Figure 1: Human endometrial comprehensive atlas
# ============================================================
print("\n--- Figure 1: Comprehensive atlas ---")

fig = plt.figure(figsize=(14, 10))
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3)

# Panel A: UMAP by celltype
ax1 = fig.add_subplot(gs[0, 0:2])
celltypes = adata.obs['celltype_manual'].values
for ct in sorted(CELLTYPE_COLORS.keys()):
    mask = celltypes == ct
    if mask.sum() > 0:
        ax1.scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
                   s=0.1, alpha=0.3, c=CELLTYPE_COLORS[ct], label=ct, rasterized=True)
ax1.set_xlabel('UMAP1')
ax1.set_ylabel('UMAP2')
ax1.set_title('A. Cell Type Annotation (314,805 cells)')
ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=6,
           markerscale=10, frameon=False)
ax1.set_xticks([])
ax1.set_yticks([])

# Panel B: UMAP by dataset
ax2 = fig.add_subplot(gs[0, 2])
datasets = adata.obs['dataset'].values
for ds, color in DATASET_COLORS.items():
    mask = datasets == ds
    if mask.sum() > 0:
        ax2.scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
                   s=0.05, alpha=0.2, c=color, label=ds.replace('GSE215968_', '').replace('GSE', 'G'),
                   rasterized=True)
ax2.set_title('B. Dataset')
ax2.legend(fontsize=5, markerscale=20, frameon=False)
ax2.set_xticks([])
ax2.set_yticks([])

# Panel C: Cell type proportions by dataset
ax3 = fig.add_subplot(gs[1, 0])
ct_prop = pd.crosstab(adata.obs['dataset'], adata.obs['celltype_manual'], normalize='index')
ct_order = ['SOX9+LGR5+_stem', 'CD133+_progenitor', 'Basalis_fibroblast',
            'Decidualized_stromal', 'Glandular_epithelium', 'Secretory_glandular',
            'Smooth_muscle', 'Lymphoid_immune', 'NKT_cells', 'Macrophages',
            'Lymphatic_endothelium', 'Erythrocytes', 'B_cells', 'Pericyte']
ct_order = [c for c in ct_order if c in ct_prop.columns]
ct_prop = ct_prop[ct_order]
ct_prop.plot(kind='barh', stacked=True, ax=ax3,
             color=[CELLTYPE_COLORS.get(c, '#888') for c in ct_order],
             legend=False, width=0.7)
ax3.set_xlabel('Proportion')
ax3.set_title('C. Cell Type Distribution')
ax3.set_yticklabels([d.replace('GSE215968_', '').replace('E-MTAB-', 'EMTAB') for d in ct_prop.index],
                     fontsize=6)

# Panel D: Lineage pie
ax4 = fig.add_subplot(gs[1, 1])
lineage_map = {
    'SOX9+LGR5+_stem': 'Stem/Progenitor', 'CD133+_progenitor': 'Stem/Progenitor',
    'Basalis_fibroblast': 'Mesenchymal', 'Decidualized_stromal': 'Mesenchymal',
    'Smooth_muscle': 'Mesenchymal', 'Pericyte': 'Mesenchymal',
    'Glandular_epithelium': 'Epithelial', 'Secretory_glandular': 'Epithelial',
    'Lymphoid_immune': 'Immune', 'NKT_cells': 'Immune',
    'Macrophages': 'Immune', 'B_cells': 'Immune', 'Erythrocytes': 'Immune',
    'Lymphatic_endothelium': 'Endothelial',
}
lineage_counts = adata.obs['celltype_manual'].map(lineage_map).value_counts()
lineage_colors = {'Mesenchymal': '#377EB8', 'Immune': '#999999', 'Epithelial': '#4DAF4A',
                  'Stem/Progenitor': '#E41A1C', 'Endothelial': '#8DA0CB'}
ax4.pie(lineage_counts.values,
        labels=[f"{l}\n({c/sum(lineage_counts)*100:.1f}%)" for l, c in lineage_counts.items()],
        colors=[lineage_colors.get(l, '#888') for l in lineage_counts.index],
        startangle=90, textprops={'fontsize': 7})
ax4.set_title('D. Lineage Distribution')

# Panel E: Cell counts bar
ax5 = fig.add_subplot(gs[1, 2])
ct_counts = adata.obs['celltype_manual'].value_counts()
ct_counts = ct_counts.reindex([c for c in ct_order if c in ct_counts.index])
bars = ax5.barh(range(len(ct_counts)), ct_counts.values,
                color=[CELLTYPE_COLORS.get(c, '#888') for c in ct_counts.index],
                edgecolor='black', linewidth=0.3)
ax5.set_yticks(range(len(ct_counts)))
ax5.set_yticklabels(ct_counts.index, fontsize=6)
ax5.invert_yaxis()
ax5.set_xlabel('Number of cells')
ax5.set_title('E. Cell Counts')
for bar, v in zip(bars, ct_counts.values):
    ax5.text(bar.get_width() + 500, bar.get_y() + bar.get_height()/2,
             f'{v:,}', va='center', fontsize=5)

plt.savefig(FIGDIR + 'Fig1_atlas_overview.png', dpi=300, bbox_inches='tight')
plt.savefig(FIGDIR + 'Fig1_atlas_overview.pdf', bbox_inches='tight')
plt.close()
print("  Fig1_atlas_overview.png/pdf")

# ============================================================
# Figure 2: Stem/progenitor cell subpopulation identification
# ============================================================
print("\n--- Figure 2: Stem cell subpopulations ---")

stem_markers = ['SOX9', 'LGR5', 'PROM1', 'CD44', 'KRT19', 'AXIN2', 'MKI67', 'TOP2A']
available_markers = [m for m in stem_markers if m in adata.var_names]

fig = plt.figure(figsize=(14, 8))
gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.4, wspace=0.35)

# Panel A-D: Feature plots
for idx, gene in enumerate(available_markers[:4]):
    ax = fig.add_subplot(gs[0, idx])
    gene_idx = list(adata.var_names).index(gene)
    x = adata.X[:, gene_idx]
    if issparse(x):
        x = x.toarray().flatten()
    else:
        x = np.asarray(x).flatten()

    # Sort for plotting (high expression on top)
    order = np.argsort(x)
    ax.scatter(adata.obsm['X_umap'][order, 0], adata.obsm['X_umap'][order, 1],
               s=0.05, c=x[order], cmap='viridis', alpha=0.3, rasterized=True)
    ax.set_title(f'{chr(65+idx)}. {gene}')
    ax.set_xticks([])
    ax.set_yticks([])

# Panel E: Stem cell UMAP zoom
ax_zoom = fig.add_subplot(gs[1, 0:2])
stem_types = ['SOX9+LGR5+_stem', 'CD133+_progenitor', 'Basalis_fibroblast']
for ct in stem_types:
    mask = adata.obs['celltype_manual'] == ct
    ax_zoom.scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
                   s=0.3, alpha=0.5, c=CELLTYPE_COLORS[ct], label=ct, rasterized=True)
# Non-stem cells in grey
other_mask = ~adata.obs['celltype_manual'].isin(stem_types)
ax_zoom.scatter(adata.obsm['X_umap'][other_mask, 0], adata.obsm['X_umap'][other_mask, 1],
               s=0.02, alpha=0.05, c='lightgrey', rasterized=True)
ax_zoom.set_title('E. Stem/Progenitor Populations')
ax_zoom.legend(markerscale=10, fontsize=7, frameon=False)
ax_zoom.set_xticks([])
ax_zoom.set_yticks([])

# Panel F: Violin plot for stem markers
ax_violin = fig.add_subplot(gs[1, 2:4])
violin_genes = [g for g in ['SOX9', 'LGR5', 'PROM1', 'CD44'] if g in adata.var_names]
stem_data = []
for ct in stem_types + ['Glandular_epithelium', 'Decidualized_stromal']:
    ct_mask = (adata.obs['celltype_manual'] == ct).values
    for gene in violin_genes:
        gene_idx = list(adata.var_names).index(gene)
        x = adata.X[ct_mask, gene_idx]
        if issparse(x):
            x = x.toarray().flatten()
        else:
            x = np.asarray(x).flatten()
        # Subsample
        if len(x) > 500:
            x = np.random.choice(x, 500, replace=False)
        for v in x:
            stem_data.append({'celltype': ct, 'gene': gene, 'expr': v})

stem_df = pd.DataFrame(stem_data)
# Simplified as heatmap (mean expression)
pivot = stem_df.groupby(['celltype', 'gene'])['expr'].mean().reset_index()
pivot = pivot.pivot(index='gene', columns='celltype', values='expr')
sns.heatmap(pivot, cmap='YlOrRd', annot=True, fmt='.2f', linewidths=0.5, ax=ax_violin,
            cbar_kws={'label': 'Mean Expression'})
ax_violin.set_title('F. Stem Marker Expression')
ax_violin.set_ylabel('')

plt.savefig(FIGDIR + 'Fig2_stem_cell_identification.png', dpi=300, bbox_inches='tight')
plt.savefig(FIGDIR + 'Fig2_stem_cell_identification.pdf', bbox_inches='tight')
plt.close()
print("  Fig2_stem_cell_identification.png/pdf")

del adata
gc.collect()

# ============================================================
# Figure 4: Stem cell depletion in AS (using pre-computed results)
# ============================================================
print("\n--- Figure 4: AS stem cell depletion ---")

fig = plt.figure(figsize=(14, 10))
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

# Panel A: Cell proportion changes (read Phase 5 results)
ax1 = fig.add_subplot(gs[0, 0])
# Manually entered data (from Phase 5)
prop_data = {
    'CD133+_progenitor': {'WOI': 7.7, 'AS': 0.68, 'padj': 0.011},
    'Secretory_glandular': {'WOI': 16.1, 'AS': 1.8, 'padj': 0.025},
    'NKT_cells': {'WOI': 4.2, 'AS': 19.1, 'padj': 0.035},
    'Macrophages': {'WOI': 0.5, 'AS': 6.8, 'padj': 0.035},
    'SOX9+LGR5+_stem': {'WOI': 6.5, 'AS': 2.4, 'padj': 0.21},
    'Decidualized_stromal': {'WOI': 48.5, 'AS': 34.2, 'padj': 0.15},
}
cts = list(prop_data.keys())
woi_vals = [prop_data[ct]['WOI'] for ct in cts]
as_vals = [prop_data[ct]['AS'] for ct in cts]
x = np.arange(len(cts))
w = 0.35
bars1 = ax1.bar(x - w/2, woi_vals, w, label='WOI Control', color='#2ca02c', alpha=0.8)
bars2 = ax1.bar(x + w/2, as_vals, w, label='AS', color='#d62728', alpha=0.8)
ax1.set_xticks(x)
ax1.set_xticklabels(cts, rotation=45, ha='right', fontsize=6)
ax1.set_ylabel('Proportion (%)')
ax1.set_title('A. Cell Proportion Changes')
ax1.legend(fontsize=7)
# Annotate significance
for i, ct in enumerate(cts):
    p = prop_data[ct]['padj']
    if p < 0.05:
        ax1.text(i, max(woi_vals[i], as_vals[i]) + 1, '*', ha='center', fontsize=12, fontweight='bold')

# Panel B: CytoTRACE changes (from Phase 6 report)
ax2 = fig.add_subplot(gs[0, 1])
cyto_data = {
    'SOX9+LGR5+': {'Ctrl': 0.334, 'AS': 0.072},
    'CD133+': {'Ctrl': 0.230, 'AS': 0.138},
    'Basalis_fib': {'Ctrl': 0.117, 'AS': 0.153},
}
cts_c = list(cyto_data.keys())
ctrl_c = [cyto_data[c]['Ctrl'] for c in cts_c]
as_c = [cyto_data[c]['AS'] for c in cts_c]
x_c = np.arange(len(cts_c))
ax2.bar(x_c - w/2, ctrl_c, w, label='Control', color='#2ca02c', alpha=0.8)
ax2.bar(x_c + w/2, as_c, w, label='AS', color='#d62728', alpha=0.8)
ax2.set_xticks(x_c)
ax2.set_xticklabels(cts_c, fontsize=7)
ax2.set_ylabel('CytoTRACE Score')
ax2.set_title('B. Stemness Score Changes')
ax2.legend(fontsize=7)
ax2.text(0, 0.35, '***', ha='center', fontsize=10)

# Panel C: DEG overview
ax3 = fig.add_subplot(gs[0, 2])
deg_counts = {
    'Decidualized_stromal': (1698, 2177),
    'Lymphatic_endo': (766, 767),
    'Smooth_muscle': (476, 253),
    'Macrophages': (271, 271),
    'Glandular_epi': (219, 290),
    'NKT_cells': (103, 228),
    'CD133+_prog': (14, 209),
    'Lymphoid_imm': (83, 124),
    'Secretory_glan': (89, 2),
    'Basalis_fib': (7, 4),
}
cts_d = list(deg_counts.keys())
up = [deg_counts[c][0] for c in cts_d]
down = [-deg_counts[c][1] for c in cts_d]
y_d = range(len(cts_d))
ax3.barh(y_d, up, color='#d62728', alpha=0.8, label='Up in AS')
ax3.barh(y_d, down, color='#2166AC', alpha=0.8, label='Down in AS')
ax3.set_yticks(y_d)
ax3.set_yticklabels(cts_d, fontsize=6)
ax3.set_xlabel('Number of DEGs')
ax3.set_title('C. Differential Expression')
ax3.axvline(0, color='black', linewidth=0.5)
ax3.legend(fontsize=6)

# Panel D-F: Pathway change heatmap
ax4 = fig.add_subplot(gs[1, :])
try:
    pw_diff = pd.read_csv('Phase_output/phase7_communication/stem_pathway_disease_diff.csv')
    pw_diff['padj'] = pd.to_numeric(pw_diff['padj'], errors='coerce')
    pw_diff['diff'] = pd.to_numeric(pw_diff['diff'], errors='coerce')

    pivot = pw_diff.pivot(index='pathway', columns='celltype', values='diff')
    pivot_p = pw_diff.pivot(index='pathway', columns='celltype', values='padj')

    # Annotation
    annot = pivot.round(3).astype(str)
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

    sns.heatmap(pivot.astype(float), cmap='RdBu_r', center=0, annot=annot, fmt='',
                linewidths=0.5, ax=ax4, annot_kws={'fontsize': 7},
                cbar_kws={'label': 'Expression Diff (AS - Control)'})
    ax4.set_title('D. Signaling Pathway Changes in Stem Cells (AS vs Control)')
    ax4.set_ylabel('')
except Exception as e:
    ax4.text(0.5, 0.5, f'Error: {e}', transform=ax4.transAxes)

plt.savefig(FIGDIR + 'Fig4_AS_stem_depletion.png', dpi=300, bbox_inches='tight')
plt.savefig(FIGDIR + 'Fig4_AS_stem_depletion.pdf', bbox_inches='tight')
plt.close()
print("  Fig4_AS_stem_depletion.png/pdf")

# ============================================================
# Figure 6: Drug targets (using Phase 8 results)
# ============================================================
print("\n--- Figure 6: Drug targets ---")

scores_df = pd.read_csv('Phase_output/phase8_drug_targets/full_target_list.csv', index_col=0)
drug_df = pd.read_csv('Phase_output/phase8_drug_targets/drug_gene_interactions.csv')

fig = plt.figure(figsize=(14, 10))
gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3)

# Panel A: Top 25 target scores
ax1 = fig.add_subplot(gs[0, 0])
top25 = scores_df.head(25)
tier_colors = {'Tier1_approved': '#2ca02c', 'Tier2_investigational': '#ff7f0e',
               'Tier3_novel': '#d62728', 'unknown': '#7f7f7f'}
bar_colors = [tier_colors.get(t, '#7f7f7f') for t in top25['druggability_tier']]
y = range(len(top25))
ax1.barh(y, top25['score'], color=bar_colors, alpha=0.85, edgecolor='black', linewidth=0.3)
ax1.set_yticks(y)
ax1.set_yticklabels(top25.index, fontsize=6)
ax1.invert_yaxis()
ax1.set_xlabel('Priority Score')
ax1.set_title('A. Top 25 Drug Targets')
legend_elements = [Patch(facecolor='#2ca02c', label='Approved drugs'),
                   Patch(facecolor='#ff7f0e', label='Investigational'),
                   Patch(facecolor='#d62728', label='Novel')]
ax1.legend(handles=legend_elements, fontsize=6, loc='lower right')

# Panel B: Evidence type heatmap
ax2 = fig.add_subplot(gs[0, 1])
# Reconstruct evidence matrix
comm_targets = json.load(open('Phase_output/phase7_communication/drug_target_candidates_from_comm.json'))
deg_files = [f for f in os.listdir('Phase_output/phase5_disease/')
             if f.startswith('deg_') and f.endswith('.csv')]
import os
all_degs_genes = {}
for f in deg_files:
    ct = f.replace('deg_', '').replace('_AS_vs_WOI.csv', '')
    df = pd.read_csv(f'Phase_output/phase5_disease/{f}', index_col=0)
    sig = df[(df['padj'] < 0.05) & (df['log2FoldChange'].abs() > 0.25)]
    all_degs_genes[ct] = set(sig.index)

ev_matrix = []
for gene in top25.head(15).index:
    ev = {
        'gene': gene,
        'DEG celltypes': sum(1 for ct, genes in all_degs_genes.items() if gene in genes),
        'Has drug': 1 if top25.loc[gene, 'n_approved_drugs'] > 0 else 0,
        'Comm ligand/receptor': 1 if gene in comm_targets.get('top_changed_ligands', []) or
                                     gene in comm_targets.get('top_changed_receptors', []) else 0,
        'Priority score': int(top25.loc[gene, 'score']),
    }
    ev_matrix.append(ev)

ev_df = pd.DataFrame(ev_matrix).set_index('gene')
sns.heatmap(ev_df, cmap='YlGnBu', annot=True, fmt='.0f', linewidths=0.3,
            ax=ax2, cbar_kws={'label': 'Score'})
ax2.set_title('B. Multi-Evidence Support')
ax2.set_ylabel('')

# Panel C: Druggability pie chart
ax3 = fig.add_subplot(gs[1, 0])
tier_counts = top25['druggability_tier'].value_counts()
tier_labels = {'Tier1_approved': 'Tier1: Approved', 'Tier2_investigational': 'Tier2: Investigational',
               'Tier3_novel': 'Tier3: Novel'}
ax3.pie(tier_counts.values,
        labels=[f"{tier_labels.get(t,t)}\n(n={c})" for t, c in tier_counts.items()],
        colors=[tier_colors.get(t, '#888') for t in tier_counts.index],
        startangle=90, textprops={'fontsize': 8},
        autopct='%1.0f%%')
ax3.set_title('C. Druggability Assessment (Top 25)')

# Panel D: Top drug-target associations
ax4 = fig.add_subplot(gs[1, 1])
# Number of approved drugs for top targets
top_drug_counts = []
for gene in top25.index:
    n = len(drug_df[(drug_df['gene'] == gene) & (drug_df['approved'] == True)]['drug'].unique())
    if n > 0:
        top_drug_counts.append({'gene': gene, 'n_approved_drugs': n})

if top_drug_counts:
    tdc = pd.DataFrame(top_drug_counts).sort_values('n_approved_drugs', ascending=True).tail(15)
    ax4.barh(range(len(tdc)), tdc['n_approved_drugs'].values,
             color='#2ca02c', alpha=0.8, edgecolor='black', linewidth=0.3)
    ax4.set_yticks(range(len(tdc)))
    ax4.set_yticklabels(tdc['gene'].values, fontsize=7)
    ax4.set_xlabel('Number of Approved Drugs')
    ax4.set_title('D. Approved Drug Availability')

plt.savefig(FIGDIR + 'Fig6_drug_targets.png', dpi=300, bbox_inches='tight')
plt.savefig(FIGDIR + 'Fig6_drug_targets.pdf', bbox_inches='tight')
plt.close()
print("  Fig6_drug_targets.png/pdf")

print("\n" + "=" * 60)
print("Publication-quality figure generation complete!")
print("=" * 60)
