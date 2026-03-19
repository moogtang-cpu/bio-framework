#!/usr/bin/env python
"""
Phase 6: Stem cell fate trajectory and regulatory network
- Step 1: Diffusion pseudotime (DPT) + PAGA trajectory
- Step 2: CytoTRACE stemness scoring (CellRank)
- Step 3: Normal vs disease trajectory comparison
- Step 4: TF activity inference (decoupler + collectri)
- Step 5: Key transcription factor identification

Note: RNA velocity requires spliced/unspliced layers, which are unavailable in the current data - skipped.
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import cellrank as cr
import decoupler as dc
from scipy import stats
from scipy.sparse import issparse
import json
import os
import warnings
warnings.filterwarnings('ignore')

FIGDIR = 'Phase_output/phase6_trajectory/figures/'
OUTDIR = 'Phase_output/phase6_trajectory/'
os.makedirs(FIGDIR, exist_ok=True)

print("=" * 60)
print("Phase 6: Stem Cell Fate Trajectory and Regulatory Network")
print("=" * 60)

# ============================================================
# 1. Load data, extract epithelial + stem cell lineage
# ============================================================
print("\n--- Loading Data ---")
adata = ad.read_h5ad('Phase_output/phase3_integration/integrated_harmony.h5ad')
print(f"Integrated atlas: {adata.shape}")

# Define cell type subsets for trajectory analysis
# Epithelial/stem cell differentiation trajectory: SOX9+LGR5+_stem -> CD133+_progenitor -> Glandular/Secretory
# Mesenchymal trajectory: Basalis_fibroblast -> Decidualized_stromal
epithelial_types = ['SOX9+LGR5+_stem', 'CD133+_progenitor',
                    'Glandular_epithelium', 'Secretory_glandular']
mesenchymal_types = ['SOX9+LGR5+_stem', 'Basalis_fibroblast',
                     'Decidualized_stromal', 'Smooth_muscle']

# Main trajectory: epithelial differentiation (including stem cells)
epi_mask = adata.obs['celltype_manual'].isin(epithelial_types)
adata_epi = adata[epi_mask].copy()
print(f"Epithelial/stem cell lineage: {adata_epi.shape[0]} cells")
for ct in epithelial_types:
    n = (adata_epi.obs['celltype_manual'] == ct).sum()
    print(f"  {ct}: {n}")

# Mesenchymal trajectory
mes_mask = adata.obs['celltype_manual'].isin(mesenchymal_types)
adata_mes = adata[mes_mask].copy()
print(f"Mesenchymal lineage: {adata_mes.shape[0]} cells")

# ============================================================
# 2. Epithelial lineage DPT trajectory analysis
# ============================================================
print("\n--- Epithelial Lineage DPT Trajectory ---")

# Use Harmony-corrected PCA
# Recompute neighbors and UMAP (epithelial subset only)
sc.pp.neighbors(adata_epi, use_rep='X_pca_harmony', n_neighbors=30, random_state=42)
sc.tl.umap(adata_epi, random_state=42)

# Leiden clustering (fine resolution)
sc.tl.leiden(adata_epi, resolution=0.8, key_added='leiden_epi', random_state=42)
n_clusters = adata_epi.obs['leiden_epi'].nunique()
print(f"  Epithelial subset Leiden clustering: {n_clusters} clusters")

# PAGA trajectory inference
sc.tl.paga(adata_epi, groups='celltype_manual')
print("  PAGA computation complete")

# Visualize PAGA
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
sc.pl.paga(adata_epi, ax=axes[0], show=False, fontsize=10, frameon=False,
           title='PAGA - Epithelial/Stem Trajectory')
sc.pl.umap(adata_epi, color='celltype_manual', ax=axes[1], show=False,
           title='UMAP - Epithelial/Stem Lineage', frameon=False)
plt.tight_layout()
plt.savefig(FIGDIR + 'paga_epithelial_trajectory.png', dpi=200, bbox_inches='tight')
plt.close()
print("  paga_epithelial_trajectory.png")

# Diffusion Pseudotime
sc.tl.diffmap(adata_epi, n_comps=15)

# Select root cell (cell with highest stemness marker expression in SOX9+LGR5+_stem)
stem_mask = adata_epi.obs['celltype_manual'] == 'SOX9+LGR5+_stem'
stem_indices = np.where(stem_mask)[0]

# Use SOX9 expression to select root
if 'SOX9' in adata_epi.var_names:
    X = adata_epi.layers.get('log_normalized', adata_epi.X)
    if issparse(X):
        sox9_expr = np.asarray(X[stem_indices, adata_epi.var_names.get_loc('SOX9')]).flatten()
    else:
        sox9_expr = X[stem_indices, adata_epi.var_names.get_loc('SOX9')]
    root_idx = stem_indices[np.argmax(sox9_expr)]
else:
    # Use the first component of diffmap
    root_idx = stem_indices[np.argmin(adata_epi.obsm['X_diffmap'][stem_indices, 0])]

adata_epi.uns['iroot'] = root_idx
sc.tl.dpt(adata_epi)
print(f"  DPT complete (root: index {root_idx}, celltype={adata_epi.obs['celltype_manual'].iloc[root_idx]})")

# DPT visualization
fig, axes = plt.subplots(1, 3, figsize=(21, 6))
sc.pl.umap(adata_epi, color='dpt_pseudotime', ax=axes[0], show=False,
           title='Diffusion Pseudotime', cmap='viridis', frameon=False)
sc.pl.umap(adata_epi, color='celltype_manual', ax=axes[1], show=False,
           title='Cell Type', frameon=False)

# DPT distribution by cell type
ct_order = ['SOX9+LGR5+_stem', 'CD133+_progenitor', 'Glandular_epithelium', 'Secretory_glandular']
dpt_data = []
for ct in ct_order:
    vals = adata_epi.obs.loc[adata_epi.obs['celltype_manual'] == ct, 'dpt_pseudotime'].dropna()
    dpt_data.append(vals)
bp = axes[2].boxplot(dpt_data, labels=[c.replace('+', '+\n') for c in ct_order],
                     patch_artist=True)
colors = ['#E64B35', '#4DBBD5', '#00A087', '#3C5488']
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)
axes[2].set_ylabel('Pseudotime')
axes[2].set_title('Pseudotime Distribution by Cell Type')
plt.tight_layout()
plt.savefig(FIGDIR + 'dpt_epithelial.png', dpi=200, bbox_inches='tight')
plt.close()
print("  dpt_epithelial.png")

# DPT statistics: median pseudotime per cell type
print("\n  Pseudotime median (epithelial lineage):")
for ct in ct_order:
    vals = adata_epi.obs.loc[adata_epi.obs['celltype_manual'] == ct, 'dpt_pseudotime'].dropna()
    print(f"    {ct}: median={vals.median():.3f}, mean={vals.mean():.3f}")

# ============================================================
# 3. CytoTRACE stemness scoring
# ============================================================
print("\n--- CytoTRACE Stemness Scoring ---")

# CellRank CytoTRACE kernel
# Should be run on the full dataset
# Use simplified CytoTRACE score: based on number of detected genes (negatively correlated with stemness)
# Standard CytoTRACE: more expressed genes -> higher stemness
cytotrace_types = epithelial_types + ['Basalis_fibroblast', 'Decidualized_stromal']
adata_ct = adata[adata.obs['celltype_manual'].isin(cytotrace_types)].copy()

# Compute number of detected genes per cell (already stored in n_genes_by_counts)
if 'n_genes_by_counts' not in adata_ct.obs.columns:
    X = adata_ct.layers.get('counts', adata_ct.X)
    if issparse(X):
        adata_ct.obs['n_genes_by_counts'] = np.asarray((X > 0).sum(axis=1)).flatten()
    else:
        adata_ct.obs['n_genes_by_counts'] = (X > 0).sum(axis=1)

# Simplified CytoTRACE score: normalized gene count
gene_counts = adata_ct.obs['n_genes_by_counts'].values.astype(float)
adata_ct.obs['cytotrace_score'] = (gene_counts - gene_counts.min()) / (gene_counts.max() - gene_counts.min())

# Attempt CellRank CytoTRACE kernel
try:
    sc.pp.neighbors(adata_ct, use_rep='X_pca_harmony', n_neighbors=30, random_state=42)
    sc.tl.umap(adata_ct, random_state=42)

    ctk = cr.kernels.CytoTRACEKernel(adata_ct)
    ctk.compute_cytotrace()
    print("  CellRank CytoTRACE kernel complete")
    has_cellrank_ct = True
except Exception as e:
    print(f"  CellRank CytoTRACE: {e}")
    print("  Using simplified CytoTRACE score")
    has_cellrank_ct = False
    sc.pp.neighbors(adata_ct, use_rep='X_pca_harmony', n_neighbors=30, random_state=42)
    sc.tl.umap(adata_ct, random_state=42)

# CytoTRACE visualization
fig, axes = plt.subplots(1, 3, figsize=(21, 6))

score_col = 'ct_pseudotime' if has_cellrank_ct and 'ct_pseudotime' in adata_ct.obs.columns else 'cytotrace_score'
sc.pl.umap(adata_ct, color=score_col, ax=axes[0], show=False,
           title='CytoTRACE Stemness Score', cmap='RdYlGn_r', frameon=False)
sc.pl.umap(adata_ct, color='celltype_manual', ax=axes[1], show=False,
           title='Cell Type', frameon=False)

# CytoTRACE by cell type
ct_labels = ['SOX9+LGR5+_stem', 'CD133+_progenitor', 'Basalis_fibroblast',
             'Glandular_epithelium', 'Secretory_glandular', 'Decidualized_stromal']
ct_scores = []
for ct in ct_labels:
    vals = adata_ct.obs.loc[adata_ct.obs['celltype_manual'] == ct, score_col].dropna()
    ct_scores.append(vals)

bp = axes[2].boxplot(ct_scores, labels=[c[:12] for c in ct_labels], patch_artist=True)
colors = ['#E64B35', '#4DBBD5', '#F39B7F', '#00A087', '#3C5488', '#8491B4']
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)
axes[2].set_ylabel('CytoTRACE Score')
axes[2].set_title('Stemness Score by Cell Type')
axes[2].tick_params(axis='x', rotation=45)
plt.tight_layout()
plt.savefig(FIGDIR + 'cytotrace_stemness.png', dpi=200, bbox_inches='tight')
plt.close()
print("  cytotrace_stemness.png")

# CytoTRACE statistics
print("\n  CytoTRACE stemness scores:")
for ct in ct_labels:
    vals = adata_ct.obs.loc[adata_ct.obs['celltype_manual'] == ct, score_col].dropna()
    print(f"    {ct}: median={vals.median():.3f}, mean={vals.mean():.3f}")

# ============================================================
# 4. Normal vs disease trajectory comparison
# ============================================================
print("\n--- Normal vs Disease Trajectory Comparison ---")

# Compare pseudotime distributions between Normal/WOI and AS in epithelial lineage
# Keep only cells with condition labels
epi_disease = adata_epi[adata_epi.obs['condition'].isin(['Normal', 'WOI Control', 'AS'])].copy()

# Merge Normal and WOI Control into "Control"
epi_disease.obs['disease_group'] = epi_disease.obs['condition'].map({
    'Normal': 'Control', 'WOI Control': 'Control', 'AS': 'AS'
})

fig, axes = plt.subplots(2, 4, figsize=(24, 10))

for idx, ct in enumerate(ct_order):
    ax = axes[idx // 4, idx % 4]
    ct_data = epi_disease[epi_disease.obs['celltype_manual'] == ct]

    for grp, color, label in [('Control', '#4DBBD5', 'Control'),
                               ('AS', '#E64B35', 'AS')]:
        vals = ct_data.obs.loc[ct_data.obs['disease_group'] == grp, 'dpt_pseudotime'].dropna()
        if len(vals) > 0:
            ax.hist(vals, bins=30, density=True, alpha=0.5, color=color, label=f'{label} (n={len(vals)})')

    ax.set_xlabel('Pseudotime')
    ax.set_ylabel('Density')
    ax.set_title(ct)
    ax.legend(fontsize=8)

    # KS test
    ctrl_vals = ct_data.obs.loc[ct_data.obs['disease_group'] == 'Control', 'dpt_pseudotime'].dropna()
    as_vals = ct_data.obs.loc[ct_data.obs['disease_group'] == 'AS', 'dpt_pseudotime'].dropna()
    if len(ctrl_vals) > 10 and len(as_vals) > 10:
        ks_stat, ks_pval = stats.ks_2samp(ctrl_vals, as_vals)
        ax.text(0.05, 0.95, f'KS p={ks_pval:.2e}', transform=ax.transAxes,
                fontsize=8, va='top')
        print(f"  {ct}: Control median={ctrl_vals.median():.3f}, AS median={as_vals.median():.3f}, KS p={ks_pval:.2e}")

# Hide unused subplots
for idx in range(len(ct_order), 8):
    axes[idx // 4, idx % 4].set_visible(False)

fig.suptitle('Pseudotime Distribution: Control vs AS (Epithelial Lineage)', fontsize=14)
plt.tight_layout()
plt.savefig(FIGDIR + 'trajectory_disease_comparison.png', dpi=200, bbox_inches='tight')
plt.close()
print("  trajectory_disease_comparison.png")

# CytoTRACE disease comparison
ct_disease = adata_ct[adata_ct.obs['condition'].isin(['Normal', 'WOI Control', 'AS'])].copy()
ct_disease.obs['disease_group'] = ct_disease.obs['condition'].map({
    'Normal': 'Control', 'WOI Control': 'Control', 'AS': 'AS'
})

print("\n  CytoTRACE disease differences:")
ct_disease_results = []
for ct in ct_labels:
    ct_data = ct_disease[ct_disease.obs['celltype_manual'] == ct]
    ctrl_vals = ct_data.obs.loc[ct_data.obs['disease_group'] == 'Control', score_col].dropna()
    as_vals = ct_data.obs.loc[ct_data.obs['disease_group'] == 'AS', score_col].dropna()
    if len(ctrl_vals) > 10 and len(as_vals) > 10:
        stat, pval = stats.mannwhitneyu(ctrl_vals, as_vals, alternative='two-sided')
        diff = as_vals.median() - ctrl_vals.median()
        print(f"    {ct}: Ctrl={ctrl_vals.median():.3f}, AS={as_vals.median():.3f}, diff={diff:+.3f}, p={pval:.2e}")
        ct_disease_results.append({
            'celltype': ct, 'ctrl_median': float(ctrl_vals.median()),
            'as_median': float(as_vals.median()), 'diff': float(diff), 'pvalue': float(pval)
        })

# ============================================================
# 5. TF activity inference (decoupler + CollecTRI)
# ============================================================
print("\n--- TF Activity Inference ---")

# Use decoupler CollecTRI/DoRothEA regulon database
# Infer TF activity on the full atlas

# Retrieve CollecTRI network
try:
    net = dc.get_collectri(organism='human', split_complexes=False)
    print(f"  CollecTRI network: {len(net)} regulatory relationships, {net['source'].nunique()} TFs")
except Exception as e:
    print(f"  CollecTRI retrieval error: {e}, trying DoRothEA")
    try:
        net = dc.get_dorothea(organism='human')
        net = net[net['confidence'].isin(['A', 'B', 'C'])]
        print(f"  DoRothEA network (ABC): {len(net)} regulatory relationships")
    except Exception as e2:
        print(f"  DoRothEA also failed: {e2}")
        net = None

if net is not None:
    # Infer TF activity on stem cell + epithelial subset
    adata_tf = adata[adata.obs['celltype_manual'].isin(cytotrace_types)].copy()

    # Use log_normalized layer
    if 'log_normalized' in adata_tf.layers:
        adata_tf.X = adata_tf.layers['log_normalized'].copy()

    # Run ULM (Univariate Linear Model) to infer TF activity
    print("  Running ULM TF activity inference...")
    dc.run_ulm(
        mat=adata_tf,
        net=net,
        source='source',
        target='target',
        weight='weight',
        verbose=True,
        use_raw=False,
    )

    # TF activities stored in adata_tf.obsm['ulm_estimate']
    tf_activities = adata_tf.obsm['ulm_estimate']
    print(f"  Inferred activity for {tf_activities.shape[1]} TFs")

    # 5a. Top TFs per cell type
    print("\n  Top active TFs per cell type:")
    tf_top = {}
    for ct in cytotrace_types:
        ct_mask = adata_tf.obs['celltype_manual'] == ct
        mean_act = tf_activities[ct_mask].mean(axis=0)
        top_tfs = mean_act.nlargest(10)
        tf_top[ct] = top_tfs.index.tolist()
        print(f"    {ct}: {', '.join(top_tfs.index[:5])}")

    # 5b. TF activity heatmap (top TFs across cell types)
    mean_tf = pd.DataFrame(index=cytotrace_types)
    for ct in cytotrace_types:
        ct_mask = adata_tf.obs['celltype_manual'] == ct
        mean_tf.loc[ct] = tf_activities[ct_mask].mean(axis=0)

    # Select top variable TFs
    tf_var = mean_tf.var(axis=0).nlargest(30).index
    plot_tf = mean_tf[tf_var]

    fig, ax = plt.subplots(figsize=(14, 6))
    sns.heatmap(plot_tf, cmap='RdBu_r', center=0, annot=True, fmt='.2f',
                linewidths=0.3, ax=ax, cbar_kws={'label': 'TF Activity (ULM)'})
    ax.set_title('Top Variable TF Activities across Cell Types')
    ax.set_ylabel('')
    plt.tight_layout()
    plt.savefig(FIGDIR + 'tf_activity_heatmap.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("  tf_activity_heatmap.png")

    # 5c. Disease-differential TF activity
    print("\n  Disease-differential TF activity (Control vs AS):")
    adata_tf.obs['disease_group'] = adata_tf.obs['condition'].map({
        'Normal': 'Control', 'WOI Control': 'Control', 'AS': 'AS',
        'AS_CD133': 'AS_CD133', 'endometrium': 'Other', 'decidua': 'Other'
    })

    disease_tf_results = {}
    for ct in ['SOX9+LGR5+_stem', 'CD133+_progenitor', 'Basalis_fibroblast',
               'Decidualized_stromal', 'Glandular_epithelium']:
        ct_data = adata_tf[
            (adata_tf.obs['celltype_manual'] == ct) &
            (adata_tf.obs['disease_group'].isin(['Control', 'AS']))
        ]

        if ct_data.shape[0] < 50:
            continue

        ctrl_mask = ct_data.obs['disease_group'] == 'Control'
        as_mask = ct_data.obs['disease_group'] == 'AS'

        ctrl_act = ct_data.obsm['ulm_estimate'][ctrl_mask.values]
        as_act = ct_data.obsm['ulm_estimate'][as_mask.values]

        if ctrl_act.shape[0] < 20 or as_act.shape[0] < 20:
            continue

        # Wilcoxon test per TF
        tf_diff = []
        for tf in tf_activities.columns:
            c_vals = ctrl_act[tf].values
            a_vals = as_act[tf].values
            if np.std(c_vals) > 0 or np.std(a_vals) > 0:
                try:
                    stat, pval = stats.mannwhitneyu(c_vals, a_vals, alternative='two-sided')
                except:
                    pval = 1.0
                diff = a_vals.mean() - c_vals.mean()
                tf_diff.append({'TF': tf, 'diff_AS_vs_Control': diff, 'pvalue': pval})

        tf_diff_df = pd.DataFrame(tf_diff)
        if len(tf_diff_df) == 0:
            continue

        # BH correction
        from statsmodels.stats.multitest import multipletests
        _, padj, _, _ = multipletests(tf_diff_df['pvalue'], method='fdr_bh')
        tf_diff_df['padj'] = padj
        tf_diff_df = tf_diff_df.sort_values('padj')

        sig_tfs = tf_diff_df[tf_diff_df['padj'] < 0.05]
        n_up = (sig_tfs['diff_AS_vs_Control'] > 0).sum()
        n_down = (sig_tfs['diff_AS_vs_Control'] < 0).sum()
        print(f"    {ct}: {n_up}↑ + {n_down}↓ TFs (padj<0.05)")
        if len(sig_tfs) > 0:
            print(f"      Top up: {', '.join(sig_tfs[sig_tfs['diff_AS_vs_Control'] > 0].head(5)['TF'].tolist())}")
            print(f"      Top down: {', '.join(sig_tfs[sig_tfs['diff_AS_vs_Control'] < 0].head(5)['TF'].tolist())}")

        disease_tf_results[ct] = tf_diff_df
        tf_diff_df.to_csv(f'{OUTDIR}tf_diff_{ct.replace("+","plus").replace("/","_")}_AS_vs_Control.csv', index=False)

    # 5d. Key TF network for stem cells
    stem_key_tfs = set()
    for ct in ['SOX9+LGR5+_stem', 'CD133+_progenitor', 'Basalis_fibroblast']:
        if ct in disease_tf_results:
            top10 = disease_tf_results[ct].head(10)['TF'].tolist()
            stem_key_tfs.update(top10)

    if stem_key_tfs:
        # TF activity comparison (Control vs AS) for stem cell top TFs
        stem_tfs = sorted(stem_key_tfs)[:20]  # up to 20
        fig, ax = plt.subplots(figsize=(12, 6))

        plot_data = []
        for ct in ['SOX9+LGR5+_stem', 'CD133+_progenitor', 'Basalis_fibroblast']:
            if ct in disease_tf_results:
                df = disease_tf_results[ct]
                for tf in stem_tfs:
                    match = df[df['TF'] == tf]
                    if len(match) > 0:
                        plot_data.append({
                            'TF': tf, 'celltype': ct,
                            'diff': match.iloc[0]['diff_AS_vs_Control'],
                            'padj': match.iloc[0]['padj'],
                        })

        if plot_data:
            plot_df = pd.DataFrame(plot_data)
            pivot = plot_df.pivot(index='TF', columns='celltype', values='diff')
            pivot = pivot.fillna(0)

            fig, ax = plt.subplots(figsize=(10, max(5, len(pivot) * 0.4)))
            sns.heatmap(pivot, cmap='RdBu_r', center=0, annot=True, fmt='.2f',
                        linewidths=0.3, ax=ax, cbar_kws={'label': 'TF Activity Diff (AS - Control)'})
            ax.set_title('Key TF Activity Changes in Stem/Progenitor Cells (AS vs Control)')
            plt.tight_layout()
            plt.savefig(FIGDIR + 'tf_stem_disease_diff.png', dpi=200, bbox_inches='tight')
            plt.close()
            print("  tf_stem_disease_diff.png")

# ============================================================
# 6. Mesenchymal trajectory analysis (supplementary)
# ============================================================
print("\n--- Mesenchymal Lineage Trajectory ---")

sc.pp.neighbors(adata_mes, use_rep='X_pca_harmony', n_neighbors=30, random_state=42)
sc.tl.umap(adata_mes, random_state=42)
sc.tl.paga(adata_mes, groups='celltype_manual')
sc.tl.diffmap(adata_mes, n_comps=15)

# Root: Basalis_fibroblast (stem cell niche component)
bf_mask = adata_mes.obs['celltype_manual'] == 'Basalis_fibroblast'
bf_indices = np.where(bf_mask)[0]
if len(bf_indices) > 0:
    adata_mes.uns['iroot'] = bf_indices[0]
    sc.tl.dpt(adata_mes)

    fig, axes = plt.subplots(1, 3, figsize=(21, 6))
    sc.pl.umap(adata_mes, color='dpt_pseudotime', ax=axes[0], show=False,
               title='Pseudotime (Mesenchymal)', cmap='viridis', frameon=False)
    sc.pl.umap(adata_mes, color='celltype_manual', ax=axes[1], show=False,
               title='Cell Type', frameon=False)
    sc.pl.paga(adata_mes, ax=axes[2], show=False, fontsize=10, frameon=False,
               title='PAGA - Mesenchymal')
    plt.tight_layout()
    plt.savefig(FIGDIR + 'dpt_mesenchymal.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("  dpt_mesenchymal.png")

    # Mesenchymal pseudotime statistics
    print("\n  Pseudotime median (mesenchymal lineage):")
    for ct in mesenchymal_types:
        vals = adata_mes.obs.loc[adata_mes.obs['celltype_manual'] == ct, 'dpt_pseudotime'].dropna()
        if len(vals) > 0:
            print(f"    {ct}: median={vals.median():.3f}")

# ============================================================
# 7. Key gene expression changes along trajectory
# ============================================================
print("\n--- Key Gene Expression along Trajectory ---")

trajectory_genes = {
    'Stem': ['SOX9', 'LGR5', 'PROM1', 'AXIN2', 'ALDH1A1'],
    'Differentiation': ['MUC1', 'FOXA2', 'PAX8', 'EPCAM', 'PRL'],
    'Wnt': ['WNT4', 'RSPO1', 'SFRP1', 'DKK1', 'CTNNB1'],
    'Niche': ['PDGFRB', 'COL1A1', 'DCN', 'VIM', 'ACTA2'],
}

available = set(adata_epi.var_names)
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

for idx, (category, genes) in enumerate(trajectory_genes.items()):
    ax = axes[idx // 2, idx % 2]
    present_genes = [g for g in genes if g in available]

    if not present_genes:
        ax.text(0.5, 0.5, f'No genes found for {category}', ha='center', va='center')
        continue

    # Sort by pseudotime
    sorted_idx = adata_epi.obs['dpt_pseudotime'].dropna().sort_values().index
    adata_sorted = adata_epi[sorted_idx]

    X = adata_sorted.layers.get('log_normalized', adata_sorted.X)
    if issparse(X):
        X = X.toarray()

    # Sliding window smoothing
    window = 200
    pt = adata_sorted.obs['dpt_pseudotime'].values

    for gene in present_genes:
        gene_idx = list(adata_sorted.var_names).index(gene)
        expr = X[:, gene_idx]
        # Rolling average
        smoothed = pd.Series(expr).rolling(window=window, min_periods=50, center=True).mean()
        ax.plot(pt, smoothed, label=gene, linewidth=1.5, alpha=0.8)

    ax.set_xlabel('Pseudotime')
    ax.set_ylabel('Expression (smoothed)')
    ax.set_title(f'{category} Genes along Trajectory')
    ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig(FIGDIR + 'trajectory_gene_expression.png', dpi=200, bbox_inches='tight')
plt.close()
print("  trajectory_gene_expression.png")

# ============================================================
# 8. Save results
# ============================================================
print("\n--- Saving Results ---")

# Save trajectory data
adata_epi.write_h5ad(OUTDIR + 'epithelial_trajectory.h5ad')
print(f"  epithelial_trajectory.h5ad ({adata_epi.shape})")

adata_mes.write_h5ad(OUTDIR + 'mesenchymal_trajectory.h5ad')
print(f"  mesenchymal_trajectory.h5ad ({adata_mes.shape})")

# Save report
phase6_report = {
    'phase': 6,
    'epithelial_trajectory': {
        'n_cells': int(adata_epi.shape[0]),
        'cell_types': epithelial_types,
        'pseudotime_medians': {
            ct: float(adata_epi.obs.loc[adata_epi.obs['celltype_manual'] == ct, 'dpt_pseudotime'].median())
            for ct in ct_order
        },
    },
    'mesenchymal_trajectory': {
        'n_cells': int(adata_mes.shape[0]),
        'cell_types': mesenchymal_types,
    },
    'cytotrace': {
        'method': 'CellRank CytoTRACE kernel' if has_cellrank_ct else 'simplified gene count',
        'score_column': score_col,
    },
    'tf_analysis': {
        'method': 'decoupler ULM + CollecTRI',
        'n_tfs': int(tf_activities.shape[1]) if net is not None else 0,
        'disease_diff_celltypes': list(disease_tf_results.keys()) if net is not None else [],
    },
    'disease_trajectory_comparison': ct_disease_results,
}

with open(OUTDIR + 'phase6_report.json', 'w') as f:
    json.dump(phase6_report, f, indent=2, default=str)
print("  phase6_report.json")

print("\n" + "=" * 60)
print("Phase 6 Complete!")
print(f"  Epithelial trajectory: {adata_epi.shape[0]} cells, DPT complete")
print(f"  Mesenchymal trajectory: {adata_mes.shape[0]} cells")
if net is not None:
    print(f"  TF activity: {tf_activities.shape[1]} TFs inferred")
print("=" * 60)
