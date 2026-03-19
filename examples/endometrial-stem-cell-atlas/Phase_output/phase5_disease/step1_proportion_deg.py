#!/usr/bin/env python
"""
Phase 5 Step 1: Disease comparison analysis — cell proportion changes + Pseudobulk DEG
Primary comparison: WOI Control vs AS (GSE215968, same study to reduce batch effects)
Secondary: RIF vs CTR (Phase 4 Visium deconvolution)
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.sparse import issparse
import json
import os
import warnings
warnings.filterwarnings('ignore')

FIGDIR = 'Phase_output/phase5_disease/figures/'
OUTDIR = 'Phase_output/phase5_disease/'
os.makedirs(FIGDIR, exist_ok=True)

# ============================================================
# 1. Load integrated atlas
# ============================================================
print("=" * 60)
print("Phase 5 Step 1: Disease differential analysis")
print("=" * 60)

print("\n--- Loading integrated atlas ---")
adata = ad.read_h5ad('Phase_output/phase3_integration/integrated_harmony.h5ad')
print(f"Integrated atlas: {adata.shape[0]} cells x {adata.shape[1]} genes")

# ============================================================
# 2. Define disease groups
# ============================================================
print("\n--- Defining disease groups ---")

# Primary analysis: WOI Control vs AS (same study GSE215968)
# Excluded: AS_CD133 (CD133+ sorted, proportions not comparable), endometrium/decidua (different tissues)
adata_disease = adata[adata.obs['condition'].isin(['WOI Control', 'AS'])].copy()
print(f"Primary analysis (WOI Control vs AS): {adata_disease.shape[0]} cells")
print(f"  WOI Control: {(adata_disease.obs['condition'] == 'WOI Control').sum()}")
print(f"  AS: {(adata_disease.obs['condition'] == 'AS').sum()}")

# Create unique sample ID (dataset + sample)
adata_disease.obs['sample_id'] = adata_disease.obs['dataset'].astype(str) + '_' + adata_disease.obs['sample'].astype(str)

# Filter very small samples (< 50 cells)
sample_counts = adata_disease.obs['sample_id'].value_counts()
valid_samples = sample_counts[sample_counts >= 50].index
n_removed = (adata_disease.obs['sample_id'].isin(sample_counts[sample_counts < 50].index)).sum()
adata_disease = adata_disease[adata_disease.obs['sample_id'].isin(valid_samples)].copy()
print(f"  After filtering samples with <50 cells: {adata_disease.shape[0]} cells, {n_removed} cells removed")

# Sample distribution
sample_cond = adata_disease.obs.groupby('sample_id')['condition'].first()
n_woi = (sample_cond == 'WOI Control').sum()
n_as = (sample_cond == 'AS').sum()
print(f"  WOI Control samples: {n_woi}, AS samples: {n_as}")

# ============================================================
# 3. Cell proportion change analysis
# ============================================================
print("\n--- Cell proportion analysis ---")

celltype_col = 'celltype_manual'
celltypes = sorted(adata_disease.obs[celltype_col].unique())

# 3a. Cell proportions per sample
prop_data = []
for sid in adata_disease.obs['sample_id'].unique():
    mask = adata_disease.obs['sample_id'] == sid
    cond = adata_disease.obs.loc[mask, 'condition'].iloc[0]
    total = mask.sum()
    for ct in celltypes:
        ct_count = (adata_disease.obs.loc[mask, celltype_col] == ct).sum()
        prop_data.append({
            'sample_id': sid,
            'condition': cond,
            'celltype': ct,
            'count': ct_count,
            'proportion': ct_count / total
        })
prop_df = pd.DataFrame(prop_data)

# 3b. Statistical testing (Wilcoxon rank-sum test per cell type)
prop_results = []
for ct in celltypes:
    ct_df = prop_df[prop_df['celltype'] == ct]
    woi_props = ct_df[ct_df['condition'] == 'WOI Control']['proportion'].values
    as_props = ct_df[ct_df['condition'] == 'AS']['proportion'].values

    # Compute means and fold change
    mean_woi = np.mean(woi_props)
    mean_as = np.mean(as_props)
    log2fc = np.log2((mean_as + 1e-6) / (mean_woi + 1e-6))

    # Wilcoxon rank-sum test
    if len(woi_props) >= 3 and len(as_props) >= 3:
        stat, pval = stats.mannwhitneyu(woi_props, as_props, alternative='two-sided')
    else:
        pval = np.nan

    prop_results.append({
        'celltype': ct,
        'mean_WOI': mean_woi,
        'mean_AS': mean_as,
        'log2FC_AS_vs_WOI': log2fc,
        'pvalue': pval,
    })

prop_results_df = pd.DataFrame(prop_results)

# BH correction
from statsmodels.stats.multitest import multipletests
valid_mask = ~prop_results_df['pvalue'].isna()
if valid_mask.sum() > 0:
    _, padj, _, _ = multipletests(prop_results_df.loc[valid_mask, 'pvalue'], method='fdr_bh')
    prop_results_df.loc[valid_mask, 'padj'] = padj
else:
    prop_results_df['padj'] = np.nan

prop_results_df = prop_results_df.sort_values('pvalue')
print("\nCell proportion changes (AS vs WOI Control):")
print(prop_results_df.to_string(index=False, float_format='%.4f'))

# 3c. Visualization - cell proportion boxplots
fig, axes = plt.subplots(2, 7, figsize=(28, 10))
axes = axes.flatten()
for idx, ct in enumerate(celltypes):
    ax = axes[idx]
    ct_df = prop_df[prop_df['celltype'] == ct]

    # Boxplot
    for i, (cond, color) in enumerate([('WOI Control', '#4DBBD5'), ('AS', '#E64B35')]):
        vals = ct_df[ct_df['condition'] == cond]['proportion'].values
        bp = ax.boxplot([vals], positions=[i], widths=0.5, patch_artist=True,
                       boxprops=dict(facecolor=color, alpha=0.6),
                       medianprops=dict(color='black'))
        ax.scatter([i] * len(vals), vals, c=color, s=20, alpha=0.7, zorder=3)

    # Annotate significance
    row = prop_results_df[prop_results_df['celltype'] == ct].iloc[0]
    if row['padj'] < 0.001:
        sig = '***'
    elif row['padj'] < 0.01:
        sig = '**'
    elif row['padj'] < 0.05:
        sig = '*'
    else:
        sig = 'ns'

    ax.set_title(f'{ct}\n{sig}', fontsize=9)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['WOI', 'AS'], fontsize=8)
    ax.set_ylabel('Proportion', fontsize=8)

# Hide extra subplots
for idx in range(len(celltypes), len(axes)):
    axes[idx].set_visible(False)

fig.suptitle('Cell Type Proportion Changes: WOI Control vs Asherman Syndrome', fontsize=14)
plt.tight_layout()
plt.savefig(FIGDIR + 'proportion_boxplot_WOI_vs_AS.png', dpi=200, bbox_inches='tight')
plt.close()
print("\n  proportion_boxplot_WOI_vs_AS.png")

# 3d. Stacked bar chart
fig, ax = plt.subplots(figsize=(14, 6))
pivot = prop_df.groupby(['condition', 'celltype'])['proportion'].mean().unstack(fill_value=0)
pivot.plot(kind='bar', stacked=True, ax=ax, colormap='tab20', edgecolor='white', linewidth=0.5)
ax.set_ylabel('Mean Proportion')
ax.set_title('Cell Type Composition: WOI Control vs AS')
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.tight_layout()
plt.savefig(FIGDIR + 'proportion_stacked_bar.png', dpi=200, bbox_inches='tight')
plt.close()
print("  proportion_stacked_bar.png")

# 3e. Log2FC forest plot
fig, ax = plt.subplots(figsize=(10, 7))
df_plot = prop_results_df.sort_values('log2FC_AS_vs_WOI')
colors = ['#E64B35' if p < 0.05 else '#999999' for p in df_plot['padj']]
ax.barh(range(len(df_plot)), df_plot['log2FC_AS_vs_WOI'], color=colors, height=0.7)
ax.set_yticks(range(len(df_plot)))
ax.set_yticklabels(df_plot['celltype'], fontsize=9)
ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
ax.set_xlabel('log2(Fold Change) AS vs WOI Control')
ax.set_title('Cell Type Proportion Changes in Asherman Syndrome')
for i, (_, row) in enumerate(df_plot.iterrows()):
    if row['padj'] < 0.05:
        ax.text(row['log2FC_AS_vs_WOI'] + 0.05 * np.sign(row['log2FC_AS_vs_WOI']),
                i, f"p={row['padj']:.2e}", va='center', fontsize=7)
plt.tight_layout()
plt.savefig(FIGDIR + 'proportion_log2fc_forest.png', dpi=200, bbox_inches='tight')
plt.close()
print("  proportion_log2fc_forest.png")

# ============================================================
# 4. Pseudobulk DEG analysis (PyDESeq2)
# ============================================================
print("\n--- Pseudobulk DEG analysis ---")

# Ensure raw counts are used
if 'counts' in adata_disease.layers:
    count_layer = 'counts'
elif issparse(adata_disease.raw.X) if adata_disease.raw is not None else False:
    count_layer = None  # use raw
else:
    count_layer = None

# Pseudobulk aggregation: sum by sample_id and celltype
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# DEG analysis for each key cell type
# Focus on stem cell-related + major cell types
focus_celltypes = [
    'SOX9+LGR5+_stem', 'CD133+_progenitor', 'Basalis_fibroblast',
    'Decidualized_stromal', 'Glandular_epithelium', 'Secretory_glandular',
    'Smooth_muscle', 'NKT_cells', 'Macrophages', 'Lymphoid_immune',
    'Lymphatic_endothelium'
]

all_deg_results = {}
deg_summary = {}

for ct in focus_celltypes:
    print(f"\n  Processing {ct}...")

    # Extract this cell type
    mask_ct = adata_disease.obs[celltype_col] == ct
    adata_ct = adata_disease[mask_ct].copy()

    # Check that both conditions have sufficient samples
    samples_per_cond = adata_ct.obs.groupby('condition')['sample_id'].nunique()
    min_samples = 3
    if any(samples_per_cond < min_samples):
        print(f"    Skipping: insufficient samples ({samples_per_cond.to_dict()})")
        continue

    # Aggregate counts by sample
    if count_layer and count_layer in adata_ct.layers:
        X = adata_ct.layers[count_layer]
    else:
        X = adata_ct.X

    if issparse(X):
        X = X.toarray()

    # Check if data is counts (integers)
    if not np.allclose(X, np.round(X)):
        # May be log-normalized, try to reverse-transform
        X = np.round(np.expm1(X) * (10000 / np.expm1(X).sum(axis=1, keepdims=True)) *
                      adata_ct.obs['total_counts'].values[:, None] / 10000)
        X = np.clip(X, 0, None)

    gene_names = adata_ct.var_names.tolist()

    # Pseudobulk: sum by sample_id
    sample_ids = adata_ct.obs['sample_id'].values
    unique_samples = sorted(set(sample_ids))

    pb_counts = np.zeros((len(unique_samples), len(gene_names)))
    pb_meta = []
    for i, sid in enumerate(unique_samples):
        mask_s = sample_ids == sid
        pb_counts[i] = X[mask_s].sum(axis=0)
        cond = adata_ct.obs.loc[mask_s, 'condition'].iloc[0]
        n_cells = mask_s.sum()
        pb_meta.append({'sample_id': sid, 'condition': cond, 'n_cells': n_cells})

    pb_counts_df = pd.DataFrame(pb_counts.astype(int),
                                index=[m['sample_id'] for m in pb_meta],
                                columns=gene_names)
    pb_meta_df = pd.DataFrame(pb_meta).set_index('sample_id')

    # Filter lowly expressed genes (at least >= 10 counts in at least 3 samples)
    gene_mask = (pb_counts_df >= 10).sum(axis=0) >= 3
    pb_counts_df = pb_counts_df.loc[:, gene_mask]

    if pb_counts_df.shape[1] < 100:
        print(f"    Skipping: only {pb_counts_df.shape[1]} genes after filtering")
        continue

    # Filter very small samples
    sample_mask = pb_meta_df['n_cells'] >= 20
    pb_counts_df = pb_counts_df.loc[sample_mask]
    pb_meta_df = pb_meta_df.loc[sample_mask]

    # Confirm both groups have samples
    if pb_meta_df['condition'].nunique() < 2:
        print(f"    Skipping: only one condition after filtering")
        continue

    n_woi_ct = (pb_meta_df['condition'] == 'WOI Control').sum()
    n_as_ct = (pb_meta_df['condition'] == 'AS').sum()
    print(f"    Pseudobulk: {pb_counts_df.shape[0]} samples ({n_woi_ct} WOI, {n_as_ct} AS), {pb_counts_df.shape[1]} genes")

    # PyDESeq2
    try:
        dds = DeseqDataSet(
            counts=pb_counts_df,
            metadata=pb_meta_df,
            design="~ condition",
        )
        dds.deseq2()

        stat_res = DeseqStats(dds, contrast=['condition', 'AS', 'WOI Control'])
        stat_res.summary()

        results_df = stat_res.results_df.copy()
        results_df = results_df.dropna(subset=['padj'])
        results_df = results_df.sort_values('padj')

        # Statistics
        sig_up = ((results_df['padj'] < 0.05) & (results_df['log2FoldChange'] > 0.25)).sum()
        sig_down = ((results_df['padj'] < 0.05) & (results_df['log2FoldChange'] < -0.25)).sum()

        print(f"    DEG results: {sig_up} up, {sig_down} down (padj<0.05, |log2FC|>0.25)")

        all_deg_results[ct] = results_df
        deg_summary[ct] = {
            'n_samples_WOI': int(n_woi_ct),
            'n_samples_AS': int(n_as_ct),
            'n_genes_tested': int(len(results_df)),
            'n_sig_up': int(sig_up),
            'n_sig_down': int(sig_down),
            'top_up': results_df[results_df['log2FoldChange'] > 0].head(10).index.tolist(),
            'top_down': results_df[results_df['log2FoldChange'] < 0].head(10).index.tolist(),
        }

        # Save full DEG table
        results_df.to_csv(f'{OUTDIR}deg_{ct.replace("+", "plus").replace("/", "_")}_AS_vs_WOI.csv')

    except Exception as e:
        print(f"    PyDESeq2 error: {e}")
        continue

# ============================================================
# 5. DEG visualization
# ============================================================
print("\n--- DEG visualization ---")

# 5a. Volcano plots (per cell type)
n_ct = len(all_deg_results)
if n_ct > 0:
    n_cols = min(4, n_ct)
    n_rows = (n_ct + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
    if n_ct == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    for idx, (ct, results_df) in enumerate(all_deg_results.items()):
        ax = axes[idx]

        # Categorize
        sig_up = (results_df['padj'] < 0.05) & (results_df['log2FoldChange'] > 0.25)
        sig_down = (results_df['padj'] < 0.05) & (results_df['log2FoldChange'] < -0.25)
        ns = ~(sig_up | sig_down)

        ax.scatter(results_df.loc[ns, 'log2FoldChange'],
                  -np.log10(results_df.loc[ns, 'padj']),
                  c='grey', s=3, alpha=0.3, rasterized=True)
        ax.scatter(results_df.loc[sig_up, 'log2FoldChange'],
                  -np.log10(results_df.loc[sig_up, 'padj']),
                  c='#E64B35', s=5, alpha=0.5, rasterized=True)
        ax.scatter(results_df.loc[sig_down, 'log2FoldChange'],
                  -np.log10(results_df.loc[sig_down, 'padj']),
                  c='#4DBBD5', s=5, alpha=0.5, rasterized=True)

        # Annotate top genes
        top_genes = pd.concat([
            results_df[sig_up].head(5),
            results_df[sig_down].head(5)
        ])
        for gene in top_genes.index:
            row = results_df.loc[gene]
            ax.annotate(gene, (row['log2FoldChange'], -np.log10(row['padj'])),
                       fontsize=6, alpha=0.8,
                       arrowprops=dict(arrowstyle='-', color='grey', lw=0.5))

        ax.axhline(-np.log10(0.05), color='grey', linestyle='--', linewidth=0.5)
        ax.axvline(0.25, color='grey', linestyle='--', linewidth=0.5)
        ax.axvline(-0.25, color='grey', linestyle='--', linewidth=0.5)
        ax.set_xlabel('log2FC (AS vs WOI)')
        ax.set_ylabel('-log10(padj)')
        ax.set_title(f'{ct}\n{chr(8593)}{sig_up.sum()} {chr(8595)}{sig_down.sum()}')

    for idx in range(len(all_deg_results), len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle('Differential Gene Expression: AS vs WOI Control (Pseudobulk)', fontsize=14)
    plt.tight_layout()
    plt.savefig(FIGDIR + 'volcano_plots_AS_vs_WOI.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("  volcano_plots_AS_vs_WOI.png")

# 5b. DEG count summary bar chart
if deg_summary:
    fig, ax = plt.subplots(figsize=(12, 5))
    ct_names = list(deg_summary.keys())
    up_counts = [deg_summary[ct]['n_sig_up'] for ct in ct_names]
    down_counts = [-deg_summary[ct]['n_sig_down'] for ct in ct_names]

    x = np.arange(len(ct_names))
    ax.bar(x, up_counts, color='#E64B35', alpha=0.8, label='Up in AS')
    ax.bar(x, down_counts, color='#4DBBD5', alpha=0.8, label='Down in AS')
    ax.set_xticks(x)
    ax.set_xticklabels(ct_names, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Number of DEGs')
    ax.set_title('DEG Counts per Cell Type (AS vs WOI Control, padj<0.05, |log2FC|>0.25)')
    ax.axhline(0, color='black', linewidth=0.5)
    ax.legend()
    plt.tight_layout()
    plt.savefig(FIGDIR + 'deg_counts_barplot.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("  deg_counts_barplot.png")

# ============================================================
# 6. Stem cell-specific analysis
# ============================================================
print("\n--- Stem cell-specific analysis ---")

# Changes in stem cell markers in AS vs WOI
stem_markers = {
    'Stemness maintenance': ['SOX9', 'LGR5', 'PROM1', 'AXIN2', 'ALDH1A1', 'TERT', 'BMI1', 'NANOG'],
    'Wnt pathway': ['WNT2', 'WNT3A', 'WNT4', 'WNT5A', 'WNT7A', 'RSPO1', 'RSPO3', 'CTNNB1',
                'LEF1', 'TCF7L2', 'SFRP1', 'SFRP4', 'DKK1', 'WIF1'],
    'Notch pathway': ['NOTCH1', 'NOTCH2', 'NOTCH3', 'JAG1', 'JAG2', 'DLL1', 'DLL4',
                  'HES1', 'HEY1', 'RBPJ'],
    'TGF-beta pathway': ['TGFB1', 'TGFB2', 'TGFB3', 'TGFBR1', 'TGFBR2', 'SMAD2', 'SMAD3',
                   'SMAD4', 'SMAD7', 'ACVR1', 'BMP2', 'BMP4', 'BMP7'],
    'Fibrosis': ['COL1A1', 'COL1A2', 'COL3A1', 'FN1', 'ACTA2', 'FAP', 'POSTN',
              'CTGF', 'SERPINE1', 'MMP2', 'MMP9', 'TIMP1', 'TIMP2'],
}

# Check which genes are present in the data
available_genes = set(adata_disease.var_names)
stem_marker_results = {}

for category, genes in stem_markers.items():
    present_genes = [g for g in genes if g in available_genes]
    stem_marker_results[category] = {
        'total': len(genes),
        'present': len(present_genes),
        'genes': present_genes
    }
    print(f"  {category}: {len(present_genes)}/{len(genes)} genes available")

# Annotate pathway-related genes in stem cell DEGs
for ct in ['SOX9+LGR5+_stem', 'CD133+_progenitor', 'Basalis_fibroblast']:
    if ct in all_deg_results:
        deg = all_deg_results[ct]
        sig_deg = deg[(deg['padj'] < 0.05) & (abs(deg['log2FoldChange']) > 0.25)]

        print(f"\n  Pathway-related DEGs in {ct}:")
        for category, info in stem_marker_results.items():
            pathway_degs = [g for g in info['genes'] if g in sig_deg.index]
            if pathway_degs:
                for g in pathway_degs:
                    row = sig_deg.loc[g]
                    direction = 'up' if row['log2FoldChange'] > 0 else 'down'
                    print(f"    {category}: {g} {direction} (log2FC={row['log2FoldChange']:.2f}, padj={row['padj']:.2e})")

# ============================================================
# 7. RIF comparison (from Phase 4 Visium deconvolution)
# ============================================================
print("\n--- RIF comparison (Visium deconvolution) ---")

deconv_df = pd.read_csv('Phase_output/phase4_spatial/deconvolution_proportions.csv', index_col=0)

# Read Visium metadata to get condition and sample (from deconvolution h5ad, index is consistent)
adata_deconv = ad.read_h5ad('Phase_output/phase4_spatial/spatial_deconvolution.h5ad', backed='r')
vis_meta = adata_deconv.obs[['condition', 'sample']].copy()

# Match condition and sample to deconv results (via shared index)
common_idx = deconv_df.index.intersection(vis_meta.index)
deconv_df = deconv_df.loc[common_idx]
deconv_df['condition'] = vis_meta.loc[common_idx, 'condition'].values
deconv_df['sample'] = vis_meta.loc[common_idx, 'sample'].values
print(f"  Matched spots: {len(common_idx)}, CTR={(deconv_df['condition']=='CTR').sum()}, RIF={(deconv_df['condition']=='RIF').sum()}")
del adata_deconv

# Clean column names (q05cell_abundance_w_sf_XXX -> XXX)
rename_cols = {}
for col in deconv_df.columns:
    if col.startswith('q05cell_abundance_w_sf_'):
        rename_cols[col] = col.replace('q05cell_abundance_w_sf_', '')
deconv_df = deconv_df.rename(columns=rename_cols)

# Get cell type columns
ct_cols = [c for c in deconv_df.columns if c not in ['condition', 'sample']]

# Aggregate by sample (mean deconvolution proportion per sample)
rif_prop = []
for sample_id in deconv_df['sample'].unique():
    mask_s = deconv_df['sample'] == sample_id
    cond = deconv_df.loc[mask_s, 'condition'].iloc[0]
    for ct in ct_cols:
        mean_val = deconv_df.loc[mask_s, ct].mean()
        rif_prop.append({
            'sample': sample_id,
            'condition': cond,
            'celltype': ct,
            'abundance': mean_val
        })
rif_prop_df = pd.DataFrame(rif_prop)

# Statistical testing
rif_results = []
for ct in ct_cols:
    ct_df = rif_prop_df[rif_prop_df['celltype'] == ct]
    ctr_vals = ct_df[ct_df['condition'] == 'CTR']['abundance'].values
    rif_vals = ct_df[ct_df['condition'] == 'RIF']['abundance'].values

    if len(ctr_vals) >= 3 and len(rif_vals) >= 3:
        stat, pval = stats.mannwhitneyu(ctr_vals, rif_vals, alternative='two-sided')
    else:
        pval = np.nan

    rif_results.append({
        'celltype': ct,
        'mean_CTR': np.mean(ctr_vals),
        'mean_RIF': np.mean(rif_vals),
        'log2FC_RIF_vs_CTR': np.log2((np.mean(rif_vals) + 1e-6) / (np.mean(ctr_vals) + 1e-6)),
        'pvalue': pval,
    })

rif_results_df = pd.DataFrame(rif_results)
valid_rif = ~rif_results_df['pvalue'].isna()
if valid_rif.sum() > 0:
    _, padj, _, _ = multipletests(rif_results_df.loc[valid_rif, 'pvalue'], method='fdr_bh')
    rif_results_df.loc[valid_rif, 'padj'] = padj

rif_results_df = rif_results_df.sort_values('pvalue')
print("\nRIF vs CTR deconvolution proportion changes:")
print(rif_results_df.to_string(index=False, float_format='%.4f'))

# ============================================================
# 8. Integrated disease comparison plot
# ============================================================
print("\n--- Integrated comparison plot ---")

# Combine AS and RIF proportion changes
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# AS proportion changes
ax = axes[0]
df_as = prop_results_df.sort_values('log2FC_AS_vs_WOI')
colors_as = ['#E64B35' if p < 0.05 else '#999999' for p in df_as['padj']]
ax.barh(range(len(df_as)), df_as['log2FC_AS_vs_WOI'], color=colors_as, height=0.7)
ax.set_yticks(range(len(df_as)))
ax.set_yticklabels(df_as['celltype'], fontsize=9)
ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
ax.set_xlabel('log2FC (AS vs WOI Control)')
ax.set_title('Asherman Syndrome')

# RIF proportion changes
ax = axes[1]
# Align cell type order
ct_order = df_as['celltype'].tolist()
rif_ordered = []
for ct in ct_order:
    match = rif_results_df[rif_results_df['celltype'] == ct]
    if len(match) > 0:
        rif_ordered.append(match.iloc[0])
    else:
        rif_ordered.append({'celltype': ct, 'log2FC_RIF_vs_CTR': 0, 'padj': 1.0})
rif_ordered_df = pd.DataFrame(rif_ordered)

colors_rif = ['#3C5488' if p < 0.05 else '#999999' for p in rif_ordered_df['padj']]
ax.barh(range(len(rif_ordered_df)), rif_ordered_df['log2FC_RIF_vs_CTR'], color=colors_rif, height=0.7)
ax.set_yticks(range(len(rif_ordered_df)))
ax.set_yticklabels(rif_ordered_df['celltype'], fontsize=9)
ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
ax.set_xlabel('log2FC (RIF vs CTR)')
ax.set_title('Recurrent Implantation Failure')

fig.suptitle('Cell Type Proportion Changes in Endometrial Diseases', fontsize=14)
plt.tight_layout()
plt.savefig(FIGDIR + 'disease_comparison_proportions.png', dpi=200, bbox_inches='tight')
plt.close()
print("  disease_comparison_proportions.png")

# ============================================================
# 9. Save results
# ============================================================
print("\n--- Saving results ---")

# Save proportion analysis
prop_results_df.to_csv(OUTDIR + 'proportion_changes_AS_vs_WOI.csv', index=False)
rif_results_df.to_csv(OUTDIR + 'proportion_changes_RIF_vs_CTR.csv', index=False)
print("  proportion_changes_AS_vs_WOI.csv")
print("  proportion_changes_RIF_vs_CTR.csv")

# Save DEG summary
phase5_report = {
    'phase': 5,
    'comparison_primary': 'WOI Control vs AS (GSE215968)',
    'comparison_secondary': 'CTR vs RIF (GSE287278 Visium deconvolution)',
    'n_cells_primary': int(adata_disease.shape[0]),
    'proportion_analysis': {
        ct: {
            'mean_WOI': float(row['mean_WOI']),
            'mean_AS': float(row['mean_AS']),
            'log2FC': float(row['log2FC_AS_vs_WOI']),
            'padj': float(row['padj']) if not np.isnan(row['padj']) else None,
        }
        for _, row in prop_results_df.iterrows()
        for ct in [row['celltype']]
    },
    'deg_summary': deg_summary,
    'rif_comparison': {
        row['celltype']: {
            'mean_CTR': float(row['mean_CTR']),
            'mean_RIF': float(row['mean_RIF']),
            'log2FC': float(row['log2FC_RIF_vs_CTR']),
            'padj': float(row['padj']) if not pd.isna(row['padj']) else None,
        }
        for _, row in rif_results_df.iterrows()
    }
}

with open(OUTDIR + 'phase5_step1_report.json', 'w') as f:
    json.dump(phase5_report, f, indent=2, default=str)
print("  phase5_step1_report.json")

print("\n" + "=" * 60)
print("Phase 5 Step 1 complete!")
print(f"  AS proportion changes: {(prop_results_df['padj'] < 0.05).sum()}/{len(prop_results_df)} significant")
print(f"  DEG analysis: {len(all_deg_results)} cell types")
total_up = sum(d.get('n_sig_up', 0) for d in deg_summary.values())
total_down = sum(d.get('n_sig_down', 0) for d in deg_summary.values())
print(f"  Total DEGs: {total_up} up + {total_down} down")
print("=" * 60)
