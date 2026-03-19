#!/usr/bin/env python
"""
Phase 6 Supplement: TF activity inference (decoupler ULM + CollecTRI from omnipath)
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import decoupler as dc
import omnipath as op
from scipy import stats
from scipy.sparse import issparse
from statsmodels.stats.multitest import multipletests
import json
import os
import warnings
warnings.filterwarnings('ignore')

FIGDIR = 'Phase_output/phase6_trajectory/figures/'
OUTDIR = 'Phase_output/phase6_trajectory/'

print("=" * 60)
print("Phase 6 Supplement: TF Activity Inference")
print("=" * 60)

# 1. Retrieve CollecTRI TF-target network
print("\n--- Retrieving CollecTRI Network ---")
collectri = op.interactions.CollecTRI.get(genesymbols=True)
print(f"  Raw: {len(collectri)} regulatory relationships")

# Use genesymbol columns
net = collectri[['source_genesymbol', 'target_genesymbol', 'is_stimulation', 'is_inhibition']].copy()
net = net.rename(columns={'source_genesymbol': 'source', 'target_genesymbol': 'target'})
# Filter out sources starting with COMPLEX
net = net[~net['source'].str.startswith('COMPLEX:')]
# Weight: stimulation=+1, inhibition=-1
net['weight'] = 1.0
net.loc[net['is_inhibition'] == True, 'weight'] = -1.0
net = net[['source', 'target', 'weight']].drop_duplicates(subset=['source', 'target'])
print(f"  Processed: {len(net)} entries, {net['source'].nunique()} TFs, {net['target'].nunique()} targets")

# 2. Load stem cell / epithelial subset
print("\n--- Loading Data ---")
cytotrace_types = ['SOX9+LGR5+_stem', 'CD133+_progenitor', 'Basalis_fibroblast',
                   'Glandular_epithelium', 'Secretory_glandular', 'Decidualized_stromal']

adata = ad.read_h5ad('Phase_output/phase3_integration/integrated_harmony.h5ad')
adata_tf = adata[adata.obs['celltype_manual'].isin(cytotrace_types)].copy()
print(f"  Subset: {adata_tf.shape}")
del adata

# Use log_normalized layer
if 'log_normalized' in adata_tf.layers:
    adata_tf.X = adata_tf.layers['log_normalized'].copy()

# 3. ULM TF activity inference
print("\n--- ULM TF Activity Inference ---")
# decoupler 2.x API: dc.mt.ulm(data, net, ...)
# net must have source/target/weight columns
dc.mt.ulm(
    data=adata_tf,
    net=net,
    verbose=True,
    raw=False,
)

tf_activities = adata_tf.obsm['score_ulm']
print(f"  Inferred activity for {tf_activities.shape[1]} TFs")

# 4. Top TFs per cell type
print("\n--- Cell Type-Specific TFs ---")
tf_top = {}
for ct in cytotrace_types:
    ct_mask = adata_tf.obs['celltype_manual'] == ct
    mean_act = tf_activities[ct_mask].mean(axis=0)
    top_tfs = mean_act.nlargest(10)
    tf_top[ct] = top_tfs.to_dict()
    print(f"  {ct} top5: {', '.join(top_tfs.index[:5].tolist())}")

# 5. TF activity heatmap
print("\n--- Visualization ---")

mean_tf_data = {}
for ct in cytotrace_types:
    ct_mask = adata_tf.obs['celltype_manual'] == ct
    mean_tf_data[ct] = tf_activities[ct_mask].mean(axis=0)
mean_tf = pd.DataFrame(mean_tf_data).T

# Top variable TFs
tf_var = mean_tf.var(axis=0).nlargest(30).index
plot_tf = mean_tf[tf_var]

fig, ax = plt.subplots(figsize=(14, 6))
sns.heatmap(plot_tf, cmap='RdBu_r', center=0, annot=True, fmt='.2f',
            linewidths=0.3, ax=ax, cbar_kws={'label': 'TF Activity (ULM)'})
ax.set_title('Top Variable TF Activities across Cell Types (CollecTRI)')
ax.set_ylabel('')
plt.tight_layout()
plt.savefig(FIGDIR + 'tf_activity_heatmap.png', dpi=200, bbox_inches='tight')
plt.close()
print("  tf_activity_heatmap.png")

# 6. Disease-differential TF activity
print("\n--- Disease-Differential TF Activity ---")
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
        print(f"  {ct}: insufficient cells, skipping")
        continue

    ctrl_mask = ct_data.obs['disease_group'] == 'Control'
    as_mask = ct_data.obs['disease_group'] == 'AS'

    ctrl_act = ct_data.obsm['score_ulm'][ctrl_mask.values]
    as_act = ct_data.obsm['score_ulm'][as_mask.values]

    if ctrl_act.shape[0] < 20 or as_act.shape[0] < 20:
        print(f"  {ct}: one group <20 cells, skipping")
        continue

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

    _, padj, _, _ = multipletests(tf_diff_df['pvalue'], method='fdr_bh')
    tf_diff_df['padj'] = padj
    tf_diff_df = tf_diff_df.sort_values('padj')

    sig_tfs = tf_diff_df[tf_diff_df['padj'] < 0.05]
    n_up = (sig_tfs['diff_AS_vs_Control'] > 0).sum()
    n_down = (sig_tfs['diff_AS_vs_Control'] < 0).sum()
    print(f"  {ct}: {n_up}↑ + {n_down}↓ TFs (padj<0.05)")
    if len(sig_tfs) > 0:
        top_up = sig_tfs[sig_tfs['diff_AS_vs_Control'] > 0].head(5)['TF'].tolist()
        top_down = sig_tfs[sig_tfs['diff_AS_vs_Control'] < 0].head(5)['TF'].tolist()
        if top_up:
            print(f"    Top up: {', '.join(top_up)}")
        if top_down:
            print(f"    Top down: {', '.join(top_down)}")

    disease_tf_results[ct] = tf_diff_df
    tf_diff_df.to_csv(f'{OUTDIR}tf_diff_{ct.replace("+","plus").replace("/","_")}_AS_vs_Control.csv', index=False)

# 7. Key TF comparison heatmap for stem cells
stem_key_tfs = set()
for ct in ['SOX9+LGR5+_stem', 'CD133+_progenitor', 'Basalis_fibroblast']:
    if ct in disease_tf_results:
        top10 = disease_tf_results[ct].head(15)['TF'].tolist()
        stem_key_tfs.update(top10)

if stem_key_tfs:
    stem_tfs = sorted(stem_key_tfs)[:25]
    plot_data = []
    for ct in ['SOX9+LGR5+_stem', 'CD133+_progenitor', 'Basalis_fibroblast',
               'Decidualized_stromal', 'Glandular_epithelium']:
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

        # Annotate significance
        sig_pivot = plot_df.pivot(index='TF', columns='celltype', values='padj')
        annot = pivot.round(2).astype(str)
        for tf in annot.index:
            for ct in annot.columns:
                if tf in sig_pivot.index and ct in sig_pivot.columns:
                    p = sig_pivot.loc[tf, ct]
                    if pd.notna(p) and p < 0.05:
                        annot.loc[tf, ct] = annot.loc[tf, ct] + '*'

        fig, ax = plt.subplots(figsize=(12, max(6, len(pivot) * 0.4)))
        sns.heatmap(pivot, cmap='RdBu_r', center=0, annot=annot, fmt='',
                    linewidths=0.3, ax=ax, cbar_kws={'label': 'TF Activity Diff (AS - Control)'})
        ax.set_title('Key TF Activity Changes in Stem/Progenitor Cells (AS vs Control)\n* = padj < 0.05')
        plt.tight_layout()
        plt.savefig(FIGDIR + 'tf_stem_disease_diff.png', dpi=200, bbox_inches='tight')
        plt.close()
        print("  tf_stem_disease_diff.png")

# 8. UMAP colored by TF activity (key TFs)
key_tfs_to_plot = ['SOX9', 'TP53', 'FOXO1', 'STAT3', 'NFKB1', 'HIF1A', 'MYC', 'JUN']
available_tfs = [tf for tf in key_tfs_to_plot if tf in tf_activities.columns]

if available_tfs:
    # UMAP required on adata_tf
    sc.pp.neighbors(adata_tf, use_rep='X_pca_harmony', n_neighbors=30, random_state=42)
    sc.tl.umap(adata_tf, random_state=42)

    n_tfs = len(available_tfs)
    n_cols = min(4, n_tfs)
    n_rows = (n_tfs + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    if n_tfs == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    for idx, tf in enumerate(available_tfs):
        adata_tf.obs[f'TF_{tf}'] = tf_activities[tf].values
        sc.pl.umap(adata_tf, color=f'TF_{tf}', ax=axes[idx], show=False,
                   title=f'{tf} Activity', cmap='RdBu_r', vcenter=0, frameon=False)

    for idx in range(len(available_tfs), len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle('Key TF Activities (ULM)', fontsize=14)
    plt.tight_layout()
    plt.savefig(FIGDIR + 'tf_umap_key_tfs.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("  tf_umap_key_tfs.png")

# 9. Update report
print("\n--- Updating Report ---")
tf_report = {
    'method': 'decoupler ULM + CollecTRI (omnipath)',
    'n_tfs': int(tf_activities.shape[1]),
    'n_regulations': int(len(net)),
    'celltype_top_tfs': {ct: list(tfs.keys())[:5] for ct, tfs in tf_top.items()},
    'disease_diff': {
        ct: {
            'n_sig_up': int((df[df['padj'] < 0.05]['diff_AS_vs_Control'] > 0).sum()),
            'n_sig_down': int((df[df['padj'] < 0.05]['diff_AS_vs_Control'] < 0).sum()),
            'top_up': df[df['diff_AS_vs_Control'] > 0].head(5)['TF'].tolist(),
            'top_down': df[df['diff_AS_vs_Control'] < 0].head(5)['TF'].tolist(),
        }
        for ct, df in disease_tf_results.items()
    }
}

# Read existing report and update
report_path = OUTDIR + 'phase6_report.json'
with open(report_path) as f:
    report = json.load(f)
report['tf_analysis'] = tf_report
with open(report_path, 'w') as f:
    json.dump(report, f, indent=2, default=str)
print("  phase6_report.json updated")

print("\n" + "=" * 60)
print("Phase 6 TF Activity Supplement Complete!")
print(f"  TFs inferred: {tf_activities.shape[1]}")
print(f"  Disease differential: {len(disease_tf_results)} cell types")
print("=" * 60)
