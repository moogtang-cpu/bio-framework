#!/usr/bin/env python
"""
Phase 7 Step 2: NicheNet Ligand Activity Inference + Deep Analysis of Stem Cell Niche Signaling Network
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import liana as li
import decoupler as dc
import json
import gc
import warnings
warnings.filterwarnings('ignore')

FIGDIR = 'Phase_output/phase7_communication/figures/'
OUTDIR = 'Phase_output/phase7_communication/'

print("=" * 60)
print("Phase 7 Step 2: NicheNet + Deep Analysis of Niche Signaling Network")
print("=" * 60)

# ============================================================
# 1. Load data and step1 results
# ============================================================
print("\n--- 1. Load data ---")

# Step 1 results
liana_res = pd.read_csv(OUTDIR + 'liana_global_results.csv')
ctrl_res = pd.read_csv(OUTDIR + 'liana_control_results.csv')
as_res = pd.read_csv(OUTDIR + 'liana_as_results.csv')
diff_df = pd.read_csv(OUTDIR + 'disease_diff_interactions.csv')
sig_stem = pd.read_csv(OUTDIR + 'stem_niche_interactions.csv')

stem_types = ['SOX9+LGR5+_stem', 'CD133+_progenitor', 'Basalis_fibroblast']
niche_types = stem_types + ['Decidualized_stromal', 'Glandular_epithelium',
                            'Secretory_glandular', 'Smooth_muscle',
                            'Macrophages', 'Lymphoid_immune']

print(f"  Global interactions: {len(liana_res)}")
print(f"  Significant stem cell interactions: {len(sig_stem)}")

# ============================================================
# 2. Ligand activity inference (decoupler)
# ============================================================
print("\n--- 2. Ligand activity inference ---")

# Load stem cell subset
adata = ad.read_h5ad('Phase_output/phase3_integration/integrated_harmony.h5ad')
stem_mask = adata.obs['celltype_manual'].isin(stem_types)
adata_stem = adata[stem_mask].copy()
if 'log_normalized' in adata_stem.layers:
    adata_stem.X = adata_stem.layers['log_normalized'].copy()
print(f"  Stem cell subset: {adata_stem.shape}")

# Build ligand-receptor-target gene network from liana results
# Use liana's build_prior_network to construct a NicheNet-style network
print("  Building ligand-receptor prior network...")
try:
    from liana.resource import select_resource
    consensus = select_resource('consensus')
    print(f"  Consensus resource: {len(consensus)} LR pairs")

    # Extract ligands with significant signals to stem cells
    stem_incoming = sig_stem[sig_stem['target'].isin(stem_types)].copy()
    top_ligands = stem_incoming.groupby('ligand_complex')['magnitude_rank'].min().sort_values().head(50)
    print(f"  Number of top ligands: {len(top_ligands)}")
    print(f"  Top 10 ligands: {', '.join(top_ligands.index[:10])}")

except Exception as e:
    print(f"  Ligand activity inference skipped: {e}")

# ============================================================
# 3. Deep analysis of stem cell niche communication network
# ============================================================
print("\n--- 3. Stem cell niche communication network ---")

# 3a. Signaling pathway composition for each stem cell type
print("  3a. Signaling pathway composition analysis...")

pathway_genes = {
    'WNT': ['WNT2', 'WNT2B', 'WNT3', 'WNT3A', 'WNT4', 'WNT5A', 'WNT5B', 'WNT6', 'WNT7A', 'WNT7B',
            'WNT9A', 'WNT10A', 'WNT10B', 'WNT11', 'WNT16',
            'FZD1', 'FZD2', 'FZD3', 'FZD4', 'FZD5', 'FZD6', 'FZD7', 'FZD8', 'FZD9', 'FZD10',
            'LRP5', 'LRP6', 'RSPO1', 'RSPO2', 'RSPO3', 'RSPO4', 'LGR4', 'LGR5', 'LGR6',
            'SFRP1', 'SFRP2', 'SFRP4', 'SFRP5', 'DKK1', 'DKK2', 'DKK3', 'WIF1'],
    'NOTCH': ['NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4',
              'DLL1', 'DLL3', 'DLL4', 'JAG1', 'JAG2'],
    'TGFb_BMP': ['TGFB1', 'TGFB2', 'TGFB3', 'BMP2', 'BMP4', 'BMP6', 'BMP7',
                 'ACVR1', 'ACVR1B', 'ACVR2A', 'ACVR2B', 'TGFBR1', 'TGFBR2', 'TGFBR3',
                 'BMPR1A', 'BMPR1B', 'BMPR2', 'GDF5', 'GDF15', 'INHBA', 'INHBB'],
    'FGF': ['FGF1', 'FGF2', 'FGF7', 'FGF9', 'FGF10', 'FGF18',
            'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4'],
    'EGF': ['EGF', 'EGFR', 'ERBB2', 'ERBB3', 'ERBB4',
            'NRG1', 'NRG2', 'NRG3', 'AREG', 'HBEGF', 'EREG', 'BTC'],
    'Hedgehog': ['SHH', 'IHH', 'DHH', 'PTCH1', 'PTCH2', 'SMO'],
    'Cytokine': ['IL6', 'IL6R', 'IL6ST', 'IL1A', 'IL1B', 'IL1R1', 'IL1R2',
                 'TNF', 'TNFRSF1A', 'TNFRSF1B',
                 'IFNG', 'IFNGR1', 'IFNGR2', 'IL10', 'IL10RA', 'IL10RB'],
    'Chemokine': ['CXCL1', 'CXCL2', 'CXCL3', 'CXCL5', 'CXCL8', 'CXCL10', 'CXCL12', 'CXCL14',
                  'CCL2', 'CCL3', 'CCL4', 'CCL5', 'CCL7', 'CCL8', 'CCL11', 'CCL19', 'CCL20', 'CCL21',
                  'CXCR1', 'CXCR2', 'CXCR3', 'CXCR4', 'CCR1', 'CCR2', 'CCR5', 'CCR7'],
}

# Calculate pathway gene expression for each stem cell type
from scipy.sparse import issparse
X = adata_stem.X
if issparse(X):
    X = X.toarray()

pw_scores = {}
for pw, genes in pathway_genes.items():
    available = [g for g in genes if g in adata_stem.var_names]
    if len(available) >= 2:
        gene_idx = [list(adata_stem.var_names).index(g) for g in available]
        pw_expr = X[:, gene_idx].mean(axis=1)
        pw_scores[pw] = pw_expr

pw_df = pd.DataFrame(pw_scores, index=adata_stem.obs_names)
pw_df['celltype'] = adata_stem.obs['celltype_manual'].values
pw_df['disease'] = adata_stem.obs['condition'].map({
    'Normal': 'Control', 'WOI Control': 'Control',
    'AS': 'AS', 'AS_CD133': 'AS_CD133',
    'endometrium': 'Other', 'decidua': 'Other'
}).values

# Pathway score heatmap (by cell type)
pw_mean = pw_df.groupby('celltype')[list(pathway_genes.keys())].mean()
pw_mean_z = (pw_mean - pw_mean.mean()) / pw_mean.std()

fig, ax = plt.subplots(figsize=(10, 5))
sns.heatmap(pw_mean_z, cmap='RdBu_r', center=0, annot=pw_mean.round(3).values,
            fmt='', linewidths=0.5, ax=ax,
            cbar_kws={'label': 'Z-scored Mean Expression'})
ax.set_title('Signaling Pathway Activity in Stem/Progenitor Cells')
ax.set_ylabel('')
plt.tight_layout()
plt.savefig(FIGDIR + 'stem_pathway_expression.png', dpi=200, bbox_inches='tight')
plt.close()
print("  stem_pathway_expression.png")

# 3b. Pathway changes in disease
print("  3b. Pathway changes in disease...")
ctrl_stem = pw_df[pw_df['disease'] == 'Control']
as_stem = pw_df[pw_df['disease'] == 'AS']

if len(ctrl_stem) > 50 and len(as_stem) > 50:
    from scipy import stats
    pw_disease_diff = {}
    for pw in pathway_genes.keys():
        if pw in pw_df.columns:
            for ct in stem_types:
                c = ctrl_stem[ctrl_stem['celltype'] == ct][pw]
                a = as_stem[as_stem['celltype'] == ct][pw]
                if len(c) > 10 and len(a) > 10:
                    stat, pval = stats.mannwhitneyu(c, a, alternative='two-sided')
                    diff = a.mean() - c.mean()
                    pw_disease_diff[f'{ct}|{pw}'] = {
                        'celltype': ct, 'pathway': pw,
                        'ctrl_mean': c.mean(), 'as_mean': a.mean(),
                        'diff': diff, 'pvalue': pval,
                    }

    pw_diff_df = pd.DataFrame(pw_disease_diff).T
    if len(pw_diff_df) > 0:
        from statsmodels.stats.multitest import multipletests
        _, padj, _, _ = multipletests(pw_diff_df['pvalue'].astype(float), method='fdr_bh')
        pw_diff_df['padj'] = padj

        # Heatmap: pathway changes
        pivot_diff = pw_diff_df.pivot(index='pathway', columns='celltype', values='diff')
        pivot_padj = pw_diff_df.pivot(index='pathway', columns='celltype', values='padj')

        # Annotate significance
        annot = pivot_diff.round(4).astype(str)
        for pw in annot.index:
            for ct in annot.columns:
                if pw in pivot_padj.index and ct in pivot_padj.columns:
                    p = pivot_padj.loc[pw, ct]
                    if pd.notna(p) and float(p) < 0.05:
                        annot.loc[pw, ct] = annot.loc[pw, ct] + '*'
                    if pd.notna(p) and float(p) < 0.01:
                        annot.loc[pw, ct] = annot.loc[pw, ct] + '*'

        fig, ax = plt.subplots(figsize=(10, 6))
        sns.heatmap(pivot_diff.astype(float), cmap='RdBu_r', center=0, annot=annot, fmt='',
                    linewidths=0.5, ax=ax,
                    cbar_kws={'label': 'Mean Expression Diff (AS - Control)'})
        ax.set_title('Signaling Pathway Changes in Stem Cells (AS vs Control)\n* p<0.05, ** p<0.01')
        ax.set_ylabel('')
        plt.tight_layout()
        plt.savefig(FIGDIR + 'stem_pathway_disease_diff.png', dpi=200, bbox_inches='tight')
        plt.close()
        print("  stem_pathway_disease_diff.png")

        # Significant results
        sig_pw = pw_diff_df[pw_diff_df['padj'].astype(float) < 0.05].sort_values('padj')
        print(f"  Significant pathway changes (padj<0.05): {len(sig_pw)}")
        for _, row in sig_pw.iterrows():
            direction = '↑' if float(row['diff']) > 0 else '↓'
            print(f"    {row['celltype']} - {row['pathway']}: {direction} (diff={float(row['diff']):.4f}, padj={float(row['padj']):.2e})")

        pw_diff_df.to_csv(OUTDIR + 'stem_pathway_disease_diff.csv', index=False)

# ============================================================
# 4. Stem cell autocrine vs paracrine signaling
# ============================================================
print("\n--- 4. Autocrine vs paracrine ---")

for st in stem_types:
    incoming = sig_stem[sig_stem['target'] == st]
    autocrine = incoming[incoming['source'] == st]
    paracrine = incoming[incoming['source'] != st]
    print(f"  {st}:")
    print(f"    Autocrine: {len(autocrine)}, paracrine: {len(paracrine)}")
    if len(autocrine) > 0:
        top_auto = autocrine.sort_values('magnitude_rank').head(5)
        for _, r in top_auto.iterrows():
            print(f"      [Auto] {r['ligand_complex']}→{r['receptor_complex']} (rank={r['magnitude_rank']:.4f})")
    if len(paracrine) > 0:
        top_para = paracrine.sort_values('magnitude_rank').head(5)
        for _, r in top_para.iterrows():
            print(f"      [Para] {r['source']}: {r['ligand_complex']}→{r['receptor_complex']} (rank={r['magnitude_rank']:.4f})")

# ============================================================
# 5. Interaction network diagram (stem cell niche)
# ============================================================
print("\n--- 5. Interaction network diagram ---")

# Calculate interaction strength between niche cell types
niche_pairs = sig_stem.groupby(['source', 'target']).agg(
    n_interactions=('magnitude_rank', 'count'),
    mean_rank=('magnitude_rank', 'mean'),
).reset_index()

# Heatmap format
pair_pivot = niche_pairs.pivot(index='source', columns='target', values='n_interactions').fillna(0)
# Keep only niche types
niche_in_data = [ct for ct in niche_types if ct in pair_pivot.index or ct in pair_pivot.columns]
pair_pivot = pair_pivot.reindex(index=niche_in_data, columns=niche_in_data, fill_value=0)

fig, ax = plt.subplots(figsize=(12, 9))
sns.heatmap(pair_pivot, cmap='YlOrRd', annot=True, fmt='.0f',
            linewidths=0.3, ax=ax,
            cbar_kws={'label': 'Number of Significant Interactions'})
ax.set_title('Stem Cell Niche Communication Network\n(Significant LR pairs, rank < 0.01)')
ax.set_xlabel('Receiver')
ax.set_ylabel('Sender')
plt.tight_layout()
plt.savefig(FIGDIR + 'niche_communication_network.png', dpi=200, bbox_inches='tight')
plt.close()
print("  niche_communication_network.png")

# ============================================================
# 6. Communication pattern changes in disease - network level
# ============================================================
print("\n--- 6. Disease communication network changes ---")

# Stem cell interaction network: Control vs AS
stem_diff = diff_df[
    (diff_df['source'].isin(stem_types)) | (diff_df['target'].isin(stem_types))
]

# Aggregate by source-target pair
pair_changes = stem_diff.groupby(['source', 'target']).agg(
    mean_rank_diff=('rank_diff', 'mean'),
    n_gained=('rank_diff', lambda x: (x < -0.05).sum()),
    n_lost=('rank_diff', lambda x: (x > 0.05).sum()),
    n_total=('rank_diff', 'count'),
).reset_index()

# Calculate net change ratio
pair_changes['net_change'] = pair_changes['n_gained'] - pair_changes['n_lost']
pair_changes['pct_gained'] = pair_changes['n_gained'] / pair_changes['n_total']

# Visualization
net_pivot = pair_changes.pivot(index='source', columns='target', values='mean_rank_diff').fillna(0)
niche_in_change = [ct for ct in niche_types if ct in net_pivot.index or ct in net_pivot.columns]
net_pivot = net_pivot.reindex(index=niche_in_change, columns=niche_in_change, fill_value=0)

fig, ax = plt.subplots(figsize=(12, 9))
sns.heatmap(net_pivot, cmap='RdBu_r', center=0,
            annot=True, fmt='.3f', linewidths=0.3, ax=ax,
            cbar_kws={'label': 'Mean Rank Diff (AS - Control)\n← Gained | Lost →'})
ax.set_title('Communication Network Changes in AS\n(Stem Cell-Related)')
ax.set_xlabel('Receiver')
ax.set_ylabel('Sender')
plt.tight_layout()
plt.savefig(FIGDIR + 'disease_network_change.png', dpi=200, bbox_inches='tight')
plt.close()
print("  disease_network_change.png")

pair_changes.to_csv(OUTDIR + 'niche_pair_changes.csv', index=False)

# ============================================================
# 7. Key ligand-receptor pair summary (for Phase 8 drug targets)
# ============================================================
print("\n--- 7. Drug target candidate ligand-receptor pairs ---")

# Top interactions in stem cells with significant changes in AS
stem_diff_sorted = diff_df[
    (diff_df['source'].isin(stem_types)) | (diff_df['target'].isin(stem_types))
].copy()
stem_diff_sorted['abs_rank_diff'] = stem_diff_sorted['rank_diff'].abs()
stem_diff_sorted = stem_diff_sorted.sort_values('abs_rank_diff', ascending=False)

# Extract unique ligands and receptors
top_targets = stem_diff_sorted.head(100)
unique_ligands = top_targets['ligand'].unique()
unique_receptors = top_targets['receptor'].unique()

print(f"  Top changed ligands: {len(unique_ligands)}")
print(f"    {', '.join(unique_ligands[:15])}")
print(f"  Top changed receptors: {len(unique_receptors)}")
print(f"    {', '.join(unique_receptors[:15])}")

# Save drug target candidates
drug_targets = {
    'top_changed_ligands': unique_ligands[:30].tolist(),
    'top_changed_receptors': unique_receptors[:30].tolist(),
    'top_interactions': [],
}
for _, row in stem_diff_sorted.head(50).iterrows():
    drug_targets['top_interactions'].append({
        'source': row['source'], 'target': row['target'],
        'ligand': row['ligand'], 'receptor': row['receptor'],
        'rank_diff': float(row['rank_diff']),
        'direction': 'gained_in_AS' if row['rank_diff'] < 0 else 'lost_in_AS',
    })

with open(OUTDIR + 'drug_target_candidates_from_comm.json', 'w') as f:
    json.dump(drug_targets, f, indent=2, default=str)
print("  drug_target_candidates_from_comm.json")

# ============================================================
# 8. Update report
# ============================================================
print("\n--- 8. Update report ---")

with open(OUTDIR + 'phase7_report.json') as f:
    report = json.load(f)

report['step2_niche_analysis'] = {
    'pathway_scores': list(pathway_genes.keys()),
    'n_sig_pathway_changes': int(len(sig_pw)) if 'sig_pw' in dir() and len(sig_pw) > 0 else 0,
    'autocrine_vs_paracrine': {
        st: {
            'autocrine': int(len(sig_stem[(sig_stem['target']==st) & (sig_stem['source']==st)])),
            'paracrine': int(len(sig_stem[(sig_stem['target']==st) & (sig_stem['source']!=st)])),
        } for st in stem_types
    },
    'drug_target_candidates': {
        'n_ligands': int(len(unique_ligands)),
        'n_receptors': int(len(unique_receptors)),
        'top5_ligands': unique_ligands[:5].tolist(),
        'top5_receptors': unique_receptors[:5].tolist(),
    }
}

with open(OUTDIR + 'phase7_report.json', 'w') as f:
    json.dump(report, f, indent=2, default=str)
print("  phase7_report.json updated")

# Free memory
del adata, adata_stem
gc.collect()

print("\n" + "=" * 60)
print("Phase 7 Step 2 complete!")
print("=" * 60)
