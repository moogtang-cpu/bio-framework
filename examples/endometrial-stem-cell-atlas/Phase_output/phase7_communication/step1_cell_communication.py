#!/usr/bin/env python
"""
Phase 7: Cell-Cell Communication and Microenvironment Analysis
Using liana (rank_aggregate multi-method consensus) + disease comparison + stem cell niche focus
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
from liana.resource import select_resource
from scipy import stats
from statsmodels.stats.multitest import multipletests
import json
import os
import gc
import warnings
warnings.filterwarnings('ignore')

FIGDIR = 'Phase_output/phase7_communication/figures/'
OUTDIR = 'Phase_output/phase7_communication/'

print("=" * 60)
print("Phase 7: Cell-Cell Communication and Microenvironment Analysis")
print("=" * 60)

# ============================================================
# 1. Load data
# ============================================================
print("\n--- 1. Load data ---")
adata = ad.read_h5ad('Phase_output/phase3_integration/integrated_harmony.h5ad')
print(f"  Total data: {adata.shape}")

# Use the log_normalized layer
if 'log_normalized' in adata.layers:
    adata.X = adata.layers['log_normalized'].copy()
    print("  Using log_normalized layer")

# Ensure var_names are unique
adata.var_names_make_unique()

# ============================================================
# 2. Global communication analysis (downsampled for efficiency)
# ============================================================
print("\n--- 2. Global communication analysis ---")

# Downsample: at most 5000 cells per cell type
np.random.seed(42)
keep_idx = []
for ct in adata.obs['celltype_manual'].unique():
    ct_idx = np.where(adata.obs['celltype_manual'] == ct)[0]
    if len(ct_idx) > 5000:
        ct_idx = np.random.choice(ct_idx, 5000, replace=False)
    keep_idx.extend(ct_idx)
keep_idx = sorted(keep_idx)
adata_sub = adata[keep_idx].copy()
print(f"  Downsampled: {adata_sub.shape[0]} cells")
for ct in sorted(adata_sub.obs['celltype_manual'].unique()):
    n = (adata_sub.obs['celltype_manual'] == ct).sum()
    print(f"    {ct}: {n}")

# liana rank_aggregate: multi-method consensus scoring
print("\n  Running liana rank_aggregate (consensus resource)...")
li.mt.rank_aggregate(
    adata_sub,
    groupby='celltype_manual',
    resource_name='consensus',
    expr_prop=0.1,       # expressed in at least 10% of cells
    verbose=True,
    use_raw=False,
    n_perms=1000,
)

# Extract results
liana_res = adata_sub.uns['liana_res'].copy()
print(f"  Total {len(liana_res)} ligand-receptor pairs")
print(f"  Significant interactions (magnitude_rank < 0.01): {(liana_res['magnitude_rank'] < 0.01).sum()}")

# Save complete results
liana_res.to_csv(OUTDIR + 'liana_global_results.csv', index=False)

# ============================================================
# 3. Stem cell niche communication focus
# ============================================================
print("\n--- 3. Stem cell niche communication ---")

stem_types = ['SOX9+LGR5+_stem', 'CD133+_progenitor', 'Basalis_fibroblast']
niche_types = stem_types + ['Decidualized_stromal', 'Glandular_epithelium',
                            'Secretory_glandular', 'Smooth_muscle',
                            'Macrophages', 'Lymphoid_immune']

# Filter: interactions involving stem cells (as source or target)
stem_lr = liana_res[
    (liana_res['source'].isin(stem_types)) |
    (liana_res['target'].isin(stem_types))
].copy()
print(f"  Interactions involving stem cells: {len(stem_lr)}")

# Significant stem cell interactions
sig_stem = stem_lr[stem_lr['magnitude_rank'] < 0.01].copy()
print(f"  Significant (rank<0.01): {len(sig_stem)}")

# Top 20 signals received by stem cells
stem_incoming = sig_stem[sig_stem['target'].isin(stem_types)].sort_values('magnitude_rank')
print(f"\n  Top signals received by stem cells (target=stem):")
for _, row in stem_incoming.head(15).iterrows():
    print(f"    {row['source']} → {row['target']}: {row['ligand_complex']}→{row['receptor_complex']} (rank={row['magnitude_rank']:.4f})")

# Top 20 signals sent by stem cells
stem_outgoing = sig_stem[sig_stem['source'].isin(stem_types)].sort_values('magnitude_rank')
print(f"\n  Top signals sent by stem cells (source=stem):")
for _, row in stem_outgoing.head(15).iterrows():
    print(f"    {row['source']} → {row['target']}: {row['ligand_complex']}→{row['receptor_complex']} (rank={row['magnitude_rank']:.4f})")

# Save stem cell niche interactions
sig_stem.to_csv(OUTDIR + 'stem_niche_interactions.csv', index=False)

# ============================================================
# 4. Key signaling pathway extraction
# ============================================================
print("\n--- 4. Key signaling pathways ---")

# Define stem cell-related pathway keywords
pathway_keywords = {
    'WNT': ['WNT', 'FZD', 'LRP', 'RSPO', 'LGR', 'DKK', 'SFRP', 'WIF'],
    'NOTCH': ['NOTCH', 'DLL', 'JAG', 'HES', 'HEY'],
    'BMP_TGFb': ['BMP', 'TGFB', 'ACVR', 'SMAD', 'GDF', 'INHB', 'NODAL'],
    'FGF': ['FGF', 'FGFR'],
    'EGF': ['EGF', 'EGFR', 'ERBB', 'NRG', 'AREG', 'HBEGF'],
    'Hedgehog': ['SHH', 'IHH', 'DHH', 'PTCH', 'SMO'],
    'Inflammation': ['TNF', 'IL1', 'IL6', 'CXCL', 'CCL', 'IFNG', 'IL10', 'NFKB'],
    'ECM_Adhesion': ['COL', 'FN1', 'LAMA', 'ITGA', 'ITGB', 'CD44'],
}

pathway_lr = {}
for pw, keywords in pathway_keywords.items():
    mask = sig_stem.apply(lambda row: any(
        kw in str(row['ligand_complex']).upper() or kw in str(row['receptor_complex']).upper()
        for kw in keywords
    ), axis=1)
    pathway_lr[pw] = sig_stem[mask]
    n = len(pathway_lr[pw])
    if n > 0:
        print(f"  {pw}: {n} significant interactions")

# ============================================================
# 5. Disease-differential communication (Control vs AS)
# ============================================================
print("\n--- 5. Disease-differential communication ---")

# Groups: WOI Control vs AS (from GSE215968)
adata_sub.obs['disease_group'] = adata_sub.obs['condition'].map({
    'Normal': 'Control', 'WOI Control': 'Control',
    'AS': 'AS', 'AS_CD133': 'AS_CD133',
    'endometrium': 'Other', 'decidua': 'Other'
})

# Run liana separately for Control and AS
ctrl_mask = adata_sub.obs['disease_group'] == 'Control'
as_mask = adata_sub.obs['disease_group'] == 'AS'

print(f"  Control: {ctrl_mask.sum()} cells")
print(f"  AS: {as_mask.sum()} cells")

# Check cell counts per cell type in each group
for group, mask in [('Control', ctrl_mask), ('AS', as_mask)]:
    print(f"\n  {group}:")
    for ct in sorted(adata_sub.obs['celltype_manual'].unique()):
        n = ((adata_sub.obs['celltype_manual'] == ct) & mask).sum()
        if n > 0:
            print(f"    {ct}: {n}")

# Run Control group
adata_ctrl = adata_sub[ctrl_mask].copy()
# Filter out cell types with too few cells
ct_counts_ctrl = adata_ctrl.obs['celltype_manual'].value_counts()
valid_cts_ctrl = ct_counts_ctrl[ct_counts_ctrl >= 30].index.tolist()
adata_ctrl = adata_ctrl[adata_ctrl.obs['celltype_manual'].isin(valid_cts_ctrl)].copy()
print(f"\n  Control valid types: {len(valid_cts_ctrl)}, cells: {adata_ctrl.shape[0]}")

print("  Running liana (Control)...")
li.mt.rank_aggregate(
    adata_ctrl,
    groupby='celltype_manual',
    resource_name='consensus',
    expr_prop=0.1,
    verbose=True,
    use_raw=False,
    n_perms=1000,
)
ctrl_res = adata_ctrl.uns['liana_res'].copy()
ctrl_res['group'] = 'Control'

# Run AS group
adata_as = adata_sub[as_mask].copy()
ct_counts_as = adata_as.obs['celltype_manual'].value_counts()
valid_cts_as = ct_counts_as[ct_counts_as >= 30].index.tolist()
adata_as = adata_as[adata_as.obs['celltype_manual'].isin(valid_cts_as)].copy()
print(f"\n  AS valid types: {len(valid_cts_as)}, cells: {adata_as.shape[0]}")

print("  Running liana (AS)...")
li.mt.rank_aggregate(
    adata_as,
    groupby='celltype_manual',
    resource_name='consensus',
    expr_prop=0.1,
    verbose=True,
    use_raw=False,
    n_perms=1000,
)
as_res = adata_as.uns['liana_res'].copy()
as_res['group'] = 'AS'

# Save
ctrl_res.to_csv(OUTDIR + 'liana_control_results.csv', index=False)
as_res.to_csv(OUTDIR + 'liana_as_results.csv', index=False)

# Merge and compare
print("\n  Disease differential analysis...")
# Create interaction IDs
for df in [ctrl_res, as_res]:
    df['interaction_id'] = (df['source'] + '|' + df['target'] + '|' +
                           df['ligand_complex'] + '|' + df['receptor_complex'])

# Find common interactions
common_ids = set(ctrl_res['interaction_id']) & set(as_res['interaction_id'])
print(f"  Common interactions: {len(common_ids)}")

# Compare magnitude_rank differences
diff_data = []
for iid in common_ids:
    c = ctrl_res[ctrl_res['interaction_id'] == iid].iloc[0]
    a = as_res[as_res['interaction_id'] == iid].iloc[0]
    diff_data.append({
        'source': c['source'], 'target': c['target'],
        'ligand': c['ligand_complex'], 'receptor': c['receptor_complex'],
        'ctrl_rank': c['magnitude_rank'], 'as_rank': a['magnitude_rank'],
        'rank_diff': a['magnitude_rank'] - c['magnitude_rank'],  # negative = stronger in AS
        'ctrl_specificity': c['specificity_rank'] if 'specificity_rank' in c else np.nan,
        'as_specificity': a['specificity_rank'] if 'specificity_rank' in a else np.nan,
    })

diff_df = pd.DataFrame(diff_data)
diff_df = diff_df.sort_values('rank_diff')

# Interactions gained in AS (rank_diff < 0 means rank is smaller/more significant in AS)
as_gained = diff_df[diff_df['rank_diff'] < -0.1].copy()
# Interactions lost in AS
as_lost = diff_df[diff_df['rank_diff'] > 0.1].copy()

print(f"  Gained in AS: {len(as_gained)}, lost: {len(as_lost)}")

# Focus on stem cell-related differences
stem_diff = diff_df[
    (diff_df['source'].isin(stem_types)) | (diff_df['target'].isin(stem_types))
].copy()
stem_gained = stem_diff[stem_diff['rank_diff'] < -0.05].sort_values('rank_diff')
stem_lost = stem_diff[stem_diff['rank_diff'] > 0.05].sort_values('rank_diff', ascending=False)

print(f"\n  Stem cell communication changes: gained {len(stem_gained)}, lost {len(stem_lost)}")

if len(stem_gained) > 0:
    print("  Top10 stem cell interactions gained in AS:")
    for _, row in stem_gained.head(10).iterrows():
        print(f"    {row['source']}→{row['target']}: {row['ligand']}→{row['receptor']} "
              f"(Ctrl:{row['ctrl_rank']:.3f} → AS:{row['as_rank']:.3f})")

if len(stem_lost) > 0:
    print("  Top10 stem cell interactions lost in AS:")
    for _, row in stem_lost.head(10).iterrows():
        print(f"    {row['source']}→{row['target']}: {row['ligand']}→{row['receptor']} "
              f"(Ctrl:{row['ctrl_rank']:.3f} → AS:{row['as_rank']:.3f})")

diff_df.to_csv(OUTDIR + 'disease_diff_interactions.csv', index=False)

# ============================================================
# 6. Visualization
# ============================================================
print("\n--- 6. Visualization ---")

# 6a. Global interaction count heatmap (source × target)
print("  6a. Interaction count heatmap...")
sig_global = liana_res[liana_res['magnitude_rank'] < 0.01]
interaction_counts = sig_global.groupby(['source', 'target']).size().reset_index(name='count')
count_pivot = interaction_counts.pivot(index='source', columns='target', values='count').fillna(0)

fig, ax = plt.subplots(figsize=(14, 10))
sns.heatmap(count_pivot, cmap='YlOrRd', annot=True, fmt='.0f',
            linewidths=0.3, ax=ax, cbar_kws={'label': 'Number of Significant Interactions'})
ax.set_title('Cell-Cell Communication Network (Significant Interactions, rank<0.01)')
ax.set_xlabel('Receiver')
ax.set_ylabel('Sender')
plt.tight_layout()
plt.savefig(FIGDIR + 'interaction_count_heatmap.png', dpi=200, bbox_inches='tight')
plt.close()
print("    interaction_count_heatmap.png")

# 6b. Stem cell niche top interactions dotplot
print("  6b. Stem cell niche dotplot...")
stem_top = sig_stem.sort_values('magnitude_rank').head(40).copy()
stem_top['lr_pair'] = stem_top['ligand_complex'] + ' → ' + stem_top['receptor_complex']
stem_top['cell_pair'] = stem_top['source'] + ' → ' + stem_top['target']
stem_top['-log10(rank)'] = -np.log10(stem_top['magnitude_rank'].clip(1e-10))

fig, ax = plt.subplots(figsize=(14, 10))
scatter = ax.scatter(
    stem_top['cell_pair'], stem_top['lr_pair'],
    s=stem_top['-log10(rank)'] * 30,
    c=stem_top['magnitude_rank'],
    cmap='viridis_r', alpha=0.8, edgecolors='k', linewidths=0.3,
)
plt.colorbar(scatter, ax=ax, label='Magnitude Rank')
ax.set_xlabel('Cell Pair (Sender → Receiver)')
ax.set_ylabel('Ligand → Receptor')
ax.set_title('Top Stem Cell Niche Interactions (rank_aggregate)')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(FIGDIR + 'stem_niche_dotplot.png', dpi=200, bbox_inches='tight')
plt.close()
print("    stem_niche_dotplot.png")

# 6c. Disease-differential communication - stem cell related
print("  6c. Disease-differential communication plot...")
stem_diff_plot = stem_diff.copy()
stem_diff_plot['interaction'] = (stem_diff_plot['source'].str[:10] + '→' +
                                 stem_diff_plot['target'].str[:10] + '\n' +
                                 stem_diff_plot['ligand'] + '→' + stem_diff_plot['receptor'])

# Take top 30 most significant changes
top_changes = pd.concat([
    stem_diff_plot.sort_values('rank_diff').head(15),
    stem_diff_plot.sort_values('rank_diff', ascending=False).head(15)
]).drop_duplicates()

if len(top_changes) > 0:
    fig, ax = plt.subplots(figsize=(12, max(6, len(top_changes) * 0.35)))
    colors = ['#d62728' if x < 0 else '#1f77b4' for x in top_changes['rank_diff']]
    y_pos = range(len(top_changes))

    ax.barh(y_pos, top_changes['rank_diff'].values, color=colors, alpha=0.8)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_changes['interaction'].values, fontsize=7)
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel('Rank Difference (AS - Control)\n← Gained in AS | Lost in AS →')
    ax.set_title('Stem Cell Communication Changes in AS')
    plt.tight_layout()
    plt.savefig(FIGDIR + 'disease_diff_stem_comm.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("    disease_diff_stem_comm.png")

# 6d. Pathway-level communication summary
print("  6d. Pathway-level communication summary...")
pathway_summary = []
for pw, pw_df in pathway_lr.items():
    if len(pw_df) == 0:
        continue
    # Count by source-target
    for direction in ['incoming', 'outgoing']:
        if direction == 'incoming':
            sub = pw_df[pw_df['target'].isin(stem_types)]
        else:
            sub = pw_df[pw_df['source'].isin(stem_types)]

        if len(sub) > 0:
            pathway_summary.append({
                'pathway': pw,
                'direction': direction,
                'n_interactions': len(sub),
                'mean_rank': sub['magnitude_rank'].mean(),
                'top_lr': sub.sort_values('magnitude_rank').iloc[0]['ligand_complex'] + '→' +
                          sub.sort_values('magnitude_rank').iloc[0]['receptor_complex'],
            })

pw_summary_df = pd.DataFrame(pathway_summary)
if len(pw_summary_df) > 0:
    fig, ax = plt.subplots(figsize=(10, 6))
    pw_pivot = pw_summary_df.pivot(index='pathway', columns='direction', values='n_interactions').fillna(0)
    pw_pivot.plot(kind='barh', ax=ax, color=['#ff7f0e', '#2ca02c'])
    ax.set_xlabel('Number of Significant Interactions')
    ax.set_title('Stem Cell Niche Pathway Communication')
    ax.legend(title='Direction', labels=['Incoming to Stem', 'Outgoing from Stem'])
    plt.tight_layout()
    plt.savefig(FIGDIR + 'pathway_communication_summary.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("    pathway_communication_summary.png")

# 6e. Communication network chord diagram (simplified as heatmap)
print("  6e. Fibrosis vs regeneration signaling...")
# Fibrosis-related interactions
fibrosis_kw = ['TGFB', 'COL', 'FN1', 'CTGF', 'PDGF', 'MMP', 'TIMP', 'LOX', 'SERPINE',
               'IL6', 'IL1B', 'TNF', 'CCL2', 'CXCL']
regen_kw = ['WNT', 'FZD', 'RSPO', 'LGR', 'NOTCH', 'DLL', 'JAG', 'EGF', 'FGF',
            'BMP', 'SHH', 'VEGF', 'HGF', 'IGF']

def classify_lr(row, keywords):
    lr = str(row['ligand_complex']).upper() + ' ' + str(row['receptor_complex']).upper()
    return any(kw in lr for kw in keywords)

# Classify within differential interactions
if len(stem_diff) > 0:
    stem_diff_c = stem_diff.copy()
    stem_diff_c['is_fibrosis'] = stem_diff_c.apply(lambda r: classify_lr(r, fibrosis_kw), axis=1)
    stem_diff_c['is_regen'] = stem_diff_c.apply(lambda r: classify_lr(r, regen_kw), axis=1)

    fib_data = stem_diff_c[stem_diff_c['is_fibrosis']]
    reg_data = stem_diff_c[stem_diff_c['is_regen']]

    print(f"\n  Fibrosis-related interactions: {len(fib_data)} (gained in AS: {(fib_data['rank_diff']<0).sum()}, lost: {(fib_data['rank_diff']>0).sum()})")
    print(f"  Regeneration-related interactions: {len(reg_data)} (gained in AS: {(reg_data['rank_diff']<0).sum()}, lost: {(reg_data['rank_diff']>0).sum()})")

    # Visualization: fibrosis vs regeneration rank_diff distribution
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for idx, (data, label, color) in enumerate([
        (fib_data, 'Fibrosis-related', '#d62728'),
        (reg_data, 'Regeneration-related', '#2ca02c')
    ]):
        if len(data) > 0:
            axes[idx].hist(data['rank_diff'], bins=20, color=color, alpha=0.7, edgecolor='black')
            axes[idx].axvline(0, color='black', linewidth=1, linestyle='--')
            axes[idx].set_xlabel('Rank Diff (AS-Control)')
            axes[idx].set_ylabel('Count')
            axes[idx].set_title(f'{label} Interactions\n(← Gained in AS | Lost in AS →)')
            med = data['rank_diff'].median()
            axes[idx].axvline(med, color='red', linewidth=1.5, linestyle='-',
                             label=f'Median={med:.3f}')
            axes[idx].legend()
        else:
            axes[idx].text(0.5, 0.5, 'No data', ha='center', va='center',
                          transform=axes[idx].transAxes)
            axes[idx].set_title(label)

    fig.suptitle('Fibrosis vs Regeneration Signaling Changes in Stem Cell Niche (AS)')
    plt.tight_layout()
    plt.savefig(FIGDIR + 'fibrosis_vs_regeneration.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("    fibrosis_vs_regeneration.png")

# ============================================================
# 7. Stem cell niche key findings summary
# ============================================================
print("\n--- 7. Key findings summary ---")

# Count interaction partners for each stem cell type
for st in stem_types:
    incoming = sig_stem[sig_stem['target'] == st]
    outgoing = sig_stem[sig_stem['source'] == st]
    print(f"\n  {st}:")
    print(f"    Incoming signals: {len(incoming)} (from {incoming['source'].nunique()} cell types)")
    print(f"    Outgoing signals: {len(outgoing)} (to {outgoing['target'].nunique()} cell types)")

    if len(incoming) > 0:
        top_senders = incoming.groupby('source').size().sort_values(ascending=False)
        print(f"    Main signal sources: {', '.join(f'{k}({v})' for k,v in top_senders.head(3).items())}")

# ============================================================
# 8. Save report
# ============================================================
print("\n--- 8. Save report ---")

# Communication statistics
report = {
    'phase': 7,
    'method': 'liana rank_aggregate (consensus resource, 4624 LR pairs)',
    'global_analysis': {
        'n_cells': int(adata_sub.shape[0]),
        'n_interactions_tested': int(len(liana_res)),
        'n_significant': int((liana_res['magnitude_rank'] < 0.01).sum()),
    },
    'stem_niche': {
        'n_stem_interactions': int(len(sig_stem)),
        'stem_incoming': int(len(sig_stem[sig_stem['target'].isin(stem_types)])),
        'stem_outgoing': int(len(sig_stem[sig_stem['source'].isin(stem_types)])),
    },
    'disease_comparison': {
        'n_common_interactions': int(len(common_ids)),
        'n_gained_in_AS': int(len(as_gained)),
        'n_lost_in_AS': int(len(as_lost)),
        'stem_gained': int(len(stem_gained)),
        'stem_lost': int(len(stem_lost)),
    },
    'pathway_communication': {
        pw: int(len(pw_df)) for pw, pw_df in pathway_lr.items() if len(pw_df) > 0
    },
}

with open(OUTDIR + 'phase7_report.json', 'w') as f:
    json.dump(report, f, indent=2, default=str)
print("  phase7_report.json")

# Free memory
del adata, adata_sub, adata_ctrl, adata_as
gc.collect()

print("\n" + "=" * 60)
print("Phase 7 Step 1 complete!")
print(f"  Global significant interactions: {report['global_analysis']['n_significant']}")
print(f"  Stem cell niche interactions: {report['stem_niche']['n_stem_interactions']}")
print(f"  Disease differences: gained {report['disease_comparison']['n_gained_in_AS']}, lost {report['disease_comparison']['n_lost_in_AS']}")
print("=" * 60)
