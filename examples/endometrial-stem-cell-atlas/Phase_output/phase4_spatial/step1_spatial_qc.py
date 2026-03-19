#!/usr/bin/env python
"""
Phase 4 Step 1: Visium spatial transcriptomics QC and preprocessing
- GSE287278 RIF Visium: 8 samples (CTR1-4, RIF1-4)
- E-MTAB-9260: spatial coordinates not yet ready, skipped for now
"""

import scanpy as sc
import squidpy as sq
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json
import os
import warnings
warnings.filterwarnings('ignore')

FIGDIR = 'Phase_output/phase4_spatial/figures/'
sc.settings.figdir = FIGDIR
os.makedirs(FIGDIR, exist_ok=True)

# ============================================================
# 1. Load all 8 Visium samples
# ============================================================
print("=" * 60)
print("Phase 4 Step 1: Spatial QC & Preprocessing")
print("=" * 60)

base = 'data/organized/GSE287278_RIF_visium/'
samples = ['CTR1', 'CTR2', 'CTR3', 'CTR4', 'RIF1', 'RIF2', 'RIF3', 'RIF4']
adatas = {}

for s in samples:
    a = ad.read_h5ad(os.path.join(base, s + '.h5ad'))
    adatas[s] = a
    print(f"  {s}: {a.shape[0]} spots, condition={a.obs['condition'].iloc[0]}")

print(f"\nTotal: {sum(a.shape[0] for a in adatas.values())} spots")

# ============================================================
# 2. QC metric calculation
# ============================================================
print("\n--- QC metric calculation ---")

qc_stats = {}
for s, adata in adatas.items():
    # Ensure X is raw counts
    adata.var_names_make_unique()

    # Calculate QC metrics
    # Mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    # Ribosomal genes
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    # Hemoglobin genes
    adata.var['hb'] = adata.var_names.str.match('^HB[^(P)]')

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=['mt', 'ribo', 'hb'],
        percent_top=None, log1p=False, inplace=True
    )

    stats = {
        'n_spots': int(adata.shape[0]),
        'n_genes': int(adata.shape[1]),
        'median_total_counts': float(np.median(adata.obs['total_counts'])),
        'median_n_genes': float(np.median(adata.obs['n_genes_by_counts'])),
        'mean_pct_mt': float(adata.obs['pct_counts_mt'].mean()),
        'mean_pct_ribo': float(adata.obs['pct_counts_ribo'].mean()),
        'mean_pct_hb': float(adata.obs['pct_counts_hb'].mean()),
        'condition': adata.obs['condition'].iloc[0],
    }
    qc_stats[s] = stats
    print(f"  {s}: spots={stats['n_spots']}, medUMI={stats['median_total_counts']:.0f}, "
          f"medGenes={stats['median_n_genes']:.0f}, mt%={stats['mean_pct_mt']:.1f}%, "
          f"ribo%={stats['mean_pct_ribo']:.1f}%, hb%={stats['mean_pct_hb']:.1f}%")

# ============================================================
# 3. QC filtering
# ============================================================
print("\n--- QC filtering ---")
# Visium QC thresholds (lenient, as spot density varies widely across tissue regions)
MIN_COUNTS = 500
MIN_GENES = 200
MAX_MT_PCT = 25  # Visium mt% is typically higher (tissue sections contain necrotic regions)

filtered_adatas = {}
for s, adata in adatas.items():
    n_before = adata.shape[0]

    # Filter low-quality spots
    mask = (
        (adata.obs['total_counts'] >= MIN_COUNTS) &
        (adata.obs['n_genes_by_counts'] >= MIN_GENES) &
        (adata.obs['pct_counts_mt'] <= MAX_MT_PCT)
    )
    adata_f = adata[mask].copy()

    # Filter lowly expressed genes (expressed in at least 5 spots)
    sc.pp.filter_genes(adata_f, min_cells=5)

    n_after = adata_f.shape[0]
    n_genes_after = adata_f.shape[1]
    pct_kept = n_after / n_before * 100

    print(f"  {s}: {n_before} -> {n_after} spots ({pct_kept:.1f}%), {n_genes_after} genes")
    qc_stats[s]['n_spots_after_qc'] = n_after
    qc_stats[s]['pct_spots_kept'] = round(pct_kept, 1)
    qc_stats[s]['n_genes_after_qc'] = n_genes_after

    filtered_adatas[s] = adata_f

# ============================================================
# 4. Normalization and feature selection
# ============================================================
print("\n--- Normalization and feature selection ---")

for s, adata in filtered_adatas.items():
    # Save raw counts
    adata.layers['counts'] = adata.X.copy()

    # Normalization (log1p 10K)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['log_normalized'] = adata.X.copy()

    print(f"  {s}: normalization complete")

# ============================================================
# 5. Merge all samples
# ============================================================
print("\n--- Merging samples ---")

adata_list = list(filtered_adatas.values())
adata_spatial = ad.concat(adata_list, join='inner')
print(f"After merge: {adata_spatial.shape[0]} spots x {adata_spatial.shape[1]} genes")

# HVG selection (using counts layer)
adata_spatial.X = adata_spatial.layers['counts'].copy()
sc.pp.highly_variable_genes(
    adata_spatial,
    flavor='seurat_v3',
    n_top_genes=3000,
    batch_key='sample'
)
n_hvg = adata_spatial.var['highly_variable'].sum()
print(f"HVG: {n_hvg} genes")

# Restore log_normalized to X
adata_spatial.X = adata_spatial.layers['log_normalized'].copy()

# PCA
sc.pp.scale(adata_spatial, max_value=10)
sc.tl.pca(adata_spatial, n_comps=30, svd_solver='randomized', random_state=42,
          use_highly_variable=True)
print(f"PCA: 30 components")

# Harmony batch correction (between samples)
import harmonypy as hm
X_pca = adata_spatial.obsm['X_pca']
meta = adata_spatial.obs[['sample']].copy()
ho = hm.run_harmony(X_pca, meta, 'sample', max_iter_harmony=20, random_state=42)
Z_corr = ho.Z_corr
if hasattr(Z_corr, 'numpy'):
    Z_corr = Z_corr.cpu().numpy()
if Z_corr.shape[0] != adata_spatial.shape[0]:
    Z_corr = Z_corr.T
adata_spatial.obsm['X_pca_harmony'] = Z_corr
print(f"Harmony batch correction complete")

# UMAP (using Harmony-corrected PCs)
sc.pp.neighbors(adata_spatial, use_rep='X_pca_harmony', n_neighbors=15, random_state=42)
sc.tl.umap(adata_spatial, random_state=42)
print(f"UMAP complete")

# Leiden clustering
for res in [0.3, 0.5, 0.8, 1.0]:
    sc.tl.leiden(adata_spatial, resolution=res, key_added=f'leiden_{res}', random_state=42)
    n_clust = adata_spatial.obs[f'leiden_{res}'].nunique()
    print(f"  Leiden res={res}: {n_clust} clusters")

# ============================================================
# 6. QC visualization
# ============================================================
print("\n--- Generating QC plots ---")

# 6a. Spatial QC plot for each sample
fig, axes = plt.subplots(2, 4, figsize=(24, 12))
for idx, (s, adata) in enumerate(filtered_adatas.items()):
    ax = axes[idx // 4, idx % 4]
    coords = adata.obsm['spatial']
    sc_plot = ax.scatter(coords[:, 0], coords[:, 1],
                        c=adata.obs['total_counts'], cmap='viridis',
                        s=8, alpha=0.8)
    ax.set_title(f'{s} ({adata.obs["condition"].iloc[0]})\n{adata.shape[0]} spots')
    ax.set_aspect('equal')
    ax.invert_yaxis()
    plt.colorbar(sc_plot, ax=ax, label='Total counts', shrink=0.6)
fig.suptitle('GSE287278 Visium - Total Counts per Spot', fontsize=16)
plt.tight_layout()
plt.savefig(FIGDIR + 'spatial_qc_total_counts.png', dpi=150, bbox_inches='tight')
plt.close()
print("  spatial_qc_total_counts.png")

# 6b. Spatial distribution of mt%
fig, axes = plt.subplots(2, 4, figsize=(24, 12))
for idx, (s, adata) in enumerate(filtered_adatas.items()):
    ax = axes[idx // 4, idx % 4]
    coords = adata.obsm['spatial']
    sc_plot = ax.scatter(coords[:, 0], coords[:, 1],
                        c=adata.obs['pct_counts_mt'], cmap='RdYlGn_r',
                        s=8, alpha=0.8, vmin=0, vmax=25)
    ax.set_title(f'{s} ({adata.obs["condition"].iloc[0]})')
    ax.set_aspect('equal')
    ax.invert_yaxis()
    plt.colorbar(sc_plot, ax=ax, label='% MT', shrink=0.6)
fig.suptitle('GSE287278 Visium - Mitochondrial %', fontsize=16)
plt.tight_layout()
plt.savefig(FIGDIR + 'spatial_qc_mt_pct.png', dpi=150, bbox_inches='tight')
plt.close()
print("  spatial_qc_mt_pct.png")

# 6c. UMAP overview
fig, axes = plt.subplots(1, 3, figsize=(21, 6))
sc.pl.umap(adata_spatial, color='condition', ax=axes[0], show=False, title='Condition')
sc.pl.umap(adata_spatial, color='sample', ax=axes[1], show=False, title='Sample')
sc.pl.umap(adata_spatial, color='leiden_0.5', ax=axes[2], show=False, title='Leiden 0.5')
plt.tight_layout()
plt.savefig(FIGDIR + 'spatial_umap_overview.png', dpi=150, bbox_inches='tight')
plt.close()
print("  spatial_umap_overview.png")

# 6d. QC violin plots
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
sc.pl.violin(adata_spatial, keys='total_counts', groupby='sample', ax=axes[0],
             show=False, rotation=45)
sc.pl.violin(adata_spatial, keys='n_genes_by_counts', groupby='sample', ax=axes[1],
             show=False, rotation=45)
sc.pl.violin(adata_spatial, keys='pct_counts_mt', groupby='sample', ax=axes[2],
             show=False, rotation=45)
plt.tight_layout()
plt.savefig(FIGDIR + 'spatial_qc_violin.png', dpi=150, bbox_inches='tight')
plt.close()
print("  spatial_qc_violin.png")

# ============================================================
# 7. Save results
# ============================================================
print("\n--- Saving results ---")

# Save merged spatial data
out_path = 'Phase_output/phase4_spatial/spatial_merged_qc.h5ad'
adata_spatial.write_h5ad(out_path)
print(f"  Merged data: {out_path} ({adata_spatial.shape})")

# Save per-sample QC data (includes spatial coordinates needed for downstream spatial analysis)
for s, adata in filtered_adatas.items():
    sample_path = f'Phase_output/phase4_spatial/{s}_qc.h5ad'
    adata.write_h5ad(sample_path)

print(f"  Per-sample data: 8 h5ad files")

# Save QC statistics
qc_report = {
    'dataset': 'GSE287278_RIF_visium',
    'platform': '10x Visium',
    'n_samples': len(samples),
    'conditions': {'CTR': 4, 'RIF': 4},
    'qc_thresholds': {
        'min_counts': MIN_COUNTS,
        'min_genes': MIN_GENES,
        'max_mt_pct': MAX_MT_PCT,
        'min_cells_per_gene': 5
    },
    'total_spots_before': sum(s['n_spots'] for s in qc_stats.values()),
    'total_spots_after': sum(s.get('n_spots_after_qc', 0) for s in qc_stats.values()),
    'n_common_genes': int(adata_spatial.shape[1]),
    'n_hvg': int(n_hvg),
    'clustering': {
        f'leiden_{r}': int(adata_spatial.obs[f'leiden_{r}'].nunique())
        for r in [0.3, 0.5, 0.8, 1.0]
    },
    'per_sample': qc_stats
}

with open('Phase_output/phase4_spatial/spatial_qc_report.json', 'w') as f:
    json.dump(qc_report, f, indent=2, default=str)
print("  spatial_qc_report.json")

print("\n" + "=" * 60)
print("Phase 4 Step 1 complete!")
print(f"Total: {qc_report['total_spots_before']} -> {qc_report['total_spots_after']} spots")
print(f"Genes: {qc_report['n_common_genes']} (HVG: {qc_report['n_hvg']})")
print("=" * 60)
