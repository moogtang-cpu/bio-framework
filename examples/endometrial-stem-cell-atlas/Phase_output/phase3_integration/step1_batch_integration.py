#!/usr/bin/env python3
"""
Phase 3 Step 1: Multi-dataset batch integration (Harmony)
Project: Construction of human endometrial basal layer stem cell atlas and discovery of
         drug intervention targets for regenerative disorders

Integration strategy:
- GSE215968_AS_vs_WOI (106,400) — Asherman's syndrome vs WOI control
- GSE215968_CD133 (69,701) — CD133+ enriched stem/progenitor cells
- GSE111976_qc (61,503) — normal menstrual cycle baseline
- E-MTAB-10287_qc (49,550) — endometrial temporal dynamics
- GSE260658_qc (endometrium + decidua subset, ~28K) — endometrium/decidua reference

Excluded:
- GSE215968_AS_pre_vs_post: contains post-treatment samples, introduces confounders, reserved for later validation
- GSE216748_organoid: in vitro culture data, reserved for Phase 8 target validation
- HECA: reference atlas, used for Phase 3 Step 3 label transfer
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import harmonypy as hm
from scipy.sparse import issparse, csr_matrix
import gc
import json
import os
import time

sc.settings.verbosity = 1
sc.settings.n_jobs = 16

BASE = "/home/moog/test/tongxue"
OUT_DIR = f"{BASE}/Phase_output/phase3_integration"
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print("Phase 3 Step 1: Multi-dataset batch integration")
print("=" * 60)

# ================================================================
# 1. Load datasets and retrieve raw counts
# ================================================================
print("\n[1/8] Loading datasets...", flush=True)
t0 = time.time()

datasets = {}

# --- GSE215968_AS_vs_WOI: back-calculate raw counts from log-normalized ---
print("  Loading GSE215968_AS_vs_WOI (back-calculating raw counts)...", flush=True)
adata = ad.read_h5ad(f"{BASE}/data/organized/GSE215968_Asherman/GSE215968_AS_vs_WOI_Control.h5ad")
total_counts = adata.obs['nCount_RNA'].values.astype(np.float64)
# Seurat LogNormalize: log1p(count/total*10000), back-calculate
if issparse(adata.X):
    x_dense = adata.X.toarray()
else:
    x_dense = adata.X.copy()
raw_counts = np.round(np.expm1(x_dense) * total_counts[:, None] / 10000).astype(np.float32)
raw_counts[raw_counts < 0] = 0
adata_new = ad.AnnData(
    X=csr_matrix(raw_counts),
    obs=pd.DataFrame({
        'dataset': 'GSE215968_AS',
        'sample': adata.obs['orig.ident'].values,
        'condition': adata.obs['Group'].values,
        'original_celltype': adata.obs['cell_subtypes'].values,
    }, index=adata.obs_names),
    var=pd.DataFrame(index=adata.var_names)
)
datasets['GSE215968_AS'] = adata_new
del adata, x_dense, raw_counts; gc.collect()
print(f"    {datasets['GSE215968_AS'].shape}", flush=True)

# --- GSE215968_CD133: back-calculate as well ---
print("  Loading GSE215968_CD133 (back-calculating raw counts)...", flush=True)
adata = ad.read_h5ad(f"{BASE}/data/organized/GSE215968_Asherman/GSE215968_CD133.h5ad")
total_counts = adata.obs['nCount_RNA'].values.astype(np.float64)
if issparse(adata.X):
    x_dense = adata.X.toarray()
else:
    x_dense = adata.X.copy()
raw_counts = np.round(np.expm1(x_dense) * total_counts[:, None] / 10000).astype(np.float32)
raw_counts[raw_counts < 0] = 0
adata_new = ad.AnnData(
    X=csr_matrix(raw_counts),
    obs=pd.DataFrame({
        'dataset': 'GSE215968_CD133',
        'sample': adata.obs['orig.ident'].values,
        'condition': 'AS_CD133',  # all AS patient CD133+ sorted cells
        'original_celltype': adata.obs['cell_type'].astype(str).values,
    }, index=adata.obs_names),
    var=pd.DataFrame(index=adata.var_names)
)
datasets['GSE215968_CD133'] = adata_new
del adata, x_dense, raw_counts; gc.collect()
print(f"    {datasets['GSE215968_CD133'].shape}", flush=True)

# --- GSE111976_qc: already has counts layer ---
print("  Loading GSE111976_qc...", flush=True)
adata = ad.read_h5ad(f"{BASE}/Phase_output/phase2_qc/GSE111976_qc.h5ad")
adata_new = ad.AnnData(
    X=adata.layers['counts'].copy(),
    obs=pd.DataFrame({
        'dataset': 'GSE111976',
        'sample': adata.obs['sample_id'].values if 'sample_id' in adata.obs else 'GSE111976',
        'condition': 'Normal',
        'original_celltype': 'unassigned',
    }, index=adata.obs_names),
    var=pd.DataFrame(index=adata.var_names)
)
datasets['GSE111976'] = adata_new
del adata; gc.collect()
print(f"    {datasets['GSE111976'].shape}", flush=True)

# --- E-MTAB-10287_qc: already has counts layer ---
print("  Loading E-MTAB-10287_qc...", flush=True)
adata = ad.read_h5ad(f"{BASE}/Phase_output/phase2_qc/E-MTAB-10287_qc.h5ad")
adata_new = ad.AnnData(
    X=adata.layers['counts'].copy(),
    obs=pd.DataFrame({
        'dataset': 'E-MTAB-10287',
        'sample': adata.obs['sample'].values,
        'condition': 'Normal',
        'original_celltype': 'unassigned',
    }, index=adata.obs_names),
    var=pd.DataFrame(index=adata.var_names)
)
datasets['E-MTAB-10287'] = adata_new
del adata; gc.collect()
print(f"    {datasets['E-MTAB-10287'].shape}", flush=True)

# --- GSE260658_qc: subset (endometrium + decidua, exclude myometrium) ---
print("  Loading GSE260658_qc (endometrium + decidua subset)...", flush=True)
adata = ad.read_h5ad(f"{BASE}/Phase_output/phase2_qc/GSE260658_qc.h5ad")
mask = adata.obs['tissue_layer'].isin(['endometrium', 'decidua'])
adata = adata[mask].copy()
adata_new = ad.AnnData(
    X=adata.layers['counts'].copy(),
    obs=pd.DataFrame({
        'dataset': 'GSE260658',
        'sample': adata.obs['sample'].values,
        'condition': adata.obs['tissue_layer'].values,
        'original_celltype': 'unassigned',
    }, index=adata.obs_names),
    var=pd.DataFrame(index=adata.var_names)
)
datasets['GSE260658'] = adata_new
del adata; gc.collect()
print(f"    {datasets['GSE260658'].shape}", flush=True)

t_load = time.time() - t0
print(f"  Loading complete ({t_load:.0f}s)", flush=True)

# ================================================================
# 2. Find common genes and merge
# ================================================================
print("\n[2/8] Gene intersection and merging...", flush=True)
common_genes = set(datasets['GSE215968_AS'].var_names)
for name, adata in datasets.items():
    common_genes &= set(adata.var_names)
common_genes = sorted(common_genes)
print(f"  Common genes: {len(common_genes)}", flush=True)

# Subset to common genes
for name in datasets:
    datasets[name] = datasets[name][:, common_genes].copy()

# Merge (obs already contains dataset column, do not use keys/label parameters)
adata_list = list(datasets.values())
adata_merged = ad.concat(adata_list, join='inner', index_unique='-')

# Ensure X is a sparse matrix
if not issparse(adata_merged.X):
    adata_merged.X = csr_matrix(adata_merged.X)

# Save raw counts to layers
adata_merged.layers['counts'] = adata_merged.X.copy()

total_cells = adata_merged.shape[0]
print(f"  After merge: {adata_merged.shape[0]} cells x {adata_merged.shape[1]} genes", flush=True)
print(f"  Dataset distribution:", flush=True)
for ds, cnt in adata_merged.obs['dataset'].value_counts().items():
    print(f"    {ds}: {cnt:,}", flush=True)

del datasets; gc.collect()

# ================================================================
# 3. Normalization and HVG selection
# ================================================================
print("\n[3/8] Normalization and HVG selection...", flush=True)

# log1p(10K) normalization
sc.pp.normalize_total(adata_merged, target_sum=1e4)
sc.pp.log1p(adata_merged)
adata_merged.layers['log_normalized'] = adata_merged.X.copy()

# Dataset-aware HVG selection (compute HVGs per dataset independently, take union of top-ranked)
sc.pp.highly_variable_genes(
    adata_merged,
    n_top_genes=3000,
    batch_key='dataset',
    flavor='seurat_v3',
    layer='counts',
    subset=False,
)
n_hvg = adata_merged.var['highly_variable'].sum()
print(f"  HVG: {n_hvg}", flush=True)

# ================================================================
# 4. PCA (optimized: scale only on HVG subset, use randomized SVD)
# ================================================================
print("\n[4/8] PCA dimensionality reduction...", flush=True)

# Extract HVG subset for PCA (avoid densifying entire genome)
adata_hvg = adata_merged[:, adata_merged.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=50, svd_solver='randomized', random_state=42)

# Write PCA results back to main object
adata_merged.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()
adata_merged.uns['pca'] = adata_hvg.uns['pca'].copy()
adata_merged.varm['PCs'] = np.zeros((adata_merged.shape[1], 50))
hvg_idx = np.where(adata_merged.var['highly_variable'])[0]
adata_merged.varm['PCs'][hvg_idx] = adata_hvg.varm['PCs']

del adata_hvg; gc.collect()

print(f"  PCA complete: 50 components (3000 HVGs)", flush=True)
print(f"  Variance explained: PC1={adata_merged.uns['pca']['variance_ratio'][0]:.3f}, "
      f"PC50={adata_merged.uns['pca']['variance_ratio'][49]:.4f}", flush=True)

# ================================================================
# 5. Harmony batch integration
# ================================================================
print("\n[5/8] Harmony batch integration...", flush=True)
t_harmony = time.time()

# Use harmonypy
ho = hm.run_harmony(
    adata_merged.obsm['X_pca'],
    adata_merged.obs,
    'dataset',
    max_iter_harmony=30,
    random_state=42,
)
# harmonypy returns a Harmony object; Z_corr shape=(n_cells, n_pcs)
Z_corr = ho.Z_corr
if hasattr(Z_corr, 'numpy'):
    Z_corr = Z_corr.cpu().numpy()  # PyTorch tensor -> numpy
if Z_corr.shape[0] != adata_merged.shape[0]:
    Z_corr = Z_corr.T  # transpose if shape is reversed
adata_merged.obsm['X_pca_harmony'] = Z_corr

t_harmony = time.time() - t_harmony
print(f"  Harmony complete ({t_harmony:.0f}s)", flush=True)

# ================================================================
# 6. Neighbor graph and UMAP
# ================================================================
print("\n[6/8] Neighbor graph and UMAP...", flush=True)

# Use Harmony-corrected PCA
sc.pp.neighbors(adata_merged, use_rep='X_pca_harmony', n_pcs=30, n_neighbors=30)
sc.tl.umap(adata_merged, min_dist=0.3, random_state=42)
print("  UMAP complete", flush=True)

# ================================================================
# 7. Multi-resolution Leiden clustering
# ================================================================
print("\n[7/8] Leiden clustering (multi-resolution)...", flush=True)

for res in [0.5, 0.8, 1.0, 1.5]:
    key = f'leiden_{res}'
    sc.tl.leiden(adata_merged, resolution=res, key_added=key, random_state=42)
    n_clusters = adata_merged.obs[key].nunique()
    print(f"  resolution={res}: {n_clusters} clusters", flush=True)

# ================================================================
# 8. Save and quality assessment
# ================================================================
print("\n[8/8] Saving results...", flush=True)

# Restore X to log-normalized (not scaled)
adata_merged.X = adata_merged.layers['log_normalized'].copy()

# Save
out_path = f"{OUT_DIR}/integrated_harmony.h5ad"
adata_merged.write_h5ad(out_path)
print(f"  Saved: {out_path}", flush=True)

# Integration quality statistics
stats = {
    'total_cells': int(adata_merged.shape[0]),
    'total_genes': int(adata_merged.shape[1]),
    'common_genes': len(common_genes),
    'n_hvg': int(n_hvg),
    'n_pcs': 50,
    'harmony_pcs_used': 30,
    'datasets': {ds: int(cnt) for ds, cnt in adata_merged.obs['dataset'].value_counts().items()},
    'clustering': {
        f'leiden_{res}': int(adata_merged.obs[f'leiden_{res}'].nunique())
        for res in [0.5, 0.8, 1.0, 1.5]
    },
    'conditions': {c: int(cnt) for c, cnt in adata_merged.obs['condition'].value_counts().items()},
    'load_time_s': round(t_load),
    'harmony_time_s': round(t_harmony),
}

with open(f"{OUT_DIR}/integration_stats.json", 'w') as f:
    json.dump(stats, f, indent=2, ensure_ascii=False)

print(f"\n{'='*60}")
print(f"Integration complete!")
print(f"  Total cells: {stats['total_cells']:,}")
print(f"  Common genes: {stats['common_genes']:,}")
print(f"  HVG: {stats['n_hvg']}")
print(f"  Leiden 0.5: {stats['clustering']['leiden_0.5']} clusters")
print(f"  Leiden 1.0: {stats['clustering']['leiden_1.0']} clusters")
print(f"{'='*60}")
