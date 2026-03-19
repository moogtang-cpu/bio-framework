#!/usr/bin/env python3
"""
Phase 3 Step 1: Multi-dataset batch integration (Harmony)
Integrate 7 endometrial tissue datasets to build a unified atlas (~783K cells)
Excluded: GSE216748 organoids (analyzed separately), HECA nuclei (snRNA-seq)
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import harmonypy as hm
import os
import json
import gc
import time
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 2
sc.settings.n_jobs = 16

OUTPUT = '/home/moog/test/tongxue/Phase_output/phase3_integration'
os.makedirs(OUTPUT, exist_ok=True)

t_start = time.time()

# =====================================================================
# 1. Dataset definitions and loading strategy
# =====================================================================
DATASETS = [
    # --- Reference atlas (already integrated, has counts layer) ---
    {
        'name': 'HECA_cells',
        'path': '/home/moog/test/tongxue/data/organized/HECA_reference/HECA_cells.h5ad',
        'has_counts': True,
        'counts_key': 'counts',
        'celltype_col': 'celltype',
        'lineage_col': 'lineage',
        'sample_col': 'sample',
        'condition': 'reference',
        'source': 'HECA',
        'extra_cols': {'stage': 'Binary Stage', 'pathology': 'Endometrial_pathology'},
    },
    # --- GSE215968 Asherman (already normalized, no counts layer) ---
    {
        'name': 'GSE215968_AS_vs_WOI',
        'path': '/home/moog/test/tongxue/data/organized/GSE215968_Asherman/GSE215968_AS_vs_WOI_Control.h5ad',
        'has_counts': False,
        'celltype_col': 'principal_cell_types',
        'sample_col': 'orig.ident',
        'condition_col': 'Group',
        'condition': 'AS_and_WOI',
        'source': 'GSE215968',
    },
    {
        'name': 'GSE215968_AS_pre_post',
        'path': '/home/moog/test/tongxue/data/organized/GSE215968_Asherman/GSE215968_AS_pre_vs_post.h5ad',
        'has_counts': False,
        'celltype_col': 'cell_type',
        'sample_col': 'orig.ident',
        'condition_col': 'Treatment_stage',
        'condition': 'AS_treatment',
        'source': 'GSE215968',
    },
    {
        'name': 'GSE215968_CD133',
        'path': '/home/moog/test/tongxue/data/organized/GSE215968_Asherman/GSE215968_CD133.h5ad',
        'has_counts': False,
        'celltype_col': 'cell_type',
        'sample_col': 'orig.ident',
        'condition': 'CD133_sorted',
        'source': 'GSE215968',
    },
    # --- Phase 2 QC datasets (have counts layer) ---
    {
        'name': 'GSE111976',
        'path': '/home/moog/test/tongxue/Phase_output/phase2_qc/GSE111976_qc.h5ad',
        'has_counts': True,
        'counts_key': 'counts',
        'celltype_col': None,
        'sample_col': 'sample_id',
        'condition': 'menstrual_cycle',
        'source': 'GSE111976',
    },
    {
        'name': 'E-MTAB-10287',
        'path': '/home/moog/test/tongxue/Phase_output/phase2_qc/E-MTAB-10287_qc.h5ad',
        'has_counts': True,
        'counts_key': 'counts',
        'celltype_col': None,
        'sample_col': 'sample',
        'condition': 'temporal',
        'source': 'E-MTAB-10287',
    },
    {
        'name': 'GSE260658',
        'path': '/home/moog/test/tongxue/Phase_output/phase2_qc/GSE260658_qc.h5ad',
        'has_counts': True,
        'counts_key': 'counts',
        'celltype_col': None,
        'sample_col': 'sample',
        'condition': 'uterus_atlas',
        'source': 'GSE260658',
        'extra_cols': {'tissue_layer': 'tissue_layer'},
    },
]

# =====================================================================
# 2. Load datasets one by one, standardize metadata, unify to log-normalized space
# =====================================================================
print("=" * 70)
print("Step 1: Load and preprocess each dataset")
print("=" * 70)

adatas = []
all_genes = []
integration_log = {}

for cfg in DATASETS:
    name = cfg['name']
    print(f"\n--- Loading {name} ---")
    t0 = time.time()

    adata = ad.read_h5ad(cfg['path'])
    print(f"  Raw: {adata.n_obs} cells x {adata.n_vars} genes")

    # --- Standardize obs metadata ---
    new_obs = pd.DataFrame(index=adata.obs_names)
    new_obs['dataset'] = name
    new_obs['source'] = cfg['source']
    new_obs['condition'] = cfg['condition']

    # Sample ID
    sample_col = cfg.get('sample_col')
    if sample_col and sample_col in adata.obs.columns:
        new_obs['sample'] = adata.obs[sample_col].astype(str).values
    else:
        new_obs['sample'] = name

    # Original cell type annotation
    ct_col = cfg.get('celltype_col')
    if ct_col and ct_col in adata.obs.columns:
        new_obs['original_celltype'] = adata.obs[ct_col].astype(str).values
    else:
        new_obs['original_celltype'] = 'unannotated'

    # Lineage annotation (if available)
    lin_col = cfg.get('lineage_col')
    if lin_col and lin_col in adata.obs.columns:
        new_obs['original_lineage'] = adata.obs[lin_col].astype(str).values
    else:
        new_obs['original_lineage'] = 'unknown'

    # Condition detail column
    cond_col = cfg.get('condition_col')
    if cond_col and cond_col in adata.obs.columns:
        new_obs['condition_detail'] = adata.obs[cond_col].astype(str).values
    else:
        new_obs['condition_detail'] = cfg['condition']

    # Extra columns
    for new_name, old_name in cfg.get('extra_cols', {}).items():
        if old_name in adata.obs.columns:
            new_obs[new_name] = adata.obs[old_name].astype(str).values

    # --- Process X matrix ---
    if cfg['has_counts']:
        counts_key = cfg.get('counts_key', 'counts')
        if counts_key in adata.layers:
            X_counts = adata.layers[counts_key]
        else:
            X_counts = adata.X

        # Ensure sparse float32
        if sp.issparse(X_counts):
            X_counts = X_counts.tocsr().astype(np.float32)
        else:
            X_counts = sp.csr_matrix(X_counts, dtype=np.float32)

        # Create new adata, re-normalize from counts
        adata_new = ad.AnnData(
            X=X_counts.copy(),
            obs=new_obs,
            var=pd.DataFrame(index=adata.var_names),
        )
        adata_new.layers['counts'] = X_counts
        sc.pp.normalize_total(adata_new, target_sum=1e4)
        sc.pp.log1p(adata_new)
        print(f"  Re-normalized from counts: normalize_total(1e4) + log1p")
    else:
        # X is already log-normalized
        X_norm = adata.X
        if sp.issparse(X_norm):
            X_norm = X_norm.tocsr().astype(np.float32)
        else:
            X_norm = sp.csr_matrix(X_norm, dtype=np.float32)

        adata_new = ad.AnnData(
            X=X_norm,
            obs=new_obs,
            var=pd.DataFrame(index=adata.var_names),
        )
        print(f"  X already normalized, using directly")

    adata_new.var_names_make_unique()
    all_genes.append(set(adata_new.var_names))

    # Ensure unique obs_names (add prefix)
    adata_new.obs_names = name + '_' + adata_new.obs_names.astype(str)

    integration_log[name] = {
        'n_cells': int(adata_new.n_obs),
        'n_genes': int(adata_new.n_vars),
        'has_counts': cfg['has_counts'],
        'condition': cfg['condition'],
        'source': cfg['source'],
    }

    adatas.append(adata_new)
    del adata, adata_new
    gc.collect()
    print(f"  Elapsed: {time.time()-t0:.1f}s")

# =====================================================================
# 3. Unify gene space (take intersection)
# =====================================================================
print(f"\n{'='*70}")
print("Step 2: Unify gene space")
print(f"{'='*70}")

common_genes = sorted(set.intersection(*all_genes))
print(f"Gene counts per dataset: {[a.n_vars for a in adatas]}")
print(f"Common gene intersection: {len(common_genes)}")

# Subset to common genes
for i in range(len(adatas)):
    adatas[i] = adatas[i][:, common_genes].copy()

# Verify that key marker genes are in the common gene set
key_markers = ['SOX9', 'LGR5', 'AXIN2', 'CD44', 'EPCAM', 'VIM', 'PECAM1',
               'PTPRC', 'MKI67', 'FOXJ1', 'MUC5B', 'KRT18', 'ACTA2',
               'PGR', 'ESR1', 'HOXA10', 'HOXA11', 'WNT7A', 'NOTCH1',
               'MMP7', 'PROM1', 'THY1', 'SUSD2', 'SSEA1']
found = [g for g in key_markers if g in common_genes]
missing = [g for g in key_markers if g not in common_genes]
print(f"Key marker genes: {len(found)}/{len(key_markers)} present in common gene set")
if missing:
    print(f"  Missing: {missing}")

# =====================================================================
# 4. Merge datasets
# =====================================================================
print(f"\n{'='*70}")
print("Step 3: Merge datasets")
print(f"{'='*70}")

adata = ad.concat(adatas, join='inner')
adata.obs_names_make_unique()
print(f"Merged: {adata.n_obs} cells x {adata.n_vars} genes")
print(f"Dataset distribution:\n{adata.obs['dataset'].value_counts().to_string()}")

del adatas
gc.collect()

# =====================================================================
# 5. HVG selection (batch-aware)
# =====================================================================
print(f"\n{'='*70}")
print("Step 4: Highly variable gene selection (batch-aware)")
print(f"{'='*70}")

# Use cell_ranger flavor adapted for log-normalized data, select per dataset batch
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    flavor='seurat',
    batch_key='dataset',
    subset=False,
)
n_hvg = adata.var['highly_variable'].sum()
print(f"Selected {n_hvg} highly variable genes (batch-aware across {adata.obs['dataset'].nunique()} datasets)")

# Save full-gene log-normalized matrix for downstream analysis
# Note: X is currently log-normalized
adata.raw = adata.copy()

# =====================================================================
# 6. Scale -> PCA (HVG only)
# =====================================================================
print(f"\n{'='*70}")
print("Step 5: Scale + PCA")
print(f"{'='*70}")

# Subset to HVG for PCA
adata_hvg = adata[:, adata.var['highly_variable']].copy()
print(f"HVG subset: {adata_hvg.n_obs} cells x {adata_hvg.n_vars} genes")

# Scale
sc.pp.scale(adata_hvg, max_value=10)
print("Scale complete")

# PCA
t_pca = time.time()
sc.tl.pca(adata_hvg, n_comps=50, svd_solver='arpack')
print(f"PCA complete: {adata_hvg.obsm['X_pca'].shape}, elapsed {time.time()-t_pca:.1f}s")

# Variance explained
cumvar = np.cumsum(adata_hvg.uns['pca']['variance_ratio'])
print(f"Cumulative variance of top 30 PCs: {cumvar[29]*100:.1f}%")
print(f"Cumulative variance of top 50 PCs: {cumvar[49]*100:.1f}%")

# Store PCA results in main adata
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata.uns['pca'] = adata_hvg.uns['pca']
del adata_hvg
gc.collect()

# =====================================================================
# 7. Harmony batch correction
# =====================================================================
print(f"\n{'='*70}")
print("Step 6: Harmony batch correction")
print(f"{'='*70}")

t_harmony = time.time()

# Use harmonypy
ho = hm.run_harmony(
    adata.obsm['X_pca'],
    adata.obs,
    'dataset',
    max_iter_harmony=20,
    verbose=True,
)

# harmonypy 0.2.0 PyTorch backend: Z_corr may be a GPU tensor
harmony_result = ho.Z_corr.T
if hasattr(harmony_result, 'cpu'):
    harmony_result = harmony_result.cpu().detach().numpy()
if harmony_result.ndim == 1:
    # If returned as 1D, try the result() method
    harmony_result = ho.result()
    if hasattr(harmony_result, 'cpu'):
        harmony_result = harmony_result.cpu().detach().numpy()
# Ensure correct shape (n_cells, n_PCs)
if harmony_result.shape[0] != adata.n_obs:
    harmony_result = harmony_result.T
adata.obsm['X_pca_harmony'] = np.array(harmony_result, dtype=np.float32)
print(f"Harmony complete: {adata.obsm['X_pca_harmony'].shape}, elapsed {time.time()-t_harmony:.1f}s")

# =====================================================================
# 8. Neighbor graph + UMAP
# =====================================================================
print(f"\n{'='*70}")
print("Step 7: Neighbors + UMAP")
print(f"{'='*70}")

t_umap = time.time()
sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_pcs=30, n_neighbors=30)
print(f"Neighbors complete, elapsed {time.time()-t_umap:.1f}s")

t_umap2 = time.time()
sc.tl.umap(adata, min_dist=0.3)
print(f"UMAP complete, elapsed {time.time()-t_umap2:.1f}s")

# =====================================================================
# 9. Leiden clustering (multi-resolution)
# =====================================================================
print(f"\n{'='*70}")
print("Step 8: Leiden clustering (multi-resolution)")
print(f"{'='*70}")

for res in [0.5, 1.0, 1.5, 2.0]:
    key = f'leiden_{res}'
    sc.tl.leiden(adata, resolution=res, key_added=key)
    n_clusters = adata.obs[key].nunique()
    print(f"  Leiden res={res}: {n_clusters} clusters")

# =====================================================================
# 10. Integration quality assessment
# =====================================================================
print(f"\n{'='*70}")
print("Step 9: Integration quality assessment")
print(f"{'='*70}")

# Dataset distribution across clusters (batch mixing)
ct = pd.crosstab(adata.obs['leiden_1.0'], adata.obs['dataset'])
# Compute dataset diversity per cluster (Shannon entropy)
from scipy.stats import entropy
cluster_entropy = ct.apply(lambda x: entropy(x / x.sum()), axis=1)
mean_entropy = cluster_entropy.mean()
max_possible_entropy = np.log(len(DATASETS))
normalized_entropy = mean_entropy / max_possible_entropy
print(f"Batch mixing score (normalized Shannon entropy): {normalized_entropy:.3f} (1.0=perfect mixing)")

# Simplified LISI metric: check whether each cluster is dominated by a single dataset
dominant_pct = ct.max(axis=1) / ct.sum(axis=1)
n_mixed = (dominant_pct < 0.9).sum()
n_total = len(dominant_pct)
print(f"Mixed clusters (no single dataset >90%): {n_mixed}/{n_total}")

# Dataset distribution across clusters
print("\nCell counts per dataset:")
print(adata.obs['dataset'].value_counts().to_string())

# =====================================================================
# 11. Save results
# =====================================================================
print(f"\n{'='*70}")
print("Step 10: Save results")
print(f"{'='*70}")

# Save full integrated object
out_path = os.path.join(OUTPUT, 'atlas_integrated.h5ad')
adata.write_h5ad(out_path)
print(f"Saved: {out_path}")
print(f"File size: {os.path.getsize(out_path)/1e9:.1f} GB")

# Save integration log
integration_log['total'] = {
    'n_cells': int(adata.n_obs),
    'n_common_genes': int(adata.n_vars),
    'n_hvg': int(n_hvg),
    'n_pcs': 50,
    'cumvar_30pc': round(float(cumvar[29] * 100), 1),
    'cumvar_50pc': round(float(cumvar[49] * 100), 1),
    'harmony_iters': ho.max_iter,
    'batch_mixing_entropy': round(float(normalized_entropy), 3),
    'n_mixed_clusters': int(n_mixed),
    'n_total_clusters': int(n_total),
    'leiden_cluster_counts': {
        f'res_{res}': int(adata.obs[f'leiden_{res}'].nunique())
        for res in [0.5, 1.0, 1.5, 2.0]
    },
    'key_markers_found': found,
    'key_markers_missing': missing,
}

log_path = os.path.join(OUTPUT, 'integration_log.json')
with open(log_path, 'w') as f:
    json.dump(integration_log, f, indent=2, ensure_ascii=False)
print(f"Log saved: {log_path}")

total_time = time.time() - t_start
print(f"\n{'='*70}")
print(f"Phase 3 Step 1 complete! Total elapsed: {total_time/60:.1f} minutes")
print(f"Integrated atlas: {adata.n_obs} cells x {adata.n_vars} genes")
print(f"{'='*70}")
