#!/usr/bin/env python3
"""
Phase 2 Step 4-5: Normalization, HVG selection, PCA dimensionality reduction
Perform standard preprocessing pipeline on 3 QC-completed datasets
"""
import scanpy as sc
import anndata as ad
import numpy as np
import os
import json
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 1

OUTPUT = '/home/moog/test/tongxue/Phase_output/phase2_qc'

# PCA parameters
N_HVG = 3000
N_PCS = 50

preprocess_summary = {}

def preprocess_dataset(adata, dataset_name, n_hvg=N_HVG, n_pcs=N_PCS):
    """Standard preprocessing pipeline: normalize -> log1p -> HVG -> scale -> PCA"""
    print(f"\n--- {dataset_name}: Normalization + HVG + PCA ---")
    print(f"  Input: {adata.n_obs} cells x {adata.n_vars} genes")

    # Preserve raw counts
    adata.layers['counts'] = adata.X.copy()

    # Normalization (library size normalization + log1p)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Save log-normalized data
    adata.layers['log_normalized'] = adata.X.copy()

    # Highly variable gene selection
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor='seurat_v3',
                                 layer='counts', subset=False)
    n_hvg_found = adata.var['highly_variable'].sum()
    print(f"  HVG: {n_hvg_found} genes selected (target={n_hvg})")

    # Scale (HVG only)
    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    sc.pp.scale(adata_hvg, max_value=10)

    # PCA
    sc.tl.pca(adata_hvg, n_comps=min(n_pcs, adata_hvg.n_obs - 1, adata_hvg.n_vars - 1))
    variance_ratio = adata_hvg.uns['pca']['variance_ratio']

    # Determine effective number of PCs (elbow rule: cutoff where variance contribution < 2%)
    cumvar = np.cumsum(variance_ratio)
    n_effective_pcs = np.searchsorted(variance_ratio < 0.02, True)
    if n_effective_pcs < 10:
        n_effective_pcs = 30  # safe default value
    n_effective_pcs = min(n_effective_pcs, n_pcs)

    print(f"  PCA: {adata_hvg.obsm['X_pca'].shape[1]} PCs computed")
    print(f"  Cumulative variance of top 30 PCs: {cumvar[29]*100:.1f}%, effective PCs: ~{n_effective_pcs}")

    # Store PCA results back into original adata
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    adata.varm['PCs'] = np.zeros((adata.n_vars, adata_hvg.obsm['X_pca'].shape[1]))
    hvg_idx = np.where(adata.var['highly_variable'])[0]
    adata.varm['PCs'][hvg_idx] = adata_hvg.varm['PCs']
    adata.uns['pca'] = adata_hvg.uns['pca']

    # Compute neighbor graph and UMAP (for QC visualization)
    sc.pp.neighbors(adata, n_pcs=min(30, n_effective_pcs), use_rep='X_pca')
    sc.tl.umap(adata)
    print(f"  UMAP complete")

    # Initial Leiden clustering (for doublet and QC assessment)
    sc.tl.leiden(adata, resolution=0.5, key_added='leiden_0.5')
    n_clusters = len(adata.obs['leiden_0.5'].unique())
    print(f"  Leiden(res=0.5): {n_clusters} clusters")

    preprocess_summary[dataset_name] = {
        'n_cells': int(adata.n_obs),
        'n_genes': int(adata.n_vars),
        'n_hvg': int(n_hvg_found),
        'n_pcs': int(adata_hvg.obsm['X_pca'].shape[1]),
        'cumvar_30pc': round(float(cumvar[29]) * 100, 1),
        'effective_pcs': int(n_effective_pcs),
        'n_clusters_leiden05': int(n_clusters),
    }

    del adata_hvg
    return adata


# Process 3 datasets
datasets = [
    ('GSE111976', 'GSE111976_qc.h5ad'),
    ('E-MTAB-10287', 'E-MTAB-10287_qc.h5ad'),
    ('GSE260658', 'GSE260658_qc.h5ad'),
]

for name, filename in datasets:
    path = os.path.join(OUTPUT, filename)
    adata = ad.read_h5ad(path)
    adata = preprocess_dataset(adata, name)
    adata.write_h5ad(path)
    print(f"  Saved: {path}")
    del adata

# Save summary
summary_path = os.path.join(OUTPUT, 'preprocess_summary.json')
with open(summary_path, 'w') as f:
    json.dump(preprocess_summary, f, indent=2, ensure_ascii=False)

print("\n" + "=" * 70)
print("Preprocessing summary")
print("=" * 70)
print(f"{'Dataset':<20} {'Cells':>7} {'HVG':>5} {'PCs':>4} {'30PC var%':>9} {'Clusters':>8}")
print("-" * 60)
for name, s in preprocess_summary.items():
    print(f"{name:<20} {s['n_cells']:>7} {s['n_hvg']:>5} {s['n_pcs']:>4} {s['cumvar_30pc']:>8.1f}% {s['n_clusters_leiden05']:>8}")
