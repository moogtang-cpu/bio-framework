#!/usr/bin/env python3
"""
Phase 3 Step 1 continued: Resume from PCA checkpoint, complete Harmony -> UMAP -> Leiden
(Re-run second half after fixing Z_corr dimension issue)
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
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

# No checkpoint file available, need to rebuild integration object
# Optimization: start directly from the merge step

t_start = time.time()

print("=" * 70)
print("Reload datasets (optimized: load and integrate directly after unification)")
print("=" * 70)

import scipy.sparse as sp

DATASETS = [
    ('HECA_cells', '/home/moog/test/tongxue/data/organized/HECA_reference/HECA_cells.h5ad',
     True, 'counts', 'celltype', 'lineage', 'sample', 'reference', 'HECA',
     {'stage': 'Binary Stage', 'pathology': 'Endometrial_pathology'}),
    ('GSE215968_AS_vs_WOI', '/home/moog/test/tongxue/data/organized/GSE215968_Asherman/GSE215968_AS_vs_WOI_Control.h5ad',
     False, None, 'principal_cell_types', None, 'orig.ident', 'AS_and_WOI', 'GSE215968',
     {'condition_detail': 'Group'}),
    ('GSE215968_AS_pre_post', '/home/moog/test/tongxue/data/organized/GSE215968_Asherman/GSE215968_AS_pre_vs_post.h5ad',
     False, None, 'cell_type', None, 'orig.ident', 'AS_treatment', 'GSE215968',
     {'condition_detail': 'Treatment_stage'}),
    ('GSE215968_CD133', '/home/moog/test/tongxue/data/organized/GSE215968_Asherman/GSE215968_CD133.h5ad',
     False, None, 'cell_type', None, 'orig.ident', 'CD133_sorted', 'GSE215968', {}),
    ('GSE111976', '/home/moog/test/tongxue/Phase_output/phase2_qc/GSE111976_qc.h5ad',
     True, 'counts', None, None, 'sample_id', 'menstrual_cycle', 'GSE111976', {}),
    ('E-MTAB-10287', '/home/moog/test/tongxue/Phase_output/phase2_qc/E-MTAB-10287_qc.h5ad',
     True, 'counts', None, None, 'sample', 'temporal', 'E-MTAB-10287', {}),
    ('GSE260658', '/home/moog/test/tongxue/Phase_output/phase2_qc/GSE260658_qc.h5ad',
     True, 'counts', None, None, 'sample', 'uterus_atlas', 'GSE260658',
     {'tissue_layer': 'tissue_layer'}),
]

adatas = []
all_genes = []
integration_log = {}

for name, path, has_counts, ck, ct_col, lin_col, samp_col, condition, source, extra in DATASETS:
    print(f"\n--- {name} ---")
    t0 = time.time()
    adata = ad.read_h5ad(path)

    # Standardize obs
    new_obs = pd.DataFrame(index=adata.obs_names)
    new_obs['dataset'] = name
    new_obs['source'] = source
    new_obs['condition'] = condition
    new_obs['sample'] = adata.obs[samp_col].astype(str).values if samp_col and samp_col in adata.obs.columns else name
    new_obs['original_celltype'] = adata.obs[ct_col].astype(str).values if ct_col and ct_col in adata.obs.columns else 'unannotated'
    new_obs['original_lineage'] = adata.obs[lin_col].astype(str).values if lin_col and lin_col in adata.obs.columns else 'unknown'
    for new_n, old_n in extra.items():
        if old_n in adata.obs.columns:
            new_obs[new_n] = adata.obs[old_n].astype(str).values

    if has_counts:
        X = adata.layers[ck] if ck and ck in adata.layers else adata.X
        if sp.issparse(X):
            X = X.tocsr().astype(np.float32)
        else:
            X = sp.csr_matrix(X, dtype=np.float32)
        adata_new = ad.AnnData(X=X.copy(), obs=new_obs, var=pd.DataFrame(index=adata.var_names))
        adata_new.layers['counts'] = X
        sc.pp.normalize_total(adata_new, target_sum=1e4)
        sc.pp.log1p(adata_new)
    else:
        X = adata.X
        if sp.issparse(X):
            X = X.tocsr().astype(np.float32)
        else:
            X = sp.csr_matrix(X, dtype=np.float32)
        adata_new = ad.AnnData(X=X, obs=new_obs, var=pd.DataFrame(index=adata.var_names))

    adata_new.var_names_make_unique()
    all_genes.append(set(adata_new.var_names))
    adata_new.obs_names = name + '_' + adata_new.obs_names.astype(str)
    integration_log[name] = {
        'n_cells': int(adata_new.n_obs),
        'n_genes': int(adata_new.n_vars),
        'has_counts': has_counts,
        'condition': condition,
        'source': source,
    }
    adatas.append(adata_new)
    del adata
    gc.collect()
    print(f"  {adata_new.n_obs} cells, elapsed {time.time()-t0:.1f}s")

# Common genes
common_genes = sorted(set.intersection(*all_genes))
print(f"\nCommon genes: {len(common_genes)}")
for i in range(len(adatas)):
    adatas[i] = adatas[i][:, common_genes].copy()

# Key markers
key_markers = ['SOX9', 'LGR5', 'AXIN2', 'CD44', 'EPCAM', 'VIM', 'PECAM1',
               'PTPRC', 'MKI67', 'FOXJ1', 'MUC5B', 'KRT18', 'ACTA2',
               'PGR', 'ESR1', 'HOXA10', 'HOXA11', 'WNT7A', 'NOTCH1',
               'MMP7', 'PROM1', 'THY1', 'SUSD2']
found = [g for g in key_markers if g in common_genes]
missing = [g for g in key_markers if g not in common_genes]
print(f"Key markers: {len(found)}/{len(key_markers)} present (missing: {missing})")

# Merge
adata = ad.concat(adatas, join='inner')
adata.obs_names_make_unique()
print(f"\nMerged: {adata.n_obs} cells x {adata.n_vars} genes")
print(adata.obs['dataset'].value_counts().to_string())
del adatas
gc.collect()

# HVG
print("\n--- HVG selection ---")
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='seurat', batch_key='dataset', subset=False)
n_hvg = adata.var['highly_variable'].sum()
print(f"HVG: {n_hvg}")

# Save raw
adata.raw = adata.copy()

# Scale + PCA
print("\n--- Scale + PCA ---")
adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)
t_pca = time.time()
sc.tl.pca(adata_hvg, n_comps=50, svd_solver='arpack')
cumvar = np.cumsum(adata_hvg.uns['pca']['variance_ratio'])
print(f"PCA: {adata_hvg.obsm['X_pca'].shape}, 30PC var={cumvar[29]*100:.1f}%, 50PC={cumvar[49]*100:.1f}%, elapsed {time.time()-t_pca:.1f}s")
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata.uns['pca'] = adata_hvg.uns['pca']
del adata_hvg
gc.collect()

# Harmony
print("\n--- Harmony ---")
t_harmony = time.time()

# Use scanpy built-in interface (more stable)
sc.external.pp.harmony_integrate(
    adata,
    key='dataset',
    basis='X_pca',
    adjusted_basis='X_pca_harmony',
    max_iter_harmony=20,
    verbose=True,
)
print(f"Harmony: {adata.obsm['X_pca_harmony'].shape}, elapsed {time.time()-t_harmony:.1f}s")

# Neighbors + UMAP
print("\n--- Neighbors + UMAP ---")
t_nn = time.time()
sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_pcs=30, n_neighbors=30)
print(f"Neighbors: elapsed {time.time()-t_nn:.1f}s")

t_umap = time.time()
sc.tl.umap(adata, min_dist=0.3)
print(f"UMAP: elapsed {time.time()-t_umap:.1f}s")

# Leiden multi-resolution
print("\n--- Leiden clustering ---")
for res in [0.5, 1.0, 1.5, 2.0]:
    key = f'leiden_{res}'
    sc.tl.leiden(adata, resolution=res, key_added=key)
    print(f"  res={res}: {adata.obs[key].nunique()} clusters")

# Integration quality assessment
print("\n--- Integration quality assessment ---")
ct = pd.crosstab(adata.obs['leiden_1.0'], adata.obs['dataset'])
from scipy.stats import entropy
cluster_entropy = ct.apply(lambda x: entropy(x / x.sum()), axis=1)
mean_entropy = cluster_entropy.mean()
max_entropy = np.log(7)  # 7 datasets
norm_entropy = mean_entropy / max_entropy
print(f"Batch mixing score (normalized Shannon entropy): {norm_entropy:.3f}")

dominant_pct = ct.max(axis=1) / ct.sum(axis=1)
n_mixed = (dominant_pct < 0.9).sum()
n_total = len(dominant_pct)
print(f"Mixed clusters (no single dataset >90%): {n_mixed}/{n_total}")

# Save
print("\n--- Save ---")
out_path = os.path.join(OUTPUT, 'atlas_integrated.h5ad')
adata.write_h5ad(out_path)
print(f"Saved: {out_path}")
print(f"Size: {os.path.getsize(out_path)/1e9:.1f} GB")

# Log
integration_log['total'] = {
    'n_cells': int(adata.n_obs),
    'n_common_genes': int(adata.n_vars),
    'n_hvg': int(n_hvg),
    'n_pcs': 50,
    'cumvar_30pc': round(float(cumvar[29] * 100), 1),
    'cumvar_50pc': round(float(cumvar[49] * 100), 1),
    'batch_mixing_entropy': round(float(norm_entropy), 3),
    'n_mixed_clusters': int(n_mixed),
    'n_total_clusters': int(n_total),
    'leiden_clusters': {f'res_{r}': int(adata.obs[f'leiden_{r}'].nunique()) for r in [0.5, 1.0, 1.5, 2.0]},
    'key_markers_found': found,
    'key_markers_missing': missing,
}
with open(os.path.join(OUTPUT, 'integration_log.json'), 'w') as f:
    json.dump(integration_log, f, indent=2, ensure_ascii=False)

total_time = time.time() - t_start
print(f"\n{'='*70}")
print(f"Phase 3 Step 1 complete! Total elapsed: {total_time/60:.1f} minutes")
print(f"Integrated atlas: {adata.n_obs} cells x {adata.n_vars} genes")
print(f"{'='*70}")
