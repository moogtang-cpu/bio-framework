#!/usr/bin/env python
"""
Phase 4 Step 2-3: Spatial clustering, region annotation, and Cell2location deconvolution
- Spatial clustering + marker gene annotation of functional/basal layers
- Cell2location: spatial deconvolution using Phase 3 scRNA-seq reference (314,805 cells, 14 cell types)
- Spatial localization validation of stem cell subpopulations
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
import gc
import warnings
warnings.filterwarnings('ignore')

FIGDIR = 'Phase_output/phase4_spatial/figures/'
os.makedirs(FIGDIR, exist_ok=True)

# ============================================================
# 1. Load QC-processed spatial data
# ============================================================
print("=" * 60)
print("Phase 4 Step 2-3: Spatial Clustering & Deconvolution")
print("=" * 60)

adata_sp = sc.read_h5ad('Phase_output/phase4_spatial/spatial_merged_qc.h5ad')
print(f"Loaded merged spatial data: {adata_sp.shape}")

# Load per-sample data (contains spatial coordinates for spatial visualization)
samples = ['CTR1','CTR2','CTR3','CTR4','RIF1','RIF2','RIF3','RIF4']
sample_adatas = {}
for s in samples:
    sample_adatas[s] = sc.read_h5ad(f'Phase_output/phase4_spatial/{s}_qc.h5ad')

# ============================================================
# 2. Spatial clustering + region annotation with marker genes
# ============================================================
print("\n--- Spatial clustering and marker gene analysis ---")

# Endometrial layer marker genes
layer_markers = {
    # Basal layer / stem cell niche
    'Basalis_stem': ['SOX9', 'LGR5', 'PROM1', 'AXIN2', 'ALDH1A1', 'SSEA1'],
    # Functional layer
    'Functionalis': ['MMP7', 'MMP11', 'LEFTY2', 'SPP1', 'PRL'],
    # Glandular
    'Glandular': ['MUC1', 'FOXA2', 'PAX8', 'EPCAM', 'KRT8', 'KRT18'],
    # Stromal
    'Stromal': ['VIM', 'DCN', 'COL1A1', 'COL3A1', 'PDGFRB'],
    # Decidualized
    'Decidualized': ['IGFBP1', 'PRL', 'FOXO1', 'WNT4'],
    # Immune
    'Immune': ['PTPRC', 'CD68', 'CD3D', 'NKG7', 'GNLY'],
    # Vascular
    'Endothelial': ['PECAM1', 'VWF', 'CDH5', 'FLT1'],
    # Smooth muscle / myometrium
    'Smooth_muscle': ['ACTA2', 'MYH11', 'DES', 'CNN1'],
}

# Check marker gene availability in merged data
available_markers = {}
for category, genes in layer_markers.items():
    present = [g for g in genes if g in adata_sp.var_names]
    available_markers[category] = present
    print(f"  {category}: {len(present)}/{len(genes)} markers available ({', '.join(present)})")

# Dotplot: all marker genes x Leiden clusters
all_markers = [g for genes in available_markers.values() for g in genes]
# Deduplicate while preserving order
seen = set()
unique_markers = []
for g in all_markers:
    if g not in seen:
        seen.add(g)
        unique_markers.append(g)

# Use log_normalized layer for marker analysis
adata_sp.X = adata_sp.layers['log_normalized'].copy()

fig, ax = plt.subplots(figsize=(20, 8))
sc.pl.dotplot(adata_sp, var_names=unique_markers, groupby='leiden_0.5',
              show=False, ax=ax)
plt.tight_layout()
plt.savefig(FIGDIR + 'spatial_dotplot_markers.png', dpi=150, bbox_inches='tight')
plt.close()
print("  spatial_dotplot_markers.png")

# ============================================================
# 3. Spatial region annotation based on marker genes
# ============================================================
print("\n--- Spatial region annotation based on marker genes ---")

# Compute mean marker gene expression per cluster
cluster_col = 'leiden_0.5'
cluster_labels = adata_sp.obs[cluster_col].unique()

# Use log_normalized data for computation
from scipy.sparse import issparse
X = adata_sp.layers['log_normalized']

marker_scores = {}
for category, genes in available_markers.items():
    if not genes:
        continue
    gene_idx = [list(adata_sp.var_names).index(g) for g in genes]
    for cl in cluster_labels:
        mask = (adata_sp.obs[cluster_col] == cl).values
        if issparse(X):
            expr = np.array(X[mask][:, gene_idx].mean(axis=0)).flatten()
        else:
            expr = X[mask][:, gene_idx].mean(axis=0)
        if cl not in marker_scores:
            marker_scores[cl] = {}
        marker_scores[cl][category] = float(np.mean(expr))

# Print top-scoring region for each cluster
print(f"\n  Cluster -> Region annotation (Leiden {cluster_col}):")
cluster_annotations = {}
for cl in sorted(cluster_labels, key=lambda x: int(x)):
    scores = marker_scores[cl]
    # Exclude Immune and Endothelial from primary region assignment (they can coexist in any layer)
    tissue_scores = {k: v for k, v in scores.items()
                     if k not in ['Immune', 'Endothelial']}
    top_region = max(tissue_scores, key=tissue_scores.get)
    top_score = tissue_scores[top_region]

    # Determine basalis vs functionalis
    basalis_score = scores.get('Basalis_stem', 0)
    func_score = scores.get('Functionalis', 0)
    gland_score = scores.get('Glandular', 0)
    stroma_score = scores.get('Stromal', 0)
    decid_score = scores.get('Decidualized', 0)
    sm_score = scores.get('Smooth_muscle', 0)

    # Annotation logic
    n_cells = int((adata_sp.obs[cluster_col] == cl).sum())
    if basalis_score > 0.3 and basalis_score > func_score:
        annotation = 'Basalis_niche'
    elif sm_score > stroma_score and sm_score > gland_score:
        annotation = 'Myometrium'
    elif decid_score > 0.3 and decid_score > stroma_score * 0.5:
        annotation = 'Decidualized_stroma'
    elif func_score > 0.2 and func_score > basalis_score:
        annotation = 'Functionalis'
    elif gland_score > stroma_score:
        annotation = 'Glandular_epithelium'
    elif stroma_score > gland_score:
        annotation = 'Stroma'
    else:
        annotation = top_region

    cluster_annotations[cl] = annotation
    print(f"    Cluster {cl} ({n_cells} spots): {annotation} "
          f"[basalis={basalis_score:.2f}, func={func_score:.2f}, "
          f"gland={gland_score:.2f}, stroma={stroma_score:.2f}]")

# Add annotations to adata
adata_sp.obs['spatial_region'] = adata_sp.obs[cluster_col].map(cluster_annotations).astype('category')
region_counts = adata_sp.obs['spatial_region'].value_counts()
print(f"\n  Region distribution:")
for region, count in region_counts.items():
    print(f"    {region}: {count} spots ({count/adata_sp.shape[0]*100:.1f}%)")

# ============================================================
# 4. Spatial region visualization
# ============================================================
print("\n--- Spatial region visualization ---")

# 4a. Spatial region map for each sample
fig, axes = plt.subplots(2, 4, figsize=(28, 14))
region_colors = {
    'Basalis_niche': '#d62728',
    'Functionalis': '#2ca02c',
    'Glandular_epithelium': '#1f77b4',
    'Stroma': '#ff7f0e',
    'Decidualized_stroma': '#9467bd',
    'Myometrium': '#8c564b',
}

for idx, (s, adata) in enumerate(sample_adatas.items()):
    ax = axes[idx // 4, idx % 4]
    # Get region annotations corresponding to this sample
    sample_mask = adata_sp.obs['sample'] == s
    sample_regions = adata_sp.obs.loc[sample_mask, 'spatial_region']

    coords = adata.obsm['spatial']
    # Match indices
    common_idx = sample_regions.index.intersection(adata.obs_names)
    if len(common_idx) > 0:
        regions_aligned = sample_regions.loc[common_idx]
        coords_aligned = adata[common_idx].obsm['spatial']

        for region in regions_aligned.unique():
            rmask = regions_aligned == region
            color = region_colors.get(region, '#aaaaaa')
            ax.scatter(coords_aligned[rmask, 0], coords_aligned[rmask, 1],
                      c=color, s=8, alpha=0.8, label=region)

    ax.set_title(f'{s} ({adata.obs["condition"].iloc[0]})')
    ax.set_aspect('equal')
    ax.invert_yaxis()
    if idx == 0:
        ax.legend(fontsize=6, loc='upper left')
    ax.set_xticks([])
    ax.set_yticks([])

fig.suptitle('GSE287278 Visium - Spatial Region Annotation', fontsize=16)
plt.tight_layout()
plt.savefig(FIGDIR + 'spatial_region_annotation.png', dpi=150, bbox_inches='tight')
plt.close()
print("  spatial_region_annotation.png")

# 4b. Spatial expression of key marker genes (per sample)
stem_markers = ['SOX9', 'LGR5', 'PROM1']
stem_available = [g for g in stem_markers if g in adata_sp.var_names]

if stem_available:
    fig, axes = plt.subplots(len(stem_available), 4, figsize=(20, 5 * len(stem_available)))
    if len(stem_available) == 1:
        axes = axes.reshape(1, -1)

    # Show only CTR samples to display normal tissue spatial distribution
    ctr_samples = ['CTR1', 'CTR2', 'CTR3', 'CTR4']
    for gi, gene in enumerate(stem_available):
        for si, s in enumerate(ctr_samples):
            ax = axes[gi, si]
            adata_s = sample_adatas[s]
            if gene in adata_s.var_names:
                # Use log_normalized
                if 'log_normalized' in adata_s.layers:
                    gene_expr = np.array(adata_s.layers['log_normalized'][:, adata_s.var_names == gene].todense()).flatten()
                else:
                    gene_expr = np.array(adata_s.X[:, adata_s.var_names == gene].todense()).flatten()
                coords = adata_s.obsm['spatial']
                sc_plot = ax.scatter(coords[:, 0], coords[:, 1],
                                   c=gene_expr, cmap='Reds', s=8, alpha=0.8,
                                   vmin=0)
                plt.colorbar(sc_plot, ax=ax, shrink=0.6)
            ax.set_title(f'{s} - {gene}')
            ax.set_aspect('equal')
            ax.invert_yaxis()
            ax.set_xticks([])
            ax.set_yticks([])

    fig.suptitle('Stem Cell Markers - Spatial Expression (CTR samples)', fontsize=14)
    plt.tight_layout()
    plt.savefig(FIGDIR + 'spatial_stem_markers_CTR.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  spatial_stem_markers_CTR.png")

    # Stem cell markers in RIF samples
    fig, axes = plt.subplots(len(stem_available), 4, figsize=(20, 5 * len(stem_available)))
    if len(stem_available) == 1:
        axes = axes.reshape(1, -1)

    rif_samples = ['RIF1', 'RIF2', 'RIF3', 'RIF4']
    for gi, gene in enumerate(stem_available):
        for si, s in enumerate(rif_samples):
            ax = axes[gi, si]
            adata_s = sample_adatas[s]
            if gene in adata_s.var_names:
                if 'log_normalized' in adata_s.layers:
                    gene_expr = np.array(adata_s.layers['log_normalized'][:, adata_s.var_names == gene].todense()).flatten()
                else:
                    gene_expr = np.array(adata_s.X[:, adata_s.var_names == gene].todense()).flatten()
                coords = adata_s.obsm['spatial']
                sc_plot = ax.scatter(coords[:, 0], coords[:, 1],
                                   c=gene_expr, cmap='Reds', s=8, alpha=0.8,
                                   vmin=0)
                plt.colorbar(sc_plot, ax=ax, shrink=0.6)
            ax.set_title(f'{s} - {gene}')
            ax.set_aspect('equal')
            ax.invert_yaxis()
            ax.set_xticks([])
            ax.set_yticks([])

    fig.suptitle('Stem Cell Markers - Spatial Expression (RIF samples)', fontsize=14)
    plt.tight_layout()
    plt.savefig(FIGDIR + 'spatial_stem_markers_RIF.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  spatial_stem_markers_RIF.png")

# ============================================================
# 5. Cell2location deconvolution
# ============================================================
print("\n--- Cell2location deconvolution ---")
print("  Loading scRNA-seq reference atlas (Phase 3)...")

import cell2location
import scvi
import torch

# Load reference atlas — only need counts layer and cell type annotations
# Check available memory first
import psutil
mem = psutil.virtual_memory()
print(f"  Available memory: {mem.available / 1e9:.1f} GB")

# To avoid memory issues, subsample for training the NB regression model
# At most 5000 cells per cell type
ref_path = 'Phase_output/phase3_integration/integrated_harmony.h5ad'
print(f"  Reading reference: {ref_path}")

# Read only required columns
adata_ref = sc.read_h5ad(ref_path, backed='r')
print(f"  Reference atlas: {adata_ref.shape}")

# Get cell type information
obs_df = adata_ref.obs[['celltype_manual', 'sample', 'dataset']].copy()
cell_types = obs_df['celltype_manual'].unique()
print(f"  Cell types: {len(cell_types)}: {list(cell_types)}")

# Subsample by cell type (at most 5000 per type)
np.random.seed(42)
sample_idx = []
for ct in cell_types:
    ct_idx = obs_df.index[obs_df['celltype_manual'] == ct].tolist()
    if len(ct_idx) > 5000:
        ct_idx = np.random.choice(ct_idx, 5000, replace=False).tolist()
    sample_idx.extend(ct_idx)
print(f"  Subsampled: {len(sample_idx)} cells (from {adata_ref.shape[0]})")

# Read subsampled subset into memory (requires counts layer)
adata_ref_sub = adata_ref[sample_idx].to_memory()
del adata_ref
gc.collect()

# Use counts layer as X
if 'counts' in adata_ref_sub.layers:
    adata_ref_sub.X = adata_ref_sub.layers['counts'].copy()
    print("  Using counts layer as X")
else:
    print("  Warning: no counts layer, using current X")

# Cell2location preparation
# Filter genes: select genes common to scRNA and spatial data
common_genes = list(set(adata_ref_sub.var_names) & set(adata_sp.var_names))
print(f"  Common genes: {len(common_genes)}")

adata_ref_sub = adata_ref_sub[:, common_genes].copy()

# Prepare spatial data (use counts)
adata_vis = adata_sp[:, common_genes].copy()
adata_vis.X = adata_vis.layers['counts'].copy()

# Fix: ensure obs_names are unique (different samples can share the same Visium barcodes)
adata_vis.obs_names_make_unique()

# Cell2location requires unique var names
adata_ref_sub.var_names_make_unique()

print(f"  Reference subset: {adata_ref_sub.shape}")
print(f"  Spatial data: {adata_vis.shape}")

# Step 5a: Train NB regression model (estimate reference gene expression signatures)
print("\n  Training NB regression model...")

# Set up cell2location model
cell2location.models.RegressionModel.setup_anndata(
    adata_ref_sub,
    batch_key='sample',
    labels_key='celltype_manual'
)

# Create model
mod_ref = cell2location.models.RegressionModel(adata_ref_sub)

# Train (GPU accelerated)
mod_ref.train(
    max_epochs=100,
    batch_size=2500,
    train_size=1,
    lr=0.002,
    accelerator='gpu',
)

# Export reference signatures
adata_ref_sub = mod_ref.export_posterior(
    adata_ref_sub,
    sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
)

# Extract estimated gene expression signatures per cell type
if 'means_per_cluster_mu_fg' in adata_ref_sub.varm:
    inf_aver = adata_ref_sub.varm['means_per_cluster_mu_fg'].copy()
    inf_aver.columns = adata_ref_sub.uns['mod']['factor_names']
    print(f"  Gene signature matrix: {inf_aver.shape}")
else:
    print("  Using alternative method to extract signatures...")
    inf_aver = pd.DataFrame(
        adata_ref_sub.varm['means_per_cluster_mu_fg'],
        index=adata_ref_sub.var_names,
        columns=adata_ref_sub.uns['mod']['factor_names']
    )

print(f"  Cell type signatures: {list(inf_aver.columns)}")

del mod_ref
gc.collect()
torch.cuda.empty_cache()

# Step 5b: Train Cell2location spatial model
print("\n  Training Cell2location spatial model...")

cell2location.models.Cell2location.setup_anndata(
    adata_vis,
    batch_key='sample'
)

mod_spatial = cell2location.models.Cell2location(
    adata_vis,
    cell_state_df=inf_aver,
    N_cells_per_location=15,  # Visium spots typically contain 15-30 cells
    detection_alpha=20
)

mod_spatial.train(
    max_epochs=5000,
    batch_size=None,
    train_size=1,
    accelerator='gpu',
)

# Export posterior estimates
adata_vis = mod_spatial.export_posterior(
    adata_vis,
    sample_kwargs={'num_samples': 1000, 'batch_size': mod_spatial.adata.n_obs}
)

# Save deconvolution results
print("\n  Saving deconvolution results...")
adata_vis.write_h5ad('Phase_output/phase4_spatial/spatial_deconvolution.h5ad')

# Extract cell type proportions per spot
cell_type_cols = [c for c in adata_vis.obsm['q05_cell_abundance_w_sf'].columns]
deconv_df = adata_vis.obsm['q05_cell_abundance_w_sf'].copy()

# Normalize to proportions
deconv_props = deconv_df.div(deconv_df.sum(axis=1), axis=0)
print(f"  Deconvolution cell types: {list(cell_type_cols)}")

# Save deconvolution proportions
deconv_props.to_csv('Phase_output/phase4_spatial/deconvolution_proportions.csv')

# ============================================================
# 6. Deconvolution visualization
# ============================================================
print("\n--- Deconvolution visualization ---")

# 6a. Spatial distribution of key cell types (CTR vs RIF)
key_types = ['SOX9+LGR5+_stem', 'CD133+_progenitor', 'Decidualized_stromal',
             'Secretory_glandular', 'Macrophages', 'Lymphoid_immune']
# Show only existing types
key_types = [ct for ct in key_types if ct in deconv_props.columns]

if key_types:
    for s_name in samples:
        fig, axes = plt.subplots(1, len(key_types), figsize=(5*len(key_types), 5))
        if len(key_types) == 1:
            axes = [axes]

        s_mask = adata_vis.obs['sample'] == s_name
        s_data = adata_vis[s_mask]
        coords = s_data.obsm['spatial']

        for ki, ct in enumerate(key_types):
            ax = axes[ki]
            values = deconv_props.loc[s_mask, ct].values
            sc_plot = ax.scatter(coords[:, 0], coords[:, 1],
                               c=values, cmap='Reds', s=10, alpha=0.8)
            ax.set_title(ct.replace('_', '\n'), fontsize=8)
            ax.set_aspect('equal')
            ax.invert_yaxis()
            ax.set_xticks([])
            ax.set_yticks([])
            plt.colorbar(sc_plot, ax=ax, shrink=0.6)

        fig.suptitle(f'{s_name} ({adata_vis.obs.loc[s_mask, "condition"].iloc[0]}) - Cell Type Deconvolution',
                    fontsize=12)
        plt.tight_layout()
        plt.savefig(FIGDIR + f'deconv_{s_name}.png', dpi=150, bbox_inches='tight')
        plt.close()

    print(f"  deconv_*.png (8 files)")

# 6b. Stem cell proportion comparison: CTR vs RIF
stem_types = [ct for ct in ['SOX9+LGR5+_stem', 'CD133+_progenitor'] if ct in deconv_props.columns]
if stem_types:
    fig, axes = plt.subplots(1, len(stem_types), figsize=(6*len(stem_types), 5))
    if len(stem_types) == 1:
        axes = [axes]

    for si, st in enumerate(stem_types):
        ax = axes[si]
        ctr_vals = deconv_props.loc[adata_vis.obs['condition'] == 'CTR', st]
        rif_vals = deconv_props.loc[adata_vis.obs['condition'] == 'RIF', st]

        ax.boxplot([ctr_vals, rif_vals], labels=['CTR', 'RIF'])
        ax.set_title(st)
        ax.set_ylabel('Cell type proportion')

        # Mann-Whitney U test
        from scipy.stats import mannwhitneyu
        stat, pval = mannwhitneyu(ctr_vals, rif_vals, alternative='two-sided')
        ax.text(0.5, 0.95, f'p={pval:.2e}', transform=ax.transAxes,
               ha='center', va='top', fontsize=10)

    fig.suptitle('Stem/Progenitor Cell Abundance: CTR vs RIF', fontsize=14)
    plt.tight_layout()
    plt.savefig(FIGDIR + 'deconv_stem_CTR_vs_RIF.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  deconv_stem_CTR_vs_RIF.png")

# ============================================================
# 7. Stem cell spatial localization statistics
# ============================================================
print("\n--- Stem cell spatial localization analysis ---")

# Compute spatial distribution of stem cells per condition
if stem_types:
    stem_spatial_stats = {}
    for cond in ['CTR', 'RIF']:
        cond_mask = adata_vis.obs['condition'] == cond
        cond_deconv = deconv_props.loc[cond_mask]

        stats = {}
        for st in stem_types:
            vals = cond_deconv[st]
            # High-enrichment spots (proportion > median + 1 MAD)
            median_val = vals.median()
            mad_val = (vals - median_val).abs().median()
            enriched = vals > (median_val + mad_val)

            stats[st] = {
                'mean_proportion': float(vals.mean()),
                'median_proportion': float(vals.median()),
                'n_enriched_spots': int(enriched.sum()),
                'pct_enriched': float(enriched.sum() / len(vals) * 100),
                'max_proportion': float(vals.max()),
            }

        stem_spatial_stats[cond] = stats

    print(f"\n  {'Cell Type':<25} {'CTR mean':>10} {'RIF mean':>10} {'CTR enriched':>14} {'RIF enriched':>14}")
    print("  " + "-" * 75)
    for st in stem_types:
        ctr = stem_spatial_stats['CTR'][st]
        rif = stem_spatial_stats['RIF'][st]
        print(f"  {st:<25} {ctr['mean_proportion']:>10.4f} {rif['mean_proportion']:>10.4f} "
              f"{ctr['n_enriched_spots']:>6} ({ctr['pct_enriched']:.1f}%) "
              f"{rif['n_enriched_spots']:>6} ({rif['pct_enriched']:.1f}%)")

# ============================================================
# 8. Update merged data and save
# ============================================================
print("\n--- Saving final results ---")

# Add region annotations back to main data
adata_sp.obs['spatial_region'] = adata_sp.obs[cluster_col].map(cluster_annotations).astype('category')
adata_sp.write_h5ad('Phase_output/phase4_spatial/spatial_merged_qc.h5ad')

# Save comprehensive report
report = {
    'phase': 4,
    'dataset': 'GSE287278_RIF_visium',
    'n_spots': int(adata_sp.shape[0]),
    'n_samples': len(samples),
    'spatial_regions': {str(k): int(v) for k, v in region_counts.items()},
    'cluster_annotations': {str(k): v for k, v in cluster_annotations.items()},
    'deconvolution': {
        'method': 'Cell2location',
        'reference': 'Phase 3 integrated atlas (314,805 cells, 14 types)',
        'n_reference_sampled': len(sample_idx),
        'n_common_genes': len(common_genes),
        'cell_types': list(cell_type_cols),
    },
}

if stem_types:
    report['stem_cell_spatial'] = stem_spatial_stats

with open('Phase_output/phase4_spatial/phase4_report.json', 'w') as f:
    json.dump(report, f, indent=2, default=str)
print("  phase4_report.json")

del mod_spatial
gc.collect()
torch.cuda.empty_cache()

print("\n" + "=" * 60)
print("Phase 4 Step 2-3 complete!")
print(f"Spatial region annotation: {len(set(cluster_annotations.values()))} categories")
print(f"Cell2location deconvolution: {len(cell_type_cols)} cell types")
print("=" * 60)
