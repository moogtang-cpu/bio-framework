#!/usr/bin/env python3
"""
Phase 2 Step 1-2: QC metric calculation and cell filtering
Perform full QC on raw datasets (GSE111976, E-MTAB-10287, GSE260658)
Validate QC metrics on pre-processed datasets (HECA, GSE215968, GSE216748)
"""
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import os
import json
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 1

ORGANIZED = '/home/moog/test/tongxue/data/organized'
OUTPUT = '/home/moog/test/tongxue/Phase_output/phase2_qc'
os.makedirs(OUTPUT, exist_ok=True)

# QC parameters
QC_PARAMS = {
    'mt_percent': 20,
    'min_genes': 200,
    'max_genes': 8000,
    'min_cells': 3,
    'min_counts': 500,
}

qc_summary = {}

def compute_qc_metrics(adata, dataset_name):
    """Compute QC metrics"""
    # Ensure gene names are strings
    adata.var_names_make_unique()

    # Mitochondrial genes
    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
    # Ribosomal genes
    adata.var['ribo'] = adata.var_names.str.upper().str.startswith(('RPS', 'RPL'))
    # Hemoglobin genes
    adata.var['hb'] = adata.var_names.str.upper().str.match('^HB[^P]')

    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'],
                                percent_top=None, log1p=False, inplace=True)
    return adata

def filter_cells(adata, dataset_name, params=QC_PARAMS):
    """Filter low-quality cells"""
    n_before = adata.n_obs

    # Filter genes (expressed in at least min_cells cells)
    sc.pp.filter_genes(adata, min_cells=params['min_cells'])

    # Filter cells
    sc.pp.filter_cells(adata, min_genes=params['min_genes'])
    sc.pp.filter_cells(adata, min_counts=params['min_counts'])

    # Filter cells with high gene count (potential multiplets)
    adata = adata[adata.obs['n_genes_by_counts'] < params['max_genes']].copy()

    # Filter high mitochondrial percentage
    adata = adata[adata.obs['pct_counts_mt'] < params['mt_percent']].copy()

    n_after = adata.n_obs
    pct_kept = 100 * n_after / n_before if n_before > 0 else 0

    print(f"  Filtered: {n_before} -> {n_after} cells ({pct_kept:.1f}% retained)")
    return adata, n_before, n_after

def qc_stats(adata, dataset_name, is_filtered=False):
    """Collect QC statistics"""
    stats = {
        'dataset': dataset_name,
        'n_cells': int(adata.n_obs),
        'n_genes': int(adata.n_vars),
        'is_filtered': is_filtered,
    }

    if 'pct_counts_mt' in adata.obs.columns:
        stats['mt_pct_median'] = float(adata.obs['pct_counts_mt'].median())
        stats['mt_pct_mean'] = float(adata.obs['pct_counts_mt'].mean())
        stats['mt_pct_95'] = float(adata.obs['pct_counts_mt'].quantile(0.95))
    elif 'percent_mito' in adata.obs.columns:
        stats['mt_pct_median'] = float(adata.obs['percent_mito'].median())
        stats['mt_pct_mean'] = float(adata.obs['percent_mito'].mean())
    elif 'mitochondria_ratio' in adata.obs.columns:
        # GSE215968 uses mitochondria_ratio (0-1 range)
        stats['mt_pct_median'] = float(adata.obs['mitochondria_ratio'].median() * 100)
        stats['mt_pct_mean'] = float(adata.obs['mitochondria_ratio'].mean() * 100)

    if 'n_genes_by_counts' in adata.obs.columns:
        stats['n_genes_median'] = int(adata.obs['n_genes_by_counts'].median())
        stats['n_genes_mean'] = int(adata.obs['n_genes_by_counts'].mean())
    elif 'n_genes' in adata.obs.columns:
        stats['n_genes_median'] = int(adata.obs['n_genes'].median())
        stats['n_genes_mean'] = int(adata.obs['n_genes'].mean())
    elif 'nFeature_RNA' in adata.obs.columns:
        stats['n_genes_median'] = int(adata.obs['nFeature_RNA'].median())
        stats['n_genes_mean'] = int(adata.obs['nFeature_RNA'].mean())

    if 'total_counts' in adata.obs.columns:
        stats['n_counts_median'] = int(adata.obs['total_counts'].median())
        stats['n_counts_mean'] = int(adata.obs['total_counts'].mean())
    elif 'n_counts' in adata.obs.columns:
        stats['n_counts_median'] = int(adata.obs['n_counts'].median())
        stats['n_counts_mean'] = int(adata.obs['n_counts'].mean())
    elif 'nCount_RNA' in adata.obs.columns:
        stats['n_counts_median'] = int(adata.obs['nCount_RNA'].median())
        stats['n_counts_mean'] = int(adata.obs['nCount_RNA'].mean())

    return stats


# =====================================================================
# Part 1: Validate QC metrics for pre-processed datasets
# =====================================================================
print("=" * 70)
print("Part 1: Validate pre-processed datasets")
print("=" * 70)

# --- HECA reference atlas ---
print("\n--- HECA_reference (cells) ---")
adata = ad.read_h5ad(os.path.join(ORGANIZED, 'HECA_reference/HECA_cells.h5ad'), backed='r')
stats = qc_stats(adata, 'HECA_cells')
qc_summary['HECA_cells'] = stats
print(f"  {stats['n_cells']} cells | mt% median: {stats.get('mt_pct_median','N/A'):.1f} | gene median: {stats.get('n_genes_median','N/A')}")
adata.file.close()

# --- GSE215968 AS vs WOI ---
print("\n--- GSE215968_AS_vs_WOI_Control ---")
adata = ad.read_h5ad(os.path.join(ORGANIZED, 'GSE215968_Asherman/GSE215968_AS_vs_WOI_Control.h5ad'), backed='r')
stats = qc_stats(adata, 'GSE215968_AS_vs_WOI')
qc_summary['GSE215968_AS_vs_WOI'] = stats
print(f"  {stats['n_cells']} cells | mt% median: {stats.get('mt_pct_median','N/A'):.1f} | gene median: {stats.get('n_genes_median','N/A')}")
adata.file.close()

# --- GSE215968 AS pre vs post ---
print("\n--- GSE215968_AS_pre_vs_post ---")
adata = ad.read_h5ad(os.path.join(ORGANIZED, 'GSE215968_Asherman/GSE215968_AS_pre_vs_post.h5ad'), backed='r')
stats = qc_stats(adata, 'GSE215968_AS_pre_vs_post')
qc_summary['GSE215968_AS_pre_vs_post'] = stats
print(f"  {stats['n_cells']} cells | mt% median: {stats.get('mt_pct_median','N/A'):.1f}")
adata.file.close()

# --- GSE215968 CD133 ---
print("\n--- GSE215968_CD133 ---")
adata = ad.read_h5ad(os.path.join(ORGANIZED, 'GSE215968_Asherman/GSE215968_CD133.h5ad'), backed='r')
stats = qc_stats(adata, 'GSE215968_CD133')
qc_summary['GSE215968_CD133'] = stats
print(f"  {stats['n_cells']} cells | mt% median: {stats.get('mt_pct_median','N/A'):.1f}")
adata.file.close()

# --- GSE216748 organoid ---
print("\n--- GSE216748_organoid ---")
adata = ad.read_h5ad(os.path.join(ORGANIZED, 'GSE216748_organoid/GSE216748_organoid.h5ad'), backed='r')
stats = qc_stats(adata, 'GSE216748_organoid')
qc_summary['GSE216748_organoid'] = stats
print(f"  {stats['n_cells']} cells | mt% median: {stats.get('mt_pct_median','N/A'):.1f}")
adata.file.close()


# =====================================================================
# Part 2: Full QC pipeline for raw datasets
# =====================================================================
print("\n" + "=" * 70)
print("Part 2: Full QC pipeline for raw datasets")
print("=" * 70)

# --- GSE111976 menstrual cycle ---
print("\n--- GSE111976_menstrual (71,032 cells) ---")
adata = ad.read_h5ad(os.path.join(ORGANIZED, 'GSE111976_menstrual/GSE111976_menstrual.h5ad'))

# Record pre-filter statistics
qc_summary['GSE111976_raw'] = {'dataset': 'GSE111976_raw', 'n_cells': adata.n_obs, 'n_genes': adata.n_vars}

adata = compute_qc_metrics(adata, 'GSE111976')
stats_pre = qc_stats(adata, 'GSE111976_pre_filter')
print(f"  Pre-filter: mt% median={stats_pre.get('mt_pct_median',0):.1f}, gene median={stats_pre.get('n_genes_median',0)}, UMI median={stats_pre.get('n_counts_median',0)}")

adata, n_before, n_after = filter_cells(adata, 'GSE111976')
stats_post = qc_stats(adata, 'GSE111976', is_filtered=True)
stats_post['n_cells_before'] = n_before
qc_summary['GSE111976'] = stats_post
print(f"  Post-filter: {n_after} cells | mt% median={stats_post.get('mt_pct_median',0):.1f}, gene median={stats_post.get('n_genes_median',0)}")

# Save filtered data
out_path = os.path.join(OUTPUT, 'GSE111976_qc.h5ad')
adata.write_h5ad(out_path)
print(f"  Saved: {out_path}")
del adata

# --- E-MTAB-10287 temporal dynamics (11 samples merged) ---
print("\n--- E-MTAB-10287_temporal (11 samples, 55,729 cells) ---")
temporal_dir = os.path.join(ORGANIZED, 'E-MTAB-10287_temporal')
h5ad_files = sorted([f for f in os.listdir(temporal_dir) if f.endswith('.h5ad')])

adatas = []
for f in h5ad_files:
    a = ad.read_h5ad(os.path.join(temporal_dir, f))
    adatas.append(a)
adata = ad.concat(adatas, join='outer')
adata.obs_names_make_unique()
print(f"  Merged: {len(h5ad_files)} samples, {adata.n_obs} cells x {adata.n_vars} genes")
del adatas

qc_summary['E-MTAB-10287_raw'] = {'dataset': 'E-MTAB-10287_raw', 'n_cells': adata.n_obs, 'n_genes': adata.n_vars}

adata = compute_qc_metrics(adata, 'E-MTAB-10287')
stats_pre = qc_stats(adata, 'E-MTAB-10287_pre_filter')
print(f"  Pre-filter: mt% median={stats_pre.get('mt_pct_median',0):.1f}, gene median={stats_pre.get('n_genes_median',0)}, UMI median={stats_pre.get('n_counts_median',0)}")

adata, n_before, n_after = filter_cells(adata, 'E-MTAB-10287')
stats_post = qc_stats(adata, 'E-MTAB-10287', is_filtered=True)
stats_post['n_cells_before'] = n_before
qc_summary['E-MTAB-10287'] = stats_post
print(f"  Post-filter: {n_after} cells | mt% median={stats_post.get('mt_pct_median',0):.1f}, gene median={stats_post.get('n_genes_median',0)}")

# Per-sample post-filter distribution
sample_counts = adata.obs['sample'].value_counts()
print(f"  Sample distribution: {dict(sample_counts)}")

out_path = os.path.join(OUTPUT, 'E-MTAB-10287_qc.h5ad')
adata.write_h5ad(out_path)
print(f"  Saved: {out_path}")
del adata

# --- GSE260658 uterus atlas (11 samples) ---
print("\n--- GSE260658_uterus_atlas (11 samples, 71,058 cells) ---")
ut_dir = os.path.join(ORGANIZED, 'GSE260658_uterus_atlas')
h5ad_files = sorted([f for f in os.listdir(ut_dir) if f.endswith('.h5ad')])

adatas = []
for f in h5ad_files:
    a = ad.read_h5ad(os.path.join(ut_dir, f))
    adatas.append(a)
adata = ad.concat(adatas, join='outer')
adata.obs_names_make_unique()
print(f"  Merged: {len(h5ad_files)} samples, {adata.n_obs} cells x {adata.n_vars} genes")
del adatas

qc_summary['GSE260658_raw'] = {'dataset': 'GSE260658_raw', 'n_cells': adata.n_obs, 'n_genes': adata.n_vars}

adata = compute_qc_metrics(adata, 'GSE260658')
stats_pre = qc_stats(adata, 'GSE260658_pre_filter')
print(f"  Pre-filter: mt% median={stats_pre.get('mt_pct_median',0):.1f}, gene median={stats_pre.get('n_genes_median',0)}, UMI median={stats_pre.get('n_counts_median',0)}")

adata, n_before, n_after = filter_cells(adata, 'GSE260658')
stats_post = qc_stats(adata, 'GSE260658', is_filtered=True)
stats_post['n_cells_before'] = n_before
qc_summary['GSE260658'] = stats_post
print(f"  Post-filter: {n_after} cells | mt% median={stats_post.get('mt_pct_median',0):.1f}, gene median={stats_post.get('n_genes_median',0)}")

# Per tissue layer and sample statistics
if 'tissue_layer' in adata.obs.columns:
    layer_counts = adata.obs['tissue_layer'].value_counts()
    print(f"  Tissue layer distribution: {dict(layer_counts)}")

sample_counts = adata.obs['sample'].value_counts()
print(f"  Sample distribution: {dict(sample_counts)}")

out_path = os.path.join(OUTPUT, 'GSE260658_qc.h5ad')
adata.write_h5ad(out_path)
print(f"  Saved: {out_path}")
del adata


# =====================================================================
# Save QC summary
# =====================================================================
print("\n" + "=" * 70)
print("QC Summary")
print("=" * 70)

summary_path = os.path.join(OUTPUT, 'qc_summary.json')
with open(summary_path, 'w') as f:
    json.dump(qc_summary, f, indent=2, ensure_ascii=False)
print(f"\nSummary saved to: {summary_path}")

# Print summary table
print(f"\n{'Dataset':<30} {'Pre-filter':>10} {'Post-filter':>11} {'Kept%':>6} {'mt% med':>7} {'gene med':>8}")
print("-" * 75)
for name, s in qc_summary.items():
    if '_raw' in name:
        continue
    n_before = s.get('n_cells_before', s['n_cells'])
    n_after = s['n_cells']
    pct = 100 * n_after / n_before if n_before > 0 else 100
    mt = s.get('mt_pct_median', -1)
    ng = s.get('n_genes_median', -1)
    mt_str = f"{mt:.1f}" if mt >= 0 else "N/A"
    ng_str = str(ng) if ng >= 0 else "N/A"
    print(f"{name:<30} {n_before:>10} {n_after:>11} {pct:>5.1f}% {mt_str:>7} {ng_str:>8}")

print("\nQC metric calculation and cell filtering complete")
