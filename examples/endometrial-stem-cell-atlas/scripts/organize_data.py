#!/usr/bin/env python3
"""
Data organization script - unify data formats, create standardized directory structure
Prepares data for subsequent Phase 2 (QC and preprocessing)
"""
import os
import sys
import json

# Run in conda environment
import anndata as ad
import scanpy as sc
import pandas as pd
import scipy.io
import scipy.sparse
import numpy as np

DATA_ROOT = "/home/moog/test/tongxue/data"
OUTPUT_DIR = os.path.join(DATA_ROOT, "organized")
os.makedirs(OUTPUT_DIR, exist_ok=True)

def save_dataset_info(name, adata, output_path, source, notes=""):
    """Save dataset metadata"""
    info = {
        "name": name,
        "source": source,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "obs_columns": list(adata.obs.columns),
        "output_path": output_path,
        "notes": notes
    }
    info_path = output_path.replace('.h5ad', '_info.json')
    with open(info_path, 'w') as f:
        json.dump(info, f, indent=2, ensure_ascii=False)
    return info

print("=" * 70)
print("Data organization - unified formatting to h5ad")
print("=" * 70)

# === 1. HECA reference atlas (already in h5ad format, create symbolic links) ===
print("\n[1/10] HECA_reference...")
heca_cells = os.path.join(DATA_ROOT, "scRNA_seq/HECA_reference/endometriumAtlasV2_cells_with_counts.h5ad")
heca_nuclei = os.path.join(DATA_ROOT, "scRNA_seq/HECA_reference/endometriumAtlasV2_nuclei.h5ad")
heca_out = os.path.join(OUTPUT_DIR, "HECA_reference")
os.makedirs(heca_out, exist_ok=True)
for src, name in [(heca_cells, "cells"), (heca_nuclei, "nuclei")]:
    dst = os.path.join(heca_out, f"HECA_{name}.h5ad")
    if not os.path.exists(dst):
        os.symlink(src, dst)
    adata = ad.read_h5ad(src, backed='r')
    save_dataset_info(f"HECA_{name}", adata, dst, "Nature Genetics 2024")
    print(f"  HECA_{name}: {adata.n_obs} cells x {adata.n_vars} genes [linked]")
    adata.file.close()

# === 2. GSE215968 Asherman syndrome ===
print("\n[2/10] GSE215968_Asherman...")
ash_out = os.path.join(OUTPUT_DIR, "GSE215968_Asherman")
os.makedirs(ash_out, exist_ok=True)
ash_files = {
    "AS_vs_WOI_Control": "GSE215968_sc_AS_vs_WOI_Control.h5ad",
    "AS_pre_vs_post": "GSE215968_sc_AS_pre_vs_post.h5ad",
    "CD133": "GSE215968_sc_CD133.h5ad"
}
for name, fname in ash_files.items():
    src = os.path.join(DATA_ROOT, f"scRNA_seq/GSE215968_Asherman/{fname}")
    dst = os.path.join(ash_out, f"GSE215968_{name}.h5ad")
    if not os.path.exists(dst):
        os.symlink(src, dst)
    adata = ad.read_h5ad(src, backed='r')
    save_dataset_info(f"GSE215968_{name}", adata, dst, "Nature Communications 2023")
    print(f"  {name}: {adata.n_obs} cells x {adata.n_vars} genes [linked]")
    adata.file.close()

# === 3. GSE216748 organoids ===
print("\n[3/10] GSE216748_organoid...")
org_out = os.path.join(OUTPUT_DIR, "GSE216748_organoid")
os.makedirs(org_out, exist_ok=True)
src = os.path.join(DATA_ROOT, "scRNA_seq/GSE216748_organoid/GSE216748_sc_AS_org_hormones.h5ad")
dst = os.path.join(org_out, "GSE216748_organoid.h5ad")
if not os.path.exists(dst):
    os.symlink(src, dst)
adata = ad.read_h5ad(src, backed='r')
save_dataset_info("GSE216748_organoid", adata, dst, "Nature Communications 2023")
print(f"  organoid: {adata.n_obs} cells x {adata.n_vars} genes [linked]")
adata.file.close()

# === 4. E-MTAB-10287 temporal dynamics (MTX -> h5ad) ===
print("\n[4/10] E-MTAB-10287_temporal...")
temporal_out = os.path.join(OUTPUT_DIR, "E-MTAB-10287_temporal")
os.makedirs(temporal_out, exist_ok=True)
temporal_dir = os.path.join(DATA_ROOT, "scRNA_seq/E-MTAB-10287_temporal")

# Read SDRF metadata
sdrf_path = os.path.join(temporal_dir, "E-MTAB-10287.sdrf.txt")
if os.path.exists(sdrf_path):
    sdrf = pd.read_csv(sdrf_path, sep='\t')

# Find all samples
mtx_files = [f.replace('.mtx', '') for f in os.listdir(temporal_dir) if f.endswith('.mtx')]
total_cells = 0
for sample in sorted(mtx_files):
    h5ad_path = os.path.join(temporal_out, f"{sample}.h5ad")
    if os.path.exists(h5ad_path):
        adata = ad.read_h5ad(h5ad_path, backed='r')
        total_cells += adata.n_obs
        adata.file.close()
        print(f"  {sample}: already exists [skipped]")
        continue

    mtx = scipy.io.mmread(os.path.join(temporal_dir, f"{sample}.mtx"))
    cells = pd.read_csv(os.path.join(temporal_dir, f"{sample}_cells.tsv"), header=0, sep='\t')
    features = pd.read_csv(os.path.join(temporal_dir, f"{sample}_features.tsv"), header=None, sep='\t')

    # MTX is in gene x cell format, needs to be transposed
    adata = ad.AnnData(X=scipy.sparse.csr_matrix(mtx.T))
    # cells_tsv first row is header, use header=0 to skip
    adata.obs_names = cells.iloc[:, 0].values.astype(str)
    adata.var_names = features[0].values
    if features.shape[1] > 1:
        adata.var['gene_name'] = features[1].values

    adata.obs['sample'] = sample
    adata.obs['dataset'] = 'E-MTAB-10287'

    adata.write_h5ad(h5ad_path)
    total_cells += adata.n_obs
    print(f"  {sample}: {adata.n_obs} cells x {adata.n_vars} genes [converted MTX->h5ad]")

print(f"  E-MTAB-10287 total: {total_cells} cells")

# === 5. GSE111976 menstrual cycle (RDS -> requires R conversion, log info first) ===
print("\n[5/10] GSE111976_menstrual...")
men_out = os.path.join(OUTPUT_DIR, "GSE111976_menstrual")
os.makedirs(men_out, exist_ok=True)
# RDS requires R to read, create a marker file first
rds_path = os.path.join(DATA_ROOT, "scRNA_seq/GSE111976_menstrual/GSE111976_ct_endo_10x.rds.gz")
csv_path = os.path.join(DATA_ROOT, "scRNA_seq/GSE111976_menstrual/GSE111976_ct.csv.gz")
umap_path = os.path.join(DATA_ROOT, "scRNA_seq/GSE111976_menstrual/GSE111976_umap_endo_10x.csv.gz")

# Read CSV file first to check format
ct = pd.read_csv(csv_path, index_col=0, nrows=5)
print(f"  ct.csv: {ct.shape[1]} columns, cols: {list(ct.columns)[:5]}...")

umap = pd.read_csv(umap_path, index_col=0, nrows=5)
print(f"  umap.csv: {umap.shape[1]} columns, cols: {list(umap.columns)}")

# RDS file needs to be converted using R
marker = os.path.join(men_out, "NEEDS_R_CONVERSION.txt")
with open(marker, 'w') as f:
    f.write(f"RDS source file: {rds_path}\n")
    f.write("Need to run in R: Seurat::SaveH5Seurat() or SeuratDisk conversion\n")
    f.write("CSV count matrix and UMAP coordinates are ready\n")
print("  RDS file requires R conversion [marked]")

# === 6. GSE260658 uterus atlas (10x MTX -> h5ad) ===
print("\n[6/10] GSE260658_uterus_atlas...")
ut_out = os.path.join(OUTPUT_DIR, "GSE260658_uterus_atlas")
os.makedirs(ut_out, exist_ok=True)
ut_dir = os.path.join(DATA_ROOT, "scRNA_seq/GSE260658_uterus_atlas")

# Group by sample
samples = {}
for f in os.listdir(ut_dir):
    if f.endswith('_matrix.mtx.gz'):
        parts = f.split('_')
        gsm = parts[0]
        sample_name = '_'.join(parts[1:-1])
        samples[sample_name] = gsm

total_cells = 0
for sample_name, gsm in sorted(samples.items()):
    h5ad_path = os.path.join(ut_out, f"{sample_name}.h5ad")
    if os.path.exists(h5ad_path):
        adata = ad.read_h5ad(h5ad_path, backed='r')
        total_cells += adata.n_obs
        adata.file.close()
        print(f"  {sample_name}: already exists [skipped]")
        continue

    adata = sc.read_10x_mtx(
        ut_dir,
        prefix=f"{gsm}_{sample_name}_",
        var_names='gene_symbols',
        cache=False
    )
    adata.obs['sample'] = sample_name
    adata.obs['gsm'] = gsm
    adata.obs['dataset'] = 'GSE260658'

    # Label tissue type
    if 'myo' in sample_name:
        adata.obs['tissue_layer'] = 'myometrium'
    elif 'endo' in sample_name:
        adata.obs['tissue_layer'] = 'endometrium'
    elif sample_name.startswith('d'):
        adata.obs['tissue_layer'] = 'decidua'
    else:
        adata.obs['tissue_layer'] = 'unknown'

    adata.write_h5ad(h5ad_path)
    total_cells += adata.n_obs
    print(f"  {sample_name}: {adata.n_obs} cells x {adata.n_vars} genes [converted 10x->h5ad]")

print(f"  GSE260658 total: {total_cells} cells")

# === 7. E-MTAB-9260 Visium (count matrix available, defer conversion until spatial coords are ready) ===
print("\n[7/10] E-MTAB-9260_visium...")
vis_out = os.path.join(OUTPUT_DIR, "E-MTAB-9260_visium")
os.makedirs(vis_out, exist_ok=True)
vis_dir = os.path.join(DATA_ROOT, "spatial/E-MTAB-9260_visium")
for f in os.listdir(vis_dir):
    if f.startswith('Visium_raw_counts'):
        src = os.path.join(vis_dir, f)
        dst = os.path.join(vis_out, f)
        if not os.path.exists(dst):
            os.symlink(src, dst)
print("  Count matrix linked, waiting for spatial coordinates download to complete before integration")

# === 8. GSE287278 RIF Visium ===
print("\n[8/10] GSE287278_RIF_visium...")
rif_out = os.path.join(OUTPUT_DIR, "GSE287278_RIF_visium")
os.makedirs(rif_out, exist_ok=True)
rif_dir = os.path.join(DATA_ROOT, "spatial/GSE287278_RIF_visium")
total_spots = 0
for sample_dir in sorted(os.listdir(rif_dir)):
    sr_path = os.path.join(rif_dir, sample_dir, f"{sample_dir}_processed_data")
    mtx_dir = os.path.join(sr_path, "filtered_feature_bc_matrix")
    spatial_dir = os.path.join(sr_path, "spatial")
    if os.path.isdir(sr_path) and os.path.isdir(mtx_dir):
        h5ad_path = os.path.join(rif_out, f"{sample_dir}.h5ad")
        if os.path.exists(h5ad_path):
            adata = ad.read_h5ad(h5ad_path, backed='r')
            total_spots += adata.n_obs
            adata.file.close()
            print(f"  {sample_dir}: already exists [skipped]")
            continue

        # Manually read 10x MTX + spatial information
        adata = sc.read_10x_mtx(mtx_dir, var_names='gene_symbols', cache=False)
        adata.obs['sample'] = sample_dir
        adata.obs['condition'] = 'CTR' if 'CTR' in sample_dir else 'RIF'
        adata.obs['dataset'] = 'GSE287278'

        # Read spatial coordinates
        positions_path = os.path.join(spatial_dir, "tissue_positions.csv")
        if os.path.exists(positions_path):
            positions = pd.read_csv(positions_path, index_col=0)
            # Ensure barcode matching
            common = adata.obs_names.intersection(positions.index)
            if len(common) > 0:
                positions = positions.loc[common]
                adata = adata[common].copy()
                adata.obsm['spatial'] = positions.iloc[:, [3, 4]].values  # pixel coordinates
                adata.obs['in_tissue'] = positions.iloc[:, 0].values

        # Read scale factors
        sf_path = os.path.join(spatial_dir, "scalefactors_json.json")
        if os.path.exists(sf_path):
            with open(sf_path) as f:
                adata.uns['spatial'] = {sample_dir: {'scalefactors': json.load(f)}}

        adata.write_h5ad(h5ad_path)
        total_spots += adata.n_obs
        print(f"  {sample_dir}: {adata.n_obs} spots x {adata.n_vars} genes [Visium MTX->h5ad]")

print(f"  GSE287278 total: {total_spots} spots")

# === 9. GSE234354 Bulk RNA-seq ===
print("\n[9/10] GSE234354_staging (Bulk)...")
bulk_out = os.path.join(OUTPUT_DIR, "bulk")
os.makedirs(bulk_out, exist_ok=True)
bulk_path = os.path.join(DATA_ROOT, "bulk_rnaseq/GSE234354_staging/GSE234354_gene_count_matrix.txt.gz")
dst = os.path.join(bulk_out, "GSE234354_gene_count_matrix.txt.gz")
if not os.path.exists(dst):
    os.symlink(bulk_path, dst)
df = pd.read_csv(bulk_path, sep='\t', nrows=5)
print(f"  GSE234354: {df.shape[1]-1} samples, cols: {list(df.columns)[:5]}...")

# === 10. GSE127918 decidualization ===
print("\n[10/10] GSE127918_decidual (Bulk)...")
for f in ["GSE127918_Biopsy_combined_dge.txt.gz", "GSE127918_Timecourse_combined_dge.txt.gz"]:
    src = os.path.join(DATA_ROOT, f"bulk_rnaseq/GSE127918_decidual/{f}")
    dst = os.path.join(bulk_out, f)
    if not os.path.exists(dst):
        os.symlink(src, dst)
    df = pd.read_csv(src, sep='\t', nrows=5)
    print(f"  {f}: {df.shape[1]} columns")

print("\n" + "=" * 70)
print("Data organization complete")
print("=" * 70)

# Output organized directory structure
print("\nOrganized directory structure:")
for d in sorted(os.listdir(OUTPUT_DIR)):
    full_path = os.path.join(OUTPUT_DIR, d)
    if os.path.isdir(full_path):
        files = os.listdir(full_path)
        h5ad_files = [f for f in files if f.endswith('.h5ad')]
        print(f"  {d}/: {len(h5ad_files)} h5ad + {len(files)-len(h5ad_files)} other files")
