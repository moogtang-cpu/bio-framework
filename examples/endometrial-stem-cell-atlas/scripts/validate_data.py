#!/usr/bin/env python3
"""
Data validation script - check the integrity and readability of all downloaded datasets
"""
import os
import sys
import gzip
import json

DATA_ROOT = "/home/moog/test/tongxue/data"

results = {}

def check_file_readable(path):
    """Check if file exists and is readable"""
    if not os.path.exists(path):
        return False, "File does not exist"
    if os.path.getsize(path) == 0:
        return False, "File is empty"
    return True, f"Size: {os.path.getsize(path) / (1024*1024):.1f} MB"

def check_h5ad(path, name):
    """Check h5ad file"""
    ok, msg = check_file_readable(path)
    if not ok:
        return {"status": "FAIL", "message": msg}
    try:
        import anndata
        if path.endswith('.gz'):
            adata = anndata.read_h5ad(path, backed='r')
        else:
            adata = anndata.read_h5ad(path, backed='r')
        info = {
            "status": "OK",
            "shape": f"{adata.n_obs} cells x {adata.n_vars} genes",
            "obs_columns": list(adata.obs.columns)[:10],
            "file_size": msg
        }
        adata.file.close()
        return info
    except Exception as e:
        return {"status": "WARN", "message": f"h5ad read error: {str(e)[:100]}", "file_size": msg}

def check_rds(path, name):
    """Check RDS file (only validate existence and size)"""
    ok, msg = check_file_readable(path)
    if not ok:
        return {"status": "FAIL", "message": msg}
    return {"status": "OK", "file_size": msg, "note": "RDS file, requires R to read"}

def check_mtx_set(directory, name):
    """Check MTX + cells + features file set"""
    mtx_files = [f for f in os.listdir(directory) if f.endswith('.mtx')]
    cell_files = [f for f in os.listdir(directory) if f.endswith('_cells.tsv') or f == 'barcodes.tsv.gz']
    feat_files = [f for f in os.listdir(directory) if f.endswith('_features.tsv') or f == 'features.tsv.gz']

    if not mtx_files:
        return {"status": "FAIL", "message": "No MTX files found"}

    # Check the triplet for each sample
    samples = set()
    for f in mtx_files:
        sample = f.replace('.mtx', '')
        samples.add(sample)

    sample_info = {}
    for s in sorted(samples):
        has_mtx = os.path.exists(os.path.join(directory, f"{s}.mtx"))
        has_cells = os.path.exists(os.path.join(directory, f"{s}_cells.tsv"))
        has_feats = os.path.exists(os.path.join(directory, f"{s}_features.tsv"))
        mtx_size = os.path.getsize(os.path.join(directory, f"{s}.mtx")) / (1024*1024) if has_mtx else 0
        sample_info[s] = {
            "mtx": has_mtx,
            "cells": has_cells,
            "features": has_feats,
            "mtx_MB": round(mtx_size, 1)
        }

    complete = all(v["mtx"] and v["cells"] and v["features"] for v in sample_info.values())
    return {
        "status": "OK" if complete else "WARN",
        "n_samples": len(samples),
        "samples": sample_info
    }

def check_tar(path, name):
    """Check tar file"""
    ok, msg = check_file_readable(path)
    if not ok:
        return {"status": "FAIL", "message": msg}
    return {"status": "OK", "file_size": msg, "note": "TAR archive, requires extraction"}

def check_gz_text(path, name):
    """Check gzip-compressed text file"""
    ok, msg = check_file_readable(path)
    if not ok:
        return {"status": "FAIL", "message": msg}
    try:
        with gzip.open(path, 'rt') as f:
            header = f.readline().strip()
            n_lines = 1
            for _ in f:
                n_lines += 1
                if n_lines > 5:
                    break
        return {"status": "OK", "file_size": msg, "header_preview": header[:100], "lines_sampled": n_lines}
    except Exception as e:
        return {"status": "WARN", "message": f"Read error: {str(e)[:100]}", "file_size": msg}

print("=" * 70)
print("Data validation report - Human endometrial basal layer stem cell atlas")
print("=" * 70)

# === 1. HECA reference atlas ===
print("\n--- 1. HECA_reference (integrated reference atlas) ---")
heca_dir = os.path.join(DATA_ROOT, "scRNA_seq/HECA_reference")
for f in ["endometriumAtlasV2_cells_with_counts.h5ad", "endometriumAtlasV2_nuclei.h5ad", "model.pt"]:
    path = os.path.join(heca_dir, f)
    if f.endswith('.h5ad'):
        result = check_h5ad(path, f)
    else:
        ok, msg = check_file_readable(path)
        result = {"status": "OK" if ok else "FAIL", "file_size": msg}
    print(f"  {f}: {result['status']} - {result.get('shape', result.get('file_size', ''))}")

# === 2. GSE215968 Asherman ===
print("\n--- 2. GSE215968_Asherman (AS patients) ---")
ash_dir = os.path.join(DATA_ROOT, "scRNA_seq/GSE215968_Asherman")
for f in os.listdir(ash_dir):
    path = os.path.join(ash_dir, f)
    if f.endswith('.h5ad.gz'):
        result = check_h5ad(path, f)
        print(f"  {f}: {result['status']} - {result.get('shape', result.get('message', result.get('file_size', '')))}")

# === 3. GSE216748 organoids ===
print("\n--- 3. GSE216748_organoid (AS organoids) ---")
org_dir = os.path.join(DATA_ROOT, "scRNA_seq/GSE216748_organoid")
for f in os.listdir(org_dir):
    path = os.path.join(org_dir, f)
    if f.endswith('.h5ad.gz'):
        result = check_h5ad(path, f)
        print(f"  {f}: {result['status']} - {result.get('shape', result.get('message', result.get('file_size', '')))}")

# === 4. E-MTAB-10287 temporal dynamics ===
print("\n--- 4. E-MTAB-10287_temporal (scRNA-seq) ---")
temporal_dir = os.path.join(DATA_ROOT, "scRNA_seq/E-MTAB-10287_temporal")
result = check_mtx_set(temporal_dir, "E-MTAB-10287")
print(f"  Status: {result['status']}, number of samples: {result.get('n_samples', 'N/A')}")
if 'samples' in result:
    for s, info in sorted(result['samples'].items()):
        complete = "complete" if all([info['mtx'], info['cells'], info['features']]) else "incomplete"
        print(f"    {s}: {complete} (MTX: {info['mtx_MB']}MB)")

# === 5. GSE111976 menstrual cycle ===
print("\n--- 5. GSE111976_menstrual (menstrual cycle atlas) ---")
men_dir = os.path.join(DATA_ROOT, "scRNA_seq/GSE111976_menstrual")
for f in sorted(os.listdir(men_dir)):
    path = os.path.join(men_dir, f)
    if f.endswith('.rds.gz'):
        result = check_rds(path, f)
    elif f.endswith('.csv.gz'):
        result = check_gz_text(path, f)
    else:
        ok, msg = check_file_readable(path)
        result = {"status": "OK" if ok else "FAIL", "file_size": msg}
    print(f"  {f}: {result['status']} - {result.get('file_size', '')}")

# === 6. GSE260658 uterus atlas ===
print("\n--- 6. GSE260658_uterus_atlas (uterine cell atlas) ---")
ut_dir = os.path.join(DATA_ROOT, "scRNA_seq/GSE260658_uterus_atlas")
for f in os.listdir(ut_dir):
    path = os.path.join(ut_dir, f)
    result = check_tar(path, f) if f.endswith('.tar') else check_file_readable(path)
    if isinstance(result, tuple):
        ok, msg = result
        result = {"status": "OK" if ok else "FAIL", "file_size": msg}
    print(f"  {f}: {result['status']} - {result.get('file_size', '')}")

# === 7. E-MTAB-9260 Visium ===
print("\n--- 7. E-MTAB-9260_visium (spatial transcriptomics) ---")
vis_dir = os.path.join(DATA_ROOT, "spatial/E-MTAB-9260_visium")
for f in sorted(os.listdir(vis_dir)):
    path = os.path.join(vis_dir, f)
    if f.endswith('.tsv.gz'):
        result = check_gz_text(path, f)
    else:
        ok, msg = check_file_readable(path)
        result = {"status": "OK" if ok else "FAIL", "file_size": msg}
    print(f"  {f}: {result['status']} - {result.get('file_size', '')}")

# === 8. GSE287278 RIF Visium ===
print("\n--- 8. GSE287278_RIF_visium (RIF spatial transcriptomics) ---")
rif_dir = os.path.join(DATA_ROOT, "spatial/GSE287278_RIF_visium")
for f in os.listdir(rif_dir):
    path = os.path.join(rif_dir, f)
    result = check_tar(path, f) if f.endswith('.tar') else check_file_readable(path)
    if isinstance(result, tuple):
        ok, msg = result
        result = {"status": "OK" if ok else "FAIL", "file_size": msg}
    print(f"  {f}: {result['status']} - {result.get('file_size', '')}")

# === 9. GSE234354 molecular staging ===
print("\n--- 9. GSE234354_staging (Bulk RNA-seq molecular staging) ---")
stag_dir = os.path.join(DATA_ROOT, "bulk_rnaseq/GSE234354_staging")
for f in os.listdir(stag_dir):
    path = os.path.join(stag_dir, f)
    result = check_gz_text(path, f) if f.endswith('.gz') else check_file_readable(path)
    if isinstance(result, tuple):
        ok, msg = result
        result = {"status": "OK" if ok else "FAIL", "file_size": msg}
    print(f"  {f}: {result['status']} - {result.get('file_size', '')}")

# === 10. GSE127918 decidualization ===
print("\n--- 10. GSE127918_decidual (decidualization pathway reference) ---")
dec_dir = os.path.join(DATA_ROOT, "bulk_rnaseq/GSE127918_decidual")
for f in sorted(os.listdir(dec_dir)):
    path = os.path.join(dec_dir, f)
    result = check_gz_text(path, f) if f.endswith('.gz') else check_file_readable(path)
    if isinstance(result, tuple):
        ok, msg = result
        result = {"status": "OK" if ok else "FAIL", "file_size": msg}
    print(f"  {f}: {result['status']} - {result.get('file_size', '')}")

# === 11. PRJNA730360 thin endometrium ===
print("\n--- 11. PRJNA730360_thin (thin endometrium, pending download) ---")
thin_dir = os.path.join(DATA_ROOT, "scRNA_seq/PRJNA730360_thin")
files = [f for f in os.listdir(thin_dir) if not f.startswith('.') and f != 'download_and_process.sh']
if not files:
    print("  Status: PENDING - raw SRA data downloading")
else:
    for f in files:
        print(f"  {f}")

print("\n" + "=" * 70)
print("Validation complete")
print("=" * 70)
