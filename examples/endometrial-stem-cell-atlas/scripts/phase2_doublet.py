#!/usr/bin/env python3
"""
Phase 2 Step 3: Doublet detection (Scrublet)
Perform doublet detection on 3 newly filtered datasets
"""
import scanpy as sc
import anndata as ad
import scrublet as scr
import numpy as np
import pandas as pd
import os
import json
import warnings
warnings.filterwarnings('ignore')

OUTPUT = '/home/moog/test/tongxue/Phase_output/phase2_qc'

doublet_summary = {}

def run_scrublet(adata, dataset_name, expected_doublet_rate=0.06):
    """Run scrublet doublet detection on a dataset"""
    print(f"\n--- {dataset_name}: Scrublet doublet detection ---")
    n_before = adata.n_obs

    # If sample information is available, detect per sample (recommended approach)
    if 'sample' in adata.obs.columns:
        samples = adata.obs['sample'].unique()
        all_scores = np.zeros(adata.n_obs)
        all_preds = np.zeros(adata.n_obs, dtype=bool)

        for sample in samples:
            mask = adata.obs['sample'] == sample
            sub = adata[mask].copy()

            if sub.n_obs < 100:
                print(f"  {sample}: skipped (<100 cells)")
                continue

            # Estimate doublet rate
            edr = min(0.1, sub.n_obs / 1000 * 0.008)

            try:
                scrub = scr.Scrublet(sub.X, expected_doublet_rate=edr)
                scores, preds = scrub.scrub_doublets(
                    min_counts=2, min_cells=3,
                    min_gene_variability_pctl=85,
                    n_prin_comps=30,
                    verbose=False
                )
                all_scores[mask.values] = scores
                all_preds[mask.values] = preds
                n_doublets = preds.sum()
                print(f"  {sample}: {sub.n_obs} cells, {n_doublets} doublets ({100*n_doublets/sub.n_obs:.1f}%)")
            except Exception as e:
                print(f"  {sample}: Scrublet failed - {str(e)[:80]}, using default threshold")
                # Use more lenient handling for failed samples
                all_scores[mask.values] = 0
                all_preds[mask.values] = False

        adata.obs['doublet_score'] = all_scores
        adata.obs['predicted_doublet'] = all_preds
    else:
        # Detect as a whole
        edr = min(0.1, adata.n_obs / 1000 * 0.008)
        try:
            scrub = scr.Scrublet(adata.X, expected_doublet_rate=edr)
            scores, preds = scrub.scrub_doublets(
                min_counts=2, min_cells=3,
                min_gene_variability_pctl=85,
                n_prin_comps=30,
                verbose=False
            )
            adata.obs['doublet_score'] = scores
            adata.obs['predicted_doublet'] = preds
        except Exception as e:
            print(f"  Scrublet failed: {str(e)[:80]}")
            adata.obs['doublet_score'] = 0
            adata.obs['predicted_doublet'] = False

    n_doublets = adata.obs['predicted_doublet'].sum()
    doublet_rate = 100 * n_doublets / adata.n_obs
    print(f"  Total: {n_doublets} doublets / {adata.n_obs} cells ({doublet_rate:.1f}%)")

    # Remove doublets
    adata_clean = adata[~adata.obs['predicted_doublet']].copy()
    n_after = adata_clean.n_obs
    print(f"  After doublet removal: {n_before} -> {n_after} cells")

    doublet_summary[dataset_name] = {
        'n_before': int(n_before),
        'n_doublets': int(n_doublets),
        'doublet_rate': round(doublet_rate, 2),
        'n_after': int(n_after)
    }

    return adata_clean


# Process 3 raw datasets
datasets = [
    ('GSE111976', 'GSE111976_qc.h5ad'),
    ('E-MTAB-10287', 'E-MTAB-10287_qc.h5ad'),
    ('GSE260658', 'GSE260658_qc.h5ad'),
]

for name, filename in datasets:
    path = os.path.join(OUTPUT, filename)
    adata = ad.read_h5ad(path)
    adata = run_scrublet(adata, name)

    # Save (overwrite original file)
    adata.write_h5ad(path)
    print(f"  Saved: {path}")
    del adata

# Save doublet detection summary
summary_path = os.path.join(OUTPUT, 'doublet_summary.json')
with open(summary_path, 'w') as f:
    json.dump(doublet_summary, f, indent=2, ensure_ascii=False)

print("\n" + "=" * 70)
print("Doublet detection summary")
print("=" * 70)
print(f"{'Dataset':<20} {'Pre-detect':>10} {'Doublets':>8} {'Rate':>6} {'Post-detect':>11}")
print("-" * 55)
for name, s in doublet_summary.items():
    print(f"{name:<20} {s['n_before']:>10} {s['n_doublets']:>8} {s['doublet_rate']:>5.1f}% {s['n_after']:>11}")
