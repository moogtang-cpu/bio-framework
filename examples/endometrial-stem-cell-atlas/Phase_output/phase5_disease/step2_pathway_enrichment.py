#!/usr/bin/env python
"""
Phase 5 Step 2: Pathway enrichment analysis (GSEA/GO/KEGG/Reactome)
Focused on Wnt/Notch/TGF-beta signaling pathways + Fibrosis vs Regeneration
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
import json
import os
import warnings
warnings.filterwarnings('ignore')

FIGDIR = 'Phase_output/phase5_disease/figures/'
OUTDIR = 'Phase_output/phase5_disease/'
os.makedirs(FIGDIR, exist_ok=True)

# ============================================================
# 1. Load DEG results
# ============================================================
print("=" * 60)
print("Phase 5 Step 2: Pathway Enrichment Analysis")
print("=" * 60)

# Load all DEG results
deg_files = [f for f in os.listdir(OUTDIR) if f.startswith('deg_') and f.endswith('.csv')]
all_degs = {}
for f in sorted(deg_files):
    ct_name = f.replace('deg_', '').replace('_AS_vs_WOI.csv', '')
    df = pd.read_csv(os.path.join(OUTDIR, f), index_col=0)
    all_degs[ct_name] = df
    sig = df[(df['padj'] < 0.05) & (abs(df['log2FoldChange']) > 0.25)]
    print(f"  {ct_name}: {len(df)} genes, {len(sig)} DEGs")

# ============================================================
# 2. GSEA preranked analysis (per cell type)
# ============================================================
print("\n--- GSEA Preranked Analysis ---")

# Use MSigDB gene sets
gene_sets = [
    'GO_Biological_Process_2023',
    'KEGG_2021_Human',
    'Reactome_2022',
    'WikiPathway_2023_Human',
]

# Custom stem cell / pathway gene sets
custom_genesets = {
    'WNT_SIGNALING_STEM': ['WNT2', 'WNT3A', 'WNT4', 'WNT5A', 'WNT7A', 'WNT7B', 'WNT10A',
                           'RSPO1', 'RSPO3', 'CTNNB1', 'LEF1', 'TCF7L2', 'AXIN2', 'LGR5',
                           'SFRP1', 'SFRP4', 'DKK1', 'DKK3', 'WIF1', 'FRZB', 'GSK3B'],
    'NOTCH_SIGNALING_STEM': ['NOTCH1', 'NOTCH2', 'NOTCH3', 'JAG1', 'JAG2', 'DLL1', 'DLL4',
                              'HES1', 'HEY1', 'HEY2', 'RBPJ', 'MAML1', 'NUMB'],
    'TGFB_SIGNALING': ['TGFB1', 'TGFB2', 'TGFB3', 'TGFBR1', 'TGFBR2', 'SMAD2', 'SMAD3',
                        'SMAD4', 'SMAD7', 'ACVR1', 'BMP2', 'BMP4', 'BMP7', 'INHBA', 'INHBB'],
    'FIBROSIS_MARKERS': ['COL1A1', 'COL1A2', 'COL3A1', 'COL5A1', 'COL6A1', 'FN1',
                          'ACTA2', 'FAP', 'POSTN', 'CTGF', 'SERPINE1', 'LOX', 'LOXL2',
                          'MMP2', 'MMP9', 'TIMP1', 'TIMP2', 'TIMP3'],
    'ENDOMETRIAL_REGENERATION': ['SOX9', 'LGR5', 'PROM1', 'ALDH1A1', 'TERT', 'BMI1',
                                  'CD44', 'SUSD2', 'THY1', 'SSEA1', 'ITGA6',
                                  'MMP7', 'MMP11', 'SPP1', 'PRL', 'IGFBP1'],
    'STEM_NICHE_FACTORS': ['RSPO1', 'RSPO3', 'WNT2', 'WNT4', 'BMP2', 'BMP4',
                            'PDGFRA', 'PDGFRB', 'CXCL12', 'SCF', 'EGF', 'FGF2', 'FGF10',
                            'IGF1', 'HGF', 'VEGFA'],
}

gsea_results = {}
gsea_top_pathways = {}

for ct_name, deg_df in all_degs.items():
    print(f"\n  Processing {ct_name}...")

    # Preranking: use -log10(pvalue) * sign(log2FC)
    rank_df = deg_df[['log2FoldChange', 'pvalue']].dropna().copy()
    rank_df['rank_score'] = -np.log10(rank_df['pvalue'].clip(lower=1e-300)) * np.sign(rank_df['log2FoldChange'])
    rank_df = rank_df.sort_values('rank_score', ascending=False)
    rnk = rank_df['rank_score']

    if len(rnk) < 100:
        print(f"    Skipping: only {len(rnk)} genes")
        continue

    # 2a. Standard GSEA (GO/KEGG/Reactome)
    ct_results = {}
    for gs in gene_sets:
        try:
            res = gp.prerank(
                rnk=rnk,
                gene_sets=gs,
                threads=8,
                min_size=15,
                max_size=500,
                permutation_num=1000,
                seed=42,
                no_plot=True,
            )
            sig_terms = res.res2d[res.res2d['FDR q-val'] < 0.05]
            ct_results[gs] = res.res2d.copy()
            print(f"    {gs}: {len(sig_terms)} enriched (FDR<0.05)")
        except Exception as e:
            print(f"    {gs}: error - {e}")

    # 2b. Custom gene set GSEA
    try:
        res_custom = gp.prerank(
            rnk=rnk,
            gene_sets=custom_genesets,
            threads=8,
            min_size=5,
            max_size=500,
            permutation_num=1000,
            seed=42,
            no_plot=True,
        )
        ct_results['Custom_Pathways'] = res_custom.res2d.copy()
        print(f"    Custom: {len(res_custom.res2d[res_custom.res2d['FDR q-val'] < 0.25])} enriched (FDR<0.25)")
    except Exception as e:
        print(f"    Custom: error - {e}")

    gsea_results[ct_name] = ct_results

    # Extract top pathways
    top_terms = []
    for gs, res_df in ct_results.items():
        top_up = res_df[res_df['NES'] > 0].head(5)
        top_down = res_df[res_df['NES'] < 0].head(5)
        for _, row in pd.concat([top_up, top_down]).iterrows():
            top_terms.append({
                'gene_set': gs,
                'term': row['Term'],
                'nes': float(row['NES']),
                'fdr': float(row['FDR q-val']),
            })
    gsea_top_pathways[ct_name] = top_terms

# ============================================================
# 3. Visualization: stem cell pathway heatmap
# ============================================================
print("\n--- Visualization ---")

# 3a. Custom pathway NES heatmap (all cell types)
custom_nes = {}
for ct_name in gsea_results:
    if 'Custom_Pathways' in gsea_results[ct_name]:
        res = gsea_results[ct_name]['Custom_Pathways']
        for _, row in res.iterrows():
            term = row['Term']
            if term not in custom_nes:
                custom_nes[term] = {}
            custom_nes[term][ct_name] = float(row['NES'])

if custom_nes:
    nes_df = pd.DataFrame(custom_nes).T
    nes_df = nes_df.fillna(0)

    fig, ax = plt.subplots(figsize=(14, 5))
    sns.heatmap(nes_df, cmap='RdBu_r', center=0, annot=True, fmt='.2f',
                linewidths=0.5, ax=ax, cbar_kws={'label': 'NES'})
    ax.set_title('Custom Pathway NES: AS vs WOI Control (GSEA prerank)')
    ax.set_ylabel('')
    plt.tight_layout()
    plt.savefig(FIGDIR + 'gsea_custom_pathways_heatmap.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("  gsea_custom_pathways_heatmap.png")

# 3b. Top GO/KEGG pathways for major cell types
focus_types = ['Decidualized_stromal', 'SOX9plusLGR5plus_stem', 'CD133plus_progenitor',
               'Basalis_fibroblast', 'Glandular_epithelium', 'NKT_cells', 'Macrophages']

for ct_name in focus_types:
    if ct_name not in gsea_results:
        continue

    # Merge results from all standard gene sets
    all_res = []
    for gs in gene_sets:
        if gs in gsea_results[ct_name]:
            res_df = gsea_results[ct_name][gs].copy()
            res_df['source'] = gs
            all_res.append(res_df)

    if not all_res:
        continue

    combined = pd.concat(all_res)
    # Ensure FDR column is numeric
    combined['FDR q-val'] = pd.to_numeric(combined['FDR q-val'], errors='coerce')
    combined['NES'] = pd.to_numeric(combined['NES'], errors='coerce')
    sig_combined = combined[combined['FDR q-val'] < 0.05].copy()

    if len(sig_combined) == 0:
        print(f"  {ct_name}: no significant pathways, skipping")
        continue

    # Top 10 up + top 10 down
    top_up = sig_combined[sig_combined['NES'] > 0].nsmallest(10, 'FDR q-val')
    top_down = sig_combined[sig_combined['NES'] < 0].nsmallest(10, 'FDR q-val')
    plot_df = pd.concat([top_up, top_down])

    if len(plot_df) < 3:
        continue

    fig, ax = plt.subplots(figsize=(12, max(4, len(plot_df) * 0.4)))
    colors = ['#E64B35' if nes > 0 else '#4DBBD5' for nes in plot_df['NES']]
    y_pos = range(len(plot_df))
    ax.barh(y_pos, plot_df['NES'], color=colors, height=0.7)
    ax.set_yticks(y_pos)

    # Truncate long term names
    labels = [t[:60] + '...' if len(t) > 60 else t for t in plot_df['Term']]
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0, color='black', linewidth=0.5)
    ax.set_xlabel('Normalized Enrichment Score (NES)')
    ax.set_title(f'{ct_name}: Top Enriched Pathways (AS vs WOI, FDR<0.05)')

    # Annotate FDR
    for i, (_, row) in enumerate(plot_df.iterrows()):
        ax.text(row['NES'] + 0.05 * np.sign(row['NES']), i,
                f"FDR={row['FDR q-val']:.2e}", va='center', fontsize=6)

    plt.tight_layout()
    plt.savefig(FIGDIR + f'gsea_top_{ct_name}.png', dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  gsea_top_{ct_name}.png")

# 3c. Wnt/Notch/TGF-beta pathway focused plot
# Extract Wnt/Notch/TGF-related pathways across all cell types
pathway_keywords = {
    'Wnt': ['WNT', 'wnt', 'Wnt', 'beta-catenin', 'frizzled'],
    'Notch': ['NOTCH', 'notch', 'Notch'],
    'TGF-beta': ['TGF', 'tgf', 'SMAD', 'smad', 'BMP', 'bmp', 'activin'],
    'Fibrosis': ['fibro', 'Fibro', 'collagen', 'Collagen', 'ECM', 'extracellular matrix'],
}

pathway_nes = {}
for ct_name in gsea_results:
    for gs, res_df in gsea_results[ct_name].items():
        for _, row in res_df.iterrows():
            term = row['Term']
            for pw_name, keywords in pathway_keywords.items():
                if any(kw in term for kw in keywords):
                    key = f"{pw_name}|{term}"
                    if key not in pathway_nes:
                        pathway_nes[key] = {}
                    pathway_nes[key][ct_name] = float(row['NES'])

if pathway_nes:
    pw_df = pd.DataFrame(pathway_nes).T
    pw_df = pw_df.fillna(0)

    # Filter to keep only pathways non-zero in at least 2 cell types
    pw_df = pw_df[(pw_df != 0).sum(axis=1) >= 2]

    if len(pw_df) > 0:
        # Truncate row names
        pw_df.index = [idx.split('|')[0] + ': ' + (idx.split('|')[1][:50] + '...' if len(idx.split('|')[1]) > 50 else idx.split('|')[1])
                       for idx in pw_df.index]

        fig_height = max(6, len(pw_df) * 0.35)
        fig, ax = plt.subplots(figsize=(14, min(fig_height, 20)))
        sns.heatmap(pw_df, cmap='RdBu_r', center=0, linewidths=0.3, ax=ax,
                    cbar_kws={'label': 'NES'}, yticklabels=True)
        ax.set_title('Wnt/Notch/TGF-β/Fibrosis Pathway NES across Cell Types (AS vs WOI)')
        ax.tick_params(axis='y', labelsize=7)
        plt.tight_layout()
        plt.savefig(FIGDIR + 'gsea_stem_pathways_heatmap.png', dpi=200, bbox_inches='tight')
        plt.close()
        print("  gsea_stem_pathways_heatmap.png")

# ============================================================
# 4. ORA analysis (over-representation of significant DEGs)
# ============================================================
print("\n--- ORA Analysis (GO Enrichment) ---")

ora_results = {}
for ct_name, deg_df in all_degs.items():
    sig_genes = deg_df[(deg_df['padj'] < 0.05) & (abs(deg_df['log2FoldChange']) > 0.25)]
    up_genes = sig_genes[sig_genes['log2FoldChange'] > 0].index.tolist()
    down_genes = sig_genes[sig_genes['log2FoldChange'] < 0].index.tolist()

    background = deg_df.index.tolist()

    ct_ora = {}
    for direction, genes, name in [(up_genes, up_genes, 'up'), (down_genes, down_genes, 'down')]:
        if len(genes) < 10:
            continue

        try:
            enr = gp.enrich(
                gene_list=genes,
                gene_sets='GO_Biological_Process_2023',
                background=background,
                no_plot=True,
            )
            sig_terms = enr.res2d[enr.res2d['Adjusted P-value'] < 0.05]
            ct_ora[name] = enr.res2d
            print(f"  {ct_name} ({name}): {len(sig_terms)} GO terms enriched")
        except Exception as e:
            print(f"  {ct_name} ({name}): error - {e}")

    ora_results[ct_name] = ct_ora

# ============================================================
# 5. Save results
# ============================================================
print("\n--- Saving Results ---")

# Save GSEA results
for ct_name, ct_results in gsea_results.items():
    for gs, res_df in ct_results.items():
        gs_short = gs.replace('_2023', '').replace('_2022', '').replace('_2021_Human', '').replace('_Human', '')
        fname = f"gsea_{ct_name}_{gs_short}.csv"
        res_df.to_csv(os.path.join(OUTDIR, fname), index=False)

# Save ORA results
for ct_name, ct_ora in ora_results.items():
    for direction, res_df in ct_ora.items():
        fname = f"ora_{ct_name}_{direction}_GO_BP.csv"
        res_df.to_csv(os.path.join(OUTDIR, fname), index=False)

# Save pathway summary
pathway_report = {
    'phase': 5,
    'step': 2,
    'analysis': 'GSEA + ORA pathway enrichment',
    'gene_sets': gene_sets + ['Custom_Pathways'],
    'custom_pathways': list(custom_genesets.keys()),
    'gsea_summary': {},
    'top_pathways': gsea_top_pathways,
}

for ct_name in gsea_results:
    pathway_report['gsea_summary'][ct_name] = {}
    for gs, res_df in gsea_results[ct_name].items():
        n_sig = (res_df['FDR q-val'] < 0.05).sum()
        pathway_report['gsea_summary'][ct_name][gs] = int(n_sig)

with open(OUTDIR + 'phase5_step2_pathway_report.json', 'w') as f:
    json.dump(pathway_report, f, indent=2, default=str)
print("  phase5_step2_pathway_report.json")

print("\n" + "=" * 60)
print("Phase 5 Step 2 Complete!")
print(f"  GSEA: {len(gsea_results)} cell types")
print(f"  ORA: {len(ora_results)} cell types")
print("=" * 60)
