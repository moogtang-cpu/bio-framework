#!/usr/bin/env python
"""
Phase 8: Drug Target Discovery and Validation
Integrates differential analysis results from Phase 5-7 for target prioritization,
drug matching, and evidence chain construction
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import requests
import json
import os
import gc
import time
import warnings
warnings.filterwarnings('ignore')

FIGDIR = 'Phase_output/phase8_drug_targets/figures/'
OUTDIR = 'Phase_output/phase8_drug_targets/'

print("=" * 60)
print("Phase 8: Drug Target Discovery and Validation")
print("=" * 60)

# ============================================================
# 1. Consolidate multi-phase evidence
# ============================================================
print("\n--- 1. Consolidate multi-phase evidence ---")

# 1a. Phase 5 DEG results
deg_files = [f for f in os.listdir('Phase_output/phase5_disease/')
             if f.startswith('deg_') and f.endswith('.csv')]

all_degs = {}
for f in deg_files:
    ct = f.replace('deg_', '').replace('_AS_vs_WOI.csv', '')
    df = pd.read_csv(f'Phase_output/phase5_disease/{f}', index_col=0)
    df['gene'] = df.index
    sig = df[(df['padj'] < 0.05) & (df['log2FoldChange'].abs() > 0.25)]
    all_degs[ct] = sig
    print(f"  DEG {ct}: {len(sig)} ({(sig['log2FoldChange']>0).sum()}↑ + {(sig['log2FoldChange']<0).sum()}↓)")

# 1b. Phase 6 TF differential results
tf_files = [f for f in os.listdir('Phase_output/phase6_trajectory/')
            if f.startswith('tf_diff_') and f.endswith('.csv')]

all_tf_diff = {}
for f in tf_files:
    ct = f.replace('tf_diff_', '').replace('_AS_vs_Control.csv', '')
    df = pd.read_csv(f'Phase_output/phase6_trajectory/{f}')
    sig = df[df['padj'] < 0.05]
    all_tf_diff[ct] = sig
    print(f"  TF diff {ct}: {len(sig)}")

# 1c. Phase 7 communication differences
comm_targets = json.load(open('Phase_output/phase7_communication/drug_target_candidates_from_comm.json'))
print(f"  Communication candidate ligands: {len(comm_targets['top_changed_ligands'])}")
print(f"  Communication candidate receptors: {len(comm_targets['top_changed_receptors'])}")

# 1d. Phase 7 pathway changes
pw_diff = pd.read_csv('Phase_output/phase7_communication/stem_pathway_disease_diff.csv')
sig_pw = pw_diff[pw_diff['padj'].astype(float) < 0.05]
print(f"  Significant pathway changes: {len(sig_pw)}")

# ============================================================
# 2. Target priority scoring
# ============================================================
print("\n--- 2. Target priority scoring ---")

# Target score = DEG evidence + TF regulation evidence + communication evidence + stemness association + pathway association
# Focus on stem cell-related genes

# Collect all candidate targets
candidate_genes = set()

# From stem cell DEGs
for ct in ['CD133plus_progenitor', 'Basalis_fibroblast']:
    if ct in all_degs:
        candidate_genes.update(all_degs[ct].head(100)['gene'].tolist())

# From decidualized stromal DEGs (largest DEG count)
if 'Decidualized_stromal' in all_degs:
    candidate_genes.update(all_degs['Decidualized_stromal'].head(200)['gene'].tolist())

# From TF differentials
for ct_key, tf_df in all_tf_diff.items():
    top_tfs = tf_df.sort_values('padj').head(30)['TF'].tolist()
    candidate_genes.update(top_tfs)

# From communication candidates
candidate_genes.update(comm_targets['top_changed_ligands'])
candidate_genes.update(comm_targets['top_changed_receptors'])

# Filter complex names (those containing _ are generally complexes)
candidate_genes = {g for g in candidate_genes if '_' not in g and pd.notna(g)}
print(f"  Total candidate targets: {len(candidate_genes)}")

# Scoring system
scores = {}
for gene in candidate_genes:
    score = 0
    evidence = []

    # DEG evidence (number of cell types with differential expression)
    deg_count = 0
    max_lfc = 0
    for ct, deg_df in all_degs.items():
        match = deg_df[deg_df['gene'] == gene]
        if len(match) > 0:
            deg_count += 1
            lfc = abs(match.iloc[0]['log2FoldChange'])
            if lfc > max_lfc:
                max_lfc = lfc
    if deg_count > 0:
        score += deg_count * 2 + min(max_lfc, 5)
        evidence.append(f"DEG_in_{deg_count}ct(maxLFC={max_lfc:.2f})")

    # Bonus for stem cell type-specific DEGs
    for ct in ['CD133plus_progenitor', 'Basalis_fibroblast']:
        if ct in all_degs:
            match = all_degs[ct][all_degs[ct]['gene'] == gene]
            if len(match) > 0:
                score += 3  # Extra points for stem cell DEGs
                evidence.append(f"stem_DEG_{ct}")

    # TF regulation evidence
    tf_count = 0
    for ct_key, tf_df in all_tf_diff.items():
        match = tf_df[tf_df['TF'] == gene]
        if len(match) > 0:
            tf_count += 1
    if tf_count > 0:
        score += tf_count * 2
        evidence.append(f"TF_diff_{tf_count}ct")

    # Communication evidence
    if gene in comm_targets['top_changed_ligands']:
        score += 3
        evidence.append("comm_ligand")
    if gene in comm_targets['top_changed_receptors']:
        score += 3
        evidence.append("comm_receptor")

    # Known stem cell / regeneration / fibrosis pathway association
    stem_genes = {'SOX9', 'LGR5', 'PROM1', 'CD44', 'AXIN2', 'WNT2', 'WNT3A', 'WNT4',
                  'WNT5A', 'WNT7A', 'RSPO1', 'RSPO3', 'NOTCH1', 'NOTCH2', 'DLL1', 'JAG1',
                  'HES1', 'BMP2', 'BMP4', 'TGFB1', 'TGFB3', 'FGF2', 'FGF7', 'FGF10',
                  'EGF', 'EGFR', 'PDGFRA', 'PDGFRB', 'MET', 'HGF', 'VEGFA', 'VEGFR2',
                  'MYC', 'TP53', 'FOXO1', 'STAT3', 'NFKB1', 'JUN', 'FOS'}
    if gene in stem_genes:
        score += 2
        evidence.append("known_stem_gene")

    fibrosis_genes = {'COL1A1', 'COL1A2', 'COL3A1', 'FN1', 'ACTA2', 'TGFB1', 'TGFB2',
                      'CTGF', 'SERPINE1', 'MMP2', 'MMP9', 'MMP14', 'TIMP1', 'TIMP2',
                      'LOXL2', 'LOX', 'PDGFA', 'PDGFB', 'IL6', 'IL1B', 'TNF', 'CCL2'}
    if gene in fibrosis_genes:
        score += 2
        evidence.append("fibrosis_gene")

    if score > 0:
        scores[gene] = {
            'score': score,
            'evidence': '; '.join(evidence),
            'n_evidence_types': len(evidence),
        }

# Ranking
scores_df = pd.DataFrame(scores).T
scores_df = scores_df.sort_values('score', ascending=False)
scores_df.index.name = 'gene'
scores_df.to_csv(OUTDIR + 'target_prioritization.csv')
print(f"  Scored targets: {len(scores_df)}")
print(f"  Top 20 targets:")
for gene, row in scores_df.head(20).iterrows():
    print(f"    {gene}: score={row['score']:.0f} ({row['evidence']})")

# ============================================================
# 3. DGIdb drug query
# ============================================================
print("\n--- 3. DGIdb drug query ---")

top_targets = scores_df.head(100).index.tolist()

# DGIdb GraphQL batch query
url = 'https://dgidb.org/api/graphql'

def query_dgidb(gene_list, batch_size=25):
    """Query DGIdb GraphQL API in batches"""
    all_results = {}
    for i in range(0, len(gene_list), batch_size):
        batch = gene_list[i:i+batch_size]
        genes_str = ', '.join(f'"{g}"' for g in batch)
        query = f'''
        {{
          genes(names: [{genes_str}]) {{
            nodes {{
              name
              interactions {{
                drug {{
                  name
                  approved
                  conceptId
                }}
                interactionScore
                interactionTypes {{
                  type
                  directionality
                }}
                publications {{
                  pmid
                }}
                sources {{
                  fullName
                }}
              }}
            }}
          }}
        }}
        '''
        try:
            r = requests.post(url, json={'query': query}, timeout=60)
            if r.status_code == 200:
                data = r.json()
                nodes = data.get('data', {}).get('genes', {}).get('nodes', [])
                for node in nodes:
                    gene_name = node['name']
                    interactions = node.get('interactions', [])
                    if interactions:
                        all_results[gene_name] = interactions
            time.sleep(1)  # Avoid rate limiting
        except Exception as e:
            print(f"  Batch {i}: Error - {e}")
    return all_results

print("  Querying DGIdb (top 100 targets)...")
dgidb_results = query_dgidb(top_targets)
print(f"  Targets with drug interactions: {len(dgidb_results)}")

# Organize drug-target relationships
drug_target_pairs = []
for gene, interactions in dgidb_results.items():
    for ix in interactions:
        drug = ix.get('drug', {})
        types = ix.get('interactionTypes', [])
        pubs = ix.get('publications', [])
        drug_target_pairs.append({
            'gene': gene,
            'drug': drug.get('name', 'Unknown'),
            'approved': drug.get('approved', False),
            'concept_id': drug.get('conceptId', ''),
            'interaction_types': ', '.join(t.get('type', '') for t in types if t.get('type')),
            'directionality': ', '.join(t.get('directionality', '') for t in types if t.get('directionality')),
            'n_publications': len(pubs),
            'score': scores_df.loc[gene, 'score'] if gene in scores_df.index else 0,
        })

drug_df = pd.DataFrame(drug_target_pairs)
if len(drug_df) > 0:
    drug_df = drug_df.sort_values(['score', 'approved'], ascending=[False, False])
    drug_df.to_csv(OUTDIR + 'drug_gene_interactions.csv', index=False)

    # Statistics
    approved_drugs = drug_df[drug_df['approved'] == True]
    unique_drugs = drug_df['drug'].nunique()
    unique_targets = drug_df['gene'].nunique()
    print(f"  Drug-target pairs: {len(drug_df)}")
    print(f"  Unique drugs: {unique_drugs}")
    print(f"  Druggable targets: {unique_targets}")
    print(f"  Approved drugs: {len(approved_drugs)} (covering {approved_drugs['gene'].nunique()} targets)")

    # Approved drugs for top targets
    print("\n  Approved drugs for top targets:")
    for gene in scores_df.head(30).index:
        gene_drugs = drug_df[(drug_df['gene'] == gene) & (drug_df['approved'] == True)]
        if len(gene_drugs) > 0:
            drug_list = gene_drugs['drug'].unique()[:5]
            print(f"    {gene} (score={scores_df.loc[gene, 'score']:.0f}): {', '.join(drug_list)}")
else:
    print("  No drug query results")

# ============================================================
# 4. Druggability assessment
# ============================================================
print("\n--- 4. Druggability assessment ---")

# Drug categories based on DGIdb
druggability = {}
for gene in scores_df.head(50).index:
    gene_info = {
        'priority_score': float(scores_df.loc[gene, 'score']),
        'evidence': scores_df.loc[gene, 'evidence'],
        'has_drug_interactions': gene in dgidb_results,
        'n_drugs': 0,
        'n_approved_drugs': 0,
        'drug_categories': [],
        'druggability_tier': 'Tier3_novel',  # default
    }

    if gene in dgidb_results:
        gene_drugs = drug_df[drug_df['gene'] == gene]
        gene_info['n_drugs'] = int(gene_drugs['drug'].nunique())
        gene_info['n_approved_drugs'] = int(gene_drugs[gene_drugs['approved'] == True]['drug'].nunique())

        # Categorize
        itypes = set()
        for _, row in gene_drugs.iterrows():
            if row['interaction_types']:
                itypes.update(row['interaction_types'].split(', '))
        gene_info['drug_categories'] = list(itypes)

        # Druggability tiering
        if gene_info['n_approved_drugs'] > 0:
            gene_info['druggability_tier'] = 'Tier1_approved'
        elif gene_info['n_drugs'] > 0:
            gene_info['druggability_tier'] = 'Tier2_investigational'

    druggability[gene] = gene_info

# Summary
tier_counts = {}
for gene, info in druggability.items():
    tier = info['druggability_tier']
    tier_counts[tier] = tier_counts.get(tier, 0) + 1
print(f"  Druggability tiers (top 50 targets):")
for tier, n in sorted(tier_counts.items()):
    print(f"    {tier}: {n}")

# ============================================================
# 5. Organoid data validation
# ============================================================
print("\n--- 5. Organoid data validation ---")

try:
    adata_org = ad.read_h5ad('data/organized/GSE216748_organoid.h5ad')
    print(f"  Organoid data: {adata_org.shape}")

    # Check expression of top targets in organoids
    top30 = scores_df.head(30).index.tolist()
    available_in_org = [g for g in top30 if g in adata_org.var_names]
    print(f"  Top 30 targets available in organoids: {len(available_in_org)}/{len(top30)}")

    if len(available_in_org) > 0 and 'log_normalized' in adata_org.layers:
        adata_org.X = adata_org.layers['log_normalized'].copy()
    elif hasattr(adata_org.X, 'toarray'):
        pass  # Use original X

    # Calculate expression levels
    from scipy.sparse import issparse
    org_expr = {}
    for gene in available_in_org:
        if gene in adata_org.var_names:
            idx = list(adata_org.var_names).index(gene)
            x = adata_org.X[:, idx]
            if issparse(x):
                x = x.toarray().flatten()
            elif hasattr(x, 'flatten'):
                x = x.flatten()
            org_expr[gene] = {
                'mean_expr': float(np.mean(x)),
                'pct_expressing': float((x > 0).mean() * 100),
            }

    # Sort by expression level
    org_df = pd.DataFrame(org_expr).T
    org_df = org_df.sort_values('mean_expr', ascending=False)
    print(f"\n  Targets with highest expression in organoids:")
    for gene, row in org_df.head(15).iterrows():
        print(f"    {gene}: mean={row['mean_expr']:.3f}, pct={row['pct_expressing']:.1f}%")

    org_df.to_csv(OUTDIR + 'organoid_validation.csv')
    del adata_org
    gc.collect()
except Exception as e:
    print(f"  Organoid validation skipped: {e}")
    org_df = pd.DataFrame()

# ============================================================
# 6. Top 10 target evidence chain construction
# ============================================================
print("\n--- 6. Top 10 target evidence chains ---")

top10_evidence = {}
for rank, (gene, row) in enumerate(scores_df.head(10).iterrows(), 1):
    chain = {
        'rank': rank,
        'gene': gene,
        'priority_score': float(row['score']),
        'evidence_summary': row['evidence'],
        'deg_evidence': {},
        'tf_evidence': {},
        'comm_evidence': {},
        'druggability': druggability.get(gene, {}),
        'organoid_expr': {},
    }

    # DEG details
    for ct, deg_df in all_degs.items():
        match = deg_df[deg_df['gene'] == gene]
        if len(match) > 0:
            chain['deg_evidence'][ct] = {
                'log2FC': float(match.iloc[0]['log2FoldChange']),
                'padj': float(match.iloc[0]['padj']),
            }

    # TF details
    for ct_key, tf_df in all_tf_diff.items():
        match = tf_df[tf_df['TF'] == gene]
        if len(match) > 0:
            chain['tf_evidence'][ct_key] = {
                'diff': float(match.iloc[0]['diff_AS_vs_Control']),
                'padj': float(match.iloc[0]['padj']),
            }

    # Communication
    if gene in comm_targets['top_changed_ligands']:
        chain['comm_evidence']['role'] = 'ligand'
    if gene in comm_targets['top_changed_receptors']:
        chain['comm_evidence']['role'] = chain['comm_evidence'].get('role', '') + '_receptor'

    # Organoid
    if gene in org_df.index:
        chain['organoid_expr'] = {
            'mean': float(org_df.loc[gene, 'mean_expr']),
            'pct': float(org_df.loc[gene, 'pct_expressing']),
        }

    top10_evidence[gene] = chain
    print(f"\n  #{rank} {gene} (score={row['score']:.0f}):")
    print(f"    Evidence: {row['evidence']}")
    if chain['druggability'].get('n_approved_drugs', 0) > 0:
        print(f"    Approved drugs: {chain['druggability']['n_approved_drugs']}")
    if chain['organoid_expr']:
        print(f"    Organoid expression: mean={chain['organoid_expr']['mean']:.3f}, pct={chain['organoid_expr']['pct']:.1f}%")

# ============================================================
# 7. Visualization
# ============================================================
print("\n--- 7. Visualization ---")

# 7a. Target priority score bar chart (top 30)
print("  7a. Target priority score plot...")
top30_scores = scores_df.head(30).copy()
top30_scores['tier'] = [druggability.get(g, {}).get('druggability_tier', 'unknown') for g in top30_scores.index]
colors = {'Tier1_approved': '#2ca02c', 'Tier2_investigational': '#ff7f0e', 'Tier3_novel': '#d62728', 'unknown': '#7f7f7f'}
bar_colors = [colors.get(t, '#7f7f7f') for t in top30_scores['tier']]

fig, ax = plt.subplots(figsize=(10, 8))
y_pos = range(len(top30_scores))
ax.barh(y_pos, top30_scores['score'], color=bar_colors, alpha=0.8, edgecolor='black', linewidth=0.3)
ax.set_yticks(y_pos)
ax.set_yticklabels(top30_scores.index)
ax.invert_yaxis()
ax.set_xlabel('Priority Score')
ax.set_title('Top 30 Drug Target Candidates for AS\n(Green=Approved drugs, Orange=Investigational, Red=Novel)')

# Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=colors['Tier1_approved'], label='Tier1: Approved drugs'),
                   Patch(facecolor=colors['Tier2_investigational'], label='Tier2: Investigational'),
                   Patch(facecolor=colors['Tier3_novel'], label='Tier3: Novel target')]
ax.legend(handles=legend_elements, loc='lower right')
plt.tight_layout()
plt.savefig(FIGDIR + 'target_priority_scores.png', dpi=200, bbox_inches='tight')
plt.close()
print("    target_priority_scores.png")

# 7b. Drug-target network heatmap
print("  7b. Drug-target network...")
if len(drug_df) > 0:
    # Top targets × approved drugs
    top_genes = scores_df.head(20).index.tolist()
    top_gene_drugs = drug_df[(drug_df['gene'].isin(top_genes)) & (drug_df['approved'] == True)]

    if len(top_gene_drugs) > 0:
        # Get most frequent drugs
        drug_counts = top_gene_drugs['drug'].value_counts()
        top_drugs = drug_counts.head(30).index.tolist()

        # Create gene-drug matrix
        pivot = top_gene_drugs[top_gene_drugs['drug'].isin(top_drugs)].groupby(
            ['gene', 'drug']).size().reset_index(name='count')
        heatmap_data = pivot.pivot(index='gene', columns='drug', values='count').fillna(0)

        if heatmap_data.shape[0] > 1 and heatmap_data.shape[1] > 1:
            fig, ax = plt.subplots(figsize=(max(12, len(top_drugs)*0.5), max(6, len(heatmap_data)*0.4)))
            sns.heatmap(heatmap_data, cmap='YlOrRd', linewidths=0.3, ax=ax,
                        cbar_kws={'label': 'Interaction Evidence'})
            ax.set_title('Drug-Target Interaction Network (Approved Drugs)')
            ax.set_ylabel('Target Gene')
            ax.set_xlabel('')
            plt.xticks(rotation=45, ha='right', fontsize=8)
            plt.tight_layout()
            plt.savefig(FIGDIR + 'drug_target_network.png', dpi=200, bbox_inches='tight')
            plt.close()
            print("    drug_target_network.png")

# 7c. Evidence type heatmap
print("  7c. Evidence heatmap...")
evidence_matrix = []
for gene in scores_df.head(25).index:
    ev = {
        'gene': gene,
        'DEG_celltypes': sum(1 for ct in all_degs if gene in all_degs[ct]['gene'].values),
        'TF_diff_celltypes': sum(1 for ct in all_tf_diff if gene in all_tf_diff[ct]['TF'].values),
        'comm_ligand': 1 if gene in comm_targets['top_changed_ligands'] else 0,
        'comm_receptor': 1 if gene in comm_targets['top_changed_receptors'] else 0,
        'has_approved_drug': 1 if druggability.get(gene, {}).get('n_approved_drugs', 0) > 0 else 0,
        'organoid_expressed': 1 if gene in org_df.index and org_df.loc[gene, 'pct_expressing'] > 10 else 0,
    }
    evidence_matrix.append(ev)

ev_df = pd.DataFrame(evidence_matrix).set_index('gene')
fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(ev_df, cmap='YlGnBu', annot=True, fmt='.0f', linewidths=0.3, ax=ax,
            cbar_kws={'label': 'Evidence Score'})
ax.set_title('Multi-Evidence Support for Top Drug Targets')
ax.set_ylabel('')
plt.tight_layout()
plt.savefig(FIGDIR + 'evidence_heatmap.png', dpi=200, bbox_inches='tight')
plt.close()
print("    evidence_heatmap.png")

# 7d. Target-pathway-phenotype association network (Sankey-like visualization)
print("  7d. Target-pathway association...")
# Simplified as grouped bar chart
pathway_links = {
    'WNT_signaling': ['WNT2', 'WNT3A', 'WNT4', 'WNT5A', 'WNT7A', 'RSPO1', 'RSPO3',
                      'FZD1', 'FZD5', 'LRP5', 'LRP6', 'LGR5', 'DKK1', 'DKK3', 'SFRP1', 'SFRP2'],
    'NOTCH_signaling': ['NOTCH1', 'NOTCH2', 'DLL1', 'DLL4', 'JAG1', 'JAG2', 'HES1'],
    'TGFb_fibrosis': ['TGFB1', 'TGFB2', 'TGFB3', 'BMP2', 'BMP4', 'COL1A1', 'COL3A1',
                      'FN1', 'ACTA2', 'SERPINE1', 'CTGF', 'LOX', 'LOXL2'],
    'Growth_factors': ['FGF2', 'FGF7', 'FGF10', 'EGF', 'EGFR', 'PDGFA', 'PDGFRA',
                       'HGF', 'MET', 'VEGFA', 'IGF1'],
    'Inflammation': ['IL6', 'IL1B', 'TNF', 'CXCL8', 'CXCL12', 'CCL2', 'NFKB1'],
    'Transcription_factors': ['SOX9', 'MYC', 'TP53', 'FOXO1', 'STAT3', 'JUN',
                              'REST', 'ZEB2', 'SNAI2', 'GRHL2'],
}

top50_genes = set(scores_df.head(50).index)
pw_target_counts = {}
for pw, genes in pathway_links.items():
    overlap = top50_genes & set(genes)
    if overlap:
        pw_target_counts[pw] = len(overlap)

fig, ax = plt.subplots(figsize=(10, 5))
pws = sorted(pw_target_counts.keys(), key=lambda x: pw_target_counts[x], reverse=True)
counts = [pw_target_counts[pw] for pw in pws]
bars = ax.barh(pws, counts, color=plt.cm.Set2(range(len(pws))), edgecolor='black', linewidth=0.3)
ax.set_xlabel('Number of Top50 Targets')
ax.set_title('Drug Targets by Signaling Pathway')
for bar, c in zip(bars, counts):
    ax.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2, str(c),
            va='center', fontsize=10)
plt.tight_layout()
plt.savefig(FIGDIR + 'targets_by_pathway.png', dpi=200, bbox_inches='tight')
plt.close()
print("    targets_by_pathway.png")

# ============================================================
# 8. Save complete report
# ============================================================
print("\n--- 8. Save report ---")

report = {
    'phase': 8,
    'target_prioritization': {
        'n_candidates': int(len(scores_df)),
        'scoring_criteria': 'DEG evidence + TF regulation + cell communication + stem gene + fibrosis gene',
        'top10': [
            {'rank': i+1, 'gene': g, 'score': float(scores_df.loc[g, 'score']),
             'evidence': scores_df.loc[g, 'evidence']}
            for i, g in enumerate(scores_df.head(10).index)
        ],
    },
    'drug_database': {
        'source': 'DGIdb GraphQL API',
        'n_targets_with_drugs': len(dgidb_results),
        'n_drug_target_pairs': len(drug_df) if len(drug_df) > 0 else 0,
        'n_approved_drugs': int(drug_df[drug_df['approved'] == True]['drug'].nunique()) if len(drug_df) > 0 else 0,
    },
    'druggability': {
        tier: count for tier, count in sorted(tier_counts.items())
    },
    'organoid_validation': {
        'n_targets_expressed': int(len(org_df[org_df['pct_expressing'] > 10])) if len(org_df) > 0 else 0,
        'top5_expressed': org_df.head(5).index.tolist() if len(org_df) > 0 else [],
    },
    'top10_evidence_chains': top10_evidence,
}

with open(OUTDIR + 'phase8_report.json', 'w') as f:
    json.dump(report, f, indent=2, default=str)
print("  phase8_report.json")

# Save complete target list (CSV)
full_target_table = scores_df.copy()
full_target_table['druggability_tier'] = [druggability.get(g, {}).get('druggability_tier', 'unknown')
                                          for g in full_target_table.index]
full_target_table['n_approved_drugs'] = [druggability.get(g, {}).get('n_approved_drugs', 0)
                                         for g in full_target_table.index]
full_target_table.to_csv(OUTDIR + 'full_target_list.csv')
print("  full_target_list.csv")

print("\n" + "=" * 60)
print("Phase 8 complete!")
print(f"  Candidate targets: {len(scores_df)}")
print(f"  With drug interactions: {len(dgidb_results)}")
print(f"  Top 10 targets: {', '.join(scores_df.head(10).index.tolist())}")
print("=" * 60)
