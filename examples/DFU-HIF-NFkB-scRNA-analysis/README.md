# HIF-1α/NF-κB Hypoxia-Inflammation Coupling in Diabetic Foot Ulcers

Integrative single-cell and bulk RNA-seq analysis investigating the HIF-1α/NF-κB coupling axis in diabetic foot ulcer (DFU) wound healing, with focus on macrophage polarization and hypoxia signaling.

## Scientific Questions

| # | Question | Addressed by |
|---|----------|-------------|
| Q1 | Which cell types are the primary carriers of HIF1A/NFKB1 co-expression? How do activation states differ across cell types? | Phase 1 |
| Q2 | How does HIF1A/NFKB1 co-expression correlate with wound healing outcomes? | Phase 1, 2 |
| Q3 | What are the common target genes co-regulated by HIF-1α and NF-κB? How are they dysregulated in diabetic wounds? | Phase 2, 4 |
| Q4 | Is HIF/NF-κB dysregulation a sustained or delayed response? How does it differ from normal acute healing? | Phase 3 |
| Q5 | How does HIF/NF-κB activation in macrophage subpopulations relate to M1/M2 polarization? | Phase 1, 4 |

## Key Findings

1. **Non-healing DFUs show M2 macrophage retention** (challenges conventional model): Healing wounds show M1-biased macrophages; non-healing wounds are locked in a dysregulated M2-dominant state, suggesting failure of normal M1→M2 dynamic transition
2. **NF-κB is the primary M1 driver** (rho = 0.728), with HIF-1α as synergistic co-activator (rho = 0.594)
3. **HIF+M1 macrophage subcluster**: 42.1% in healing vs. 9.9% in non-healing wounds
4. **HE-Fibroblasts** (MMP1+/MMP3+/HIF1A+): 4.2-fold enriched in healing DFU
5. **26 HIF+NF-κB co-regulated target genes** identified, 7 validated across ≥3 layers (VEGFA, IL6, IL1B, PTGS2, BNIP3, ADM, ANGPTL4)
6. **Hypoxia score**: significantly higher in healing wounds (Macrophage p = 2.9e-194)

## Datasets

All raw data are publicly available from GEO. Download and place them under `data/` before running the scripts.

| GEO ID | Type | Samples | Description | Used in |
|--------|------|---------|-------------|---------|
| [GSE165816](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165816) | scRNA-seq | 50 | DFU foot/forearm skin + PBMCs (Healing/Non-healing/Diabetic/Healthy) | Phase 1, 4 |
| [GSE134431](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134431) | Bulk RNA-seq | 21 | DFU patients vs. diabetic foot skin | Phase 2 |
| [GSE199939](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199939) | Bulk RNA-seq | 21 | DFU vs. non-diabetic foot skin (TPM) | Phase 2 |
| [GSE80178](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80178) | Microarray | 12 | DFU vs. diabetic/non-diabetic foot skin | Phase 2 |
| [GSE28914](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28914) | Microarray | 25 | Normal skin re-epithelialization time series (D0/D3/D7) | Phase 3 |
| [GSE50425](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50425) | Microarray | 12 | Late-stage skin graft donor site healing (D0/D14/D21) | Phase 3 |
| [GSE147890](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147890) | Microarray | 26 | Humanized skin mouse model (Control/Diabetic, 0h/24h) | Phase 3 |
| [GSE15949](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15949) | Microarray | 2 | Macrophage hypoxia response (<0.5% O2) | Phase 4 |
| [GSE16099](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16099) | Microarray | 4 | Macrophage hypoxia 18h (0.1% O2) | Phase 4 |

## Analysis Pipeline

```
Phase 1: Single-Cell Core Analysis (GSE165816)
  ├── QC & filtering → Integration (Harmony) → Clustering
  ├── Cell type annotation (6 major types)
  ├── HIF1A/NFKB1 expression & co-expression analysis
  ├── Pathway scoring (HIF, NF-κB, M1/M2)
  ├── Healing vs. Non-healing differential analysis
  ├── Fibroblast subclustering (HE-Fibroblast identification)
  ├── Macrophage subclustering (M1/M2 polarization)
  └── Cell-cell communication (CellChat)

Phase 2: Bulk RNA-seq Multi-Cohort Validation (GSE134431, GSE199939, GSE80178)
  ├── Target gene expression validation
  ├── WGCNA co-expression network
  ├── Pathway enrichment
  └── Immune infiltration estimation

Phase 3: Temporal Dynamics (GSE28914, GSE50425, GSE147890)
  ├── Normal wound healing time course (D0→D3→D7→D14→D21)
  ├── Diabetic vs. control wound response (0h→24h)
  └── HIF/NF-κB temporal activation patterns

Phase 4: Macrophage Mechanism Dissection (GSE15949, GSE16099, GSE165816)
  ├── Hypoxia-responsive gene signature construction
  ├── HIF-1α/NF-κB co-regulated target gene network
  ├── Macrophage hypoxia scoring in scRNA-seq
  ├── Subcluster-level polarization-hypoxia coupling
  └── Pseudotime trajectory analysis
```

## Directory Structure

```
DFU-HIF-NFkB-scRNA-analysis/
├── README.md
├── data/
│   └── data_info.txt                        # Dataset descriptions and sample groupings
├── docs/
│   └── bio_analyze.md                       # Analysis design document
├── project_skills/                          # Topic configuration and workflow definitions
│   ├── TOPIC.yml                            # Research topic, scientific questions, dataset specs
│   ├── phases/phase_index.yml               # 4-phase analysis plan with substeps
│   └── runtime/
│       ├── phase_summaries.yml              # Per-phase key findings
│       └── workflow_state.yml               # Execution state tracker
├── scripts/                                 # R analysis scripts (run sequentially)
│   ├── phase1_1_qc.R                       # scRNA-seq QC and filtering
│   ├── phase1_1b_integrate.R               # Normalization, Harmony integration, clustering
│   ├── phase1_2_annotation.R               # Cell type annotation (marker-based)
│   ├── phase1_2b_apply_annotation.R        # Refined annotation application
│   ├── phase1_3_5_expression_pathway.R     # Expression analysis + pathway scoring + DE
│   ├── phase1_6_8_subpop.R                 # Fibroblast/macrophage subclustering + CellChat
│   ├── phase2_bulk.R                       # Bulk RNA-seq validation + WGCNA
│   ├── phase3_temporal.R                   # Temporal dynamics analysis
│   ├── phase4_macrophage.R                 # Macrophage mechanism dissection
│   ├── step4_figures.R                     # Publication figure generation
│   └── step4_figures_fix.R                 # Additional figures (Fig3, 5, 7 + supplements)
└── Phase_output/                            # All analysis results
    ├── Phase1/                              # scRNA-seq results
    │   ├── QC/                              # Quality control metrics and plots
    │   ├── clustering/                      # UMAP plots, markers, cell type tables
    │   ├── expression/                      # HIF1A/NFKB1 co-expression analysis
    │   ├── pathway/                         # Pathway scores (HIF, NF-κB, M1/M2)
    │   ├── differential/                    # DEGs per cell type (Healing vs. Non-healing)
    │   ├── fibroblast/                      # HE-Fibroblast subclustering
    │   ├── macrophage/                      # M1/M2 polarization analysis
    │   └── cellchat/                        # Ligand-receptor communication
    ├── Phase2/                              # Bulk RNA-seq validation
    ├── Phase3/                              # Temporal dynamics DEGs
    ├── Phase4/                              # Macrophage hypoxia signatures
    ├── publication_figures/                  # Final figures (Fig1-7, FigS1-S4, TableS1-S5)
    ├── final_report/technical_report.md     # Comprehensive technical report
    ├── manuscript/manuscript_draft.md       # Manuscript draft
    └── research_summary/figure_plan.yml    # Figure design plan
```

## Requirements

- **R** >= 4.4
- **Key R packages**: Seurat (v5), harmony, limma, WGCNA, CellChat, slingshot, ComplexHeatmap, pheatmap, ggplot2, dplyr, openxlsx
- **Platform annotation**: GPL5104 (for GSE16099, auto-downloaded)

## Reproduction

1. Download raw data from GEO (see Datasets table) and place files under `data/` following the structure in `data/data_info.txt`
2. Run scripts sequentially:
   ```bash
   Rscript scripts/phase1_1_qc.R
   Rscript scripts/phase1_1b_integrate.R
   Rscript scripts/phase1_2_annotation.R
   Rscript scripts/phase1_2b_apply_annotation.R
   Rscript scripts/phase1_3_5_expression_pathway.R
   Rscript scripts/phase1_6_8_subpop.R
   Rscript scripts/phase2_bulk.R
   Rscript scripts/phase3_temporal.R
   Rscript scripts/phase4_macrophage.R
   Rscript scripts/step4_figures.R
   Rscript scripts/step4_figures_fix.R
   ```
3. Results are generated in `Phase_output/`

## Notes

- The `Phase_output/` directory contains pre-computed results (CSV tables, PNG/PDF figures). Large intermediate files (Seurat `.rds` objects, ~7GB) are excluded from this repository but will be regenerated by the scripts.
- scRNA-seq data (GSE165816) requires ~16GB RAM for integration; WGCNA in Phase 2 benefits from multi-core processing.
- `data_info.txt` contains detailed sample grouping information for all 9 datasets.

## License

This analysis is provided for research and educational purposes.
