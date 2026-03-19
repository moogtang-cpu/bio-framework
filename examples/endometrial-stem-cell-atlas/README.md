# Endometrial Stem Cell Atlas

> **Human Endometrial Basalis Stem Cell Atlas Construction and Drug Intervention Target Discovery for Regenerative Disorders**

## Project Overview

- **Research Type**: Computational biology / Basic research
- **Organism**: Homo sapiens
- **Tissue**: Endometrium (basalis layer focus)
- **Disease**: Asherman's syndrome (intrauterine adhesions), thin endometrium
- **Scale**: 314,805 cells from 5 integrated scRNA-seq datasets + 8 Visium spatial samples
- **Target Journal**: Cell Reports

## Analysis Pipeline (8 Phases)

| Phase | Name | Key Methods |
|-------|------|-------------|
| 1 | Environment Setup & Data Preparation | Conda, GEO/ArrayExpress download, h5ad conversion |
| 2 | Quality Control & Preprocessing | Scrublet, mt% filtering, HVG selection, PCA |
| 3 | Multi-Dataset Integration & Atlas Construction | Harmony (GPU), Leiden clustering, HECA reference mapping, manual annotation |
| 4 | Spatial Transcriptomics Integration | Visium QC, Cell2location deconvolution |
| 5 | Disease Comparative Analysis | PyDESeq2 pseudobulk DEG, GSEA (GO/KEGG/Reactome/WikiPathway) |
| 6 | Stem Cell Trajectory & Regulatory Networks | DPT, PAGA, CytoTRACE, decoupler ULM + CollecTRI |
| 7 | Cell-Cell Communication & Microenvironment | LIANA rank_aggregate, pathway scoring |
| 8 | Drug Target Discovery & Validation | Multi-evidence scoring, DGIdb, organoid validation |

## Key Findings

1. **14 cell types** identified including **SOX9+LGR5+ stem cells** (3.0%) and **CD133+ progenitors** (3.6%)
2. **Hierarchical stem cell lineage**: SOX9+LGR5+ → CD133+ → Glandular/Secretory epithelium
3. **AS stem cell depletion**: CD133+ progenitors ~11-fold reduced (padj=0.011); SOX9+ stemness dramatically decreased (CytoTRACE 0.334→0.072)
4. **Niche remodeling**: Basalis fibroblasts transformed from niche supporters to fibrosis drivers
5. **463 drug target candidates** with 60 having known drug interactions; MET (21 approved drugs), CXCL8 (43 approved drugs) as top druggable targets

## Directory Structure

```
endometrial-stem-cell-atlas/
├── README.md
├── project_skills/
│   ├── TOPIC.yml                  # Project topic definition (questions, hypotheses, datasets)
│   ├── phases/
│   │   └── phase_index.yml        # 8-phase workflow definition
│   └── runtime/
│       ├── phase_summaries.yml    # Detailed summary of each phase
│       ├── workflow_state.yml     # Workflow progress state
│       └── findings.yml           # Issues discovered and resolutions
├── env/
│   └── environment.yml            # Conda environment specification
├── scripts/                       # Analysis scripts (Phase 1-3)
│   ├── organize_data.py
│   ├── validate_data.py
│   ├── phase2_qc_metrics.py
│   ├── phase2_doublet.py
│   ├── phase2_normalize_pca.py
│   └── phase3_batch_integration*.py
└── Phase_output/                  # Analysis outputs (Phase 3-8 + post-analysis)
    ├── phase3_integration/        # Integration results + figures
    ├── phase4_spatial/            # Spatial deconvolution results
    ├── phase5_disease/            # DEG tables, GSEA results
    ├── phase6_trajectory/         # Trajectory + TF activity results
    ├── phase7_communication/      # LIANA cell communication results
    ├── phase8_drug_targets/       # Drug target prioritization
    ├── publication_figures/        # 6 publication-ready figures (PNG + PDF)
    ├── manuscript/                # Manuscript draft
    ├── research_summary/          # Research summary and figure plan
    └── final_report/              # Technical report
```

## Datasets Used

| ID | Type | Cells/Spots | Source |
|----|------|-------------|--------|
| HECA_integrated | scRNA+snRNA+Spatial | 625,773 | Nature Genetics 2024 |
| GSE215968 | scRNA-seq | 299,351 | Nature Communications 2023 |
| GSE111976 | scRNA-seq | 73,181 | Nature Medicine 2020 |
| E-MTAB-10287 | scRNA-seq | 55,729 | Nature Genetics 2021 |
| GSE260658 | scRNA-seq | 71,058 | GEO |
| GSE287278 | 10x Visium | 10,163 | Scientific Data 2025 |
| GSE216748 | scRNA-seq (organoid) | 49,609 | Nature Communications 2023 |

> **Note**: Raw data files (.h5ad) are not included due to size. All datasets are publicly available from GEO/ArrayExpress using the accession numbers listed above.

## Publication

The manuscript draft is available at `Phase_output/manuscript/manuscript_draft.md`. Publication-ready figures are in `Phase_output/publication_figures/`.
