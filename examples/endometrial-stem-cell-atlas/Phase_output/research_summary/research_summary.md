
## Human Endometrial Basalis Stem Cell Atlas and Drug Intervention Targets for Regenerative Disorders

### Key Findings

1. **High-Resolution Stem Cell Atlas**: Integrated 5 scRNA-seq datasets (314,805 cells) to construct a comprehensive human endometrial atlas, identifying 14 cell types including two key stem/progenitor subpopulations:
   - SOX9+LGR5+_stem (9,489 cells, 3.0%): Highest stemness stem cell population
   - CD133+_progenitor (11,466 cells, 3.6%): PROM1+ marked progenitor cells

2. **Stem Cell Lineage Hierarchy**: DPT trajectory analysis reveals two differentiation lineages:
   - Epithelial lineage: SOX9+LGR5+_stem → CD133+_progenitor → Glandular/Secretory epithelium
   - Mesenchymal lineage: Basalis_fibroblast → Decidualized_stromal → Smooth_muscle
   CytoTRACE scoring confirms SOX9+LGR5+ as the highest stemness population (0.342)

3. **Asherman's Syndrome Stem Cell Depletion**: Systematic stem cell dysfunction discovered in AS:
   - CD133+ progenitor cell count dramatically reduced (7.7%→0.68%, padj=0.011)
   - SOX9+ stem cell stemness dramatically decreased (CytoTRACE 0.334→0.072, p=1.78e-48)
   - All key signaling pathways (WNT/NOTCH/TGFb/FGF/EGF/Hedgehog) globally downregulated in SOX9+ stem cells
   - Basalis_fibroblast inversely upregulated (fibrotic reprogramming)

4. **Niche Microenvironment Remodeling**: Fundamental changes in stem cell niche communication network in AS:
   - Net enhancement of fibrotic signals (2,167↑ vs 1,433↓)
   - ECM/Adhesion pathways dominant (299 significant interactions)
   - Basalis_fibroblast transformed from niche supporter to fibrosis driver

5. **Drug Target Discovery**: Integrated DEG/TF/communication three-layer evidence to identify 463 candidate targets:
   - Top targets: IER3, ELF3, CXCL8, MET, PGR, etc.
   - 60/100 targets have known drug interactions, 24 have approved drugs
   - MET (cabozantinib, etc.) and PGR (progesterone, etc.) have highest translational potential
   - Organoid validation: 30/30 targets detectable in organoids

### Spatial Validation
RIF Visium spatial analysis validated stem cell enrichment in the Basalis_niche region (SOX9+ 4.31%, CD133+ 9-13%). Cell2location deconvolution shows reduced SOX9+ stem cells in RIF (CTR=3.26% vs RIF=2.48%, p<0.001).


## Figure Plan

### Fig1: Comprehensive Human Endometrial Atlas Construction
- A: Analysis workflow diagram (5 datasets → 314,805 cells → 14 cell types)
- B: UMAP by celltype_manual (14 cell types, Nature-style color scheme)
- C: UMAP by dataset (5 colors)
- D: Marker dotplot (top 3 markers × 14 types)
- E: Cell type proportion bar plot (by dataset)
- F: Lineage distribution pie chart (Mesenchymal/Immune/Epithelial/Stem/Endothelial)

### Fig2: Stem/Progenitor Subpopulation Identification and Molecular Characterization
- A: Stem cell subpopulation UMAP zoom-in (SOX9+LGR5+, CD133+, Basalis_fib)
- B: SOX9/LGR5/PROM1 feature plot on UMAP
- C: Violin plot: SOX9, LGR5, PROM1, CD44, KRT19, AXIN2 by stem types
- D: Stem cell marker heatmap (stem vs non-stem)
- E: CytoTRACE score UMAP
- F: CytoTRACE score violin by celltype

### Fig3: Stem Cell Spatial Localization and Differentiation Trajectory
- A: Visium spatial region annotation (Stroma/Functionalis/Basalis_niche)
- B: Cell2location deconvolution — SOX9+LGR5+ spatial distribution
- C: Cell2location deconvolution — CD133+_progenitor spatial distribution
- D: PAGA connectivity graph (epithelial lineage)
- E: DPT pseudotime UMAP (epithelial lineage)
- F: DPT distribution violin by celltype

### Fig4: Stem Cell Depletion and Dysfunction in Asherman's Syndrome
- A: Cell proportion changes bar plot (Control vs AS, significance annotated)
- B: CytoTRACE Control vs AS violin (SOX9+, CD133+)
- C: DPT distribution Control vs AS (KS test)
- D: Volcano plot — CD133+ DEGs
- E: Volcano plot — Decidualized_stromal DEGs
- F: GSEA custom pathways heatmap (WNT/NOTCH/TGFb/Fibrosis × celltypes)

### Fig5: Transcription Factor Regulation and Cell Communication Network Remodeling
- A: TF activity heatmap (top variable TFs × celltypes)
- B: Disease differential TF heatmap (stem cell key TFs, AS vs Control)
- C: Communication interaction count heatmap (sender × receiver)
- D: Stem cell pathway score heatmap (8 pathways × 3 stem types)
- E: Pathway disease change heatmap (AS vs Control, padj annotated)
- F: Fibrosis vs regeneration signal distribution plot

### Fig6: Drug Target Prioritization and Druggability Assessment
- A: Top30 target priority scores (three-color Tier annotation)
- B: Multi-evidence support heatmap (DEG/TF/Comm/Drug/Organoid)
- C: Drug-target network diagram (top approved drugs)
- D: Target-pathway association bar plot
- E: Organoid expression validation (top15 targets)
- F: Core mechanism model diagram (schematic)


### Supplementary Figures
- S1: QC metrics across datasets
- S2: Batch integration quality (ASW, UMAP before/after)
- S3: HECA reference mapping confidence
- S4: RIF Visium deconvolution all cell types
- S5: Mesenchymal trajectory PAGA/DPT
- S6: Full DEG heatmap all cell types
- S7: Complete GSEA results
- S8: TF activity UMAP for key TFs
- S9: All pathway communication results
- S10: Full drug target list

### Supplementary Tables
- ST1: Dataset metadata and QC summary
- ST2: Cell type marker genes
- ST3: DEG full results per cell type
- ST4: GSEA/ORA pathway enrichment results
- ST5: TF activity differential results
- ST6: Cell communication results (liana)
- ST7: Drug target prioritization full list
- ST8: DGIdb drug-target interactions
