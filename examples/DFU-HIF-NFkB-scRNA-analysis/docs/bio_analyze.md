# HIF-1α and NF-κB Hypoxia-Inflammation Coupling Research
---

## Scientific Questions Directly Answerable by Bioinformatics

```
┌─────────────────────────────────────────────────────────────────────┐
│              Scientific Questions Answerable by Bioinformatics      │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  Q1: Which cell types are the main carriers of HIF1A/NFKB1         │
│      co-expression in diabetic wounds?                             │
│      What are the differences in activation status and             │
│      correlation of the two pathways across different cell types?  │
│                                                                     │
│  Q2: How is the co-expression pattern of HIF1A/NFKB1 associated    │
│      with wound healing outcomes?                                  │
│      In which cell types are pathway activity differences          │
│      manifested between healing vs non-healing groups?             │
│                                                                     │
│  Q3: What are the common target genes co-regulated by HIF-1α       │
│      and NF-κB?                                                    │
│      How are these common target genes dysregulated in diabetic    │
│      wounds?                                                       │
│                                                                     │
│  Q4: Is the HIF/NF-κB dysregulation in diabetic wounds sustained   │
│      or delayed in response?                                       │
│      How does it differ from the temporal dynamics pattern of      │
│      normal acute healing?                                         │
│                                                                     │
│  Q5: How is the HIF/NF-κB activation status in macrophage          │
│      subpopulations associated with M1/M2 polarization?            │
│      What are the transcriptional reprogramming features of        │
│      macrophages under hypoxic conditions?                         │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

---

---

## Bioinformatics Analysis Phase Design

```
┌─────────────────────────────────────────────────────────────────────┐
│                                                                     │
│  ══════════════════════════════════════════════════════════════    │
│  Phase 1: Single-cell Core Analysis                                │
│  ══════════════════════════════════════════════════════════════    │
│                                                                     │
│  Dataset: GSE165816 (n=50)                                         │
│  Answers questions: Q1, Q2, Q5                                     │
│                                                                     │
│  ┌─────────────────────────────────────────────────────────────┐   │
│  │ 1.1 Data Preprocessing and QC                               │   │
│  │     • Cell filtering (nFeature/nCount/percent.mt)           │   │
│  │     • Batch effect correction (Harmony/CCA)                 │   │
│  │     • Dimensionality reduction and clustering (PCA/UMAP)    │   │
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 1.2 Cell Type Annotation                                    │   │
│  │     • Marker gene identification                            │   │
│  │     • Major cell types: Fibroblasts, Macrophages,           │   │
│  │       Endothelial cells, Keratinocytes, T cells,            │   │
│  │       Neutrophils, etc.                                     │   │
│  │     • Cell proportion statistics for each group             │   │
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 1.3 HIF1A/NFKB1 Expression Analysis                         │   │
│  │     • Expression levels in each cell type (VlnPlot/DotPlot) │   │
│  │     • Proportion of HIF1A-NFKB1 co-expressing cells         │   │
│  │     • Expression correlation analysis (Pearson/Spearman)    │   │
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 1.4 Pathway Activity Scoring                                │   │
│  │     • HIF pathway score (AUCell/ssGSEA)                     │   │
│  │     • NF-κB pathway score                                   │   │
│  │     • Coupling score: HIF_score × NFKB_score correlation    │   │
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 1.5 Healing vs Non-healing Differential Analysis            │   │
│  │     • Cell proportion differences                           │   │
│  │     • Pathway activity differences                          │   │
│  │     • Differential genes in each cell type                  │   │
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 1.6 Fibroblast Subpopulation Analysis                       │   │
│  │     • Re-clustering to identify HE-Fibro subpopulation      │   │
│  │     • HIF1A⁺ feature validation (MMP1/MMP3/CHI3L1/TNFAIP6)  │   │
│  │     • HE-Fibro proportion in healing vs non-healing         │   │
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 1.7 Macrophage Analysis                                     │   │
│  │     • M1/M2 polarization scoring                            │   │
│  │     • HIF/NF-κB association with polarization status        │   │
│  │     • Validation: M1↑ in healing vs M2↑ in non-healing     │   │
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 1.8 Cell Communication Analysis                             │   │
│  │     • CellChat/CellPhoneDB                                  │   │
│  │     • Ligand-receptor interactions in HIF/NF-κB high cells  │   │
│  │     • Communication differences: healing vs non-healing     │   │
│  └─────────────────────────────────────────────────────────────┘   │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────┐
│                                                                     │
│  ══════════════════════════════════════════════════════════════    │
│  Phase 2: Bulk RNA-seq Multi-cohort Validation                     │
│  ══════════════════════════════════════════════════════════════    │
│                                                                     │
│  Datasets: GSE134431 (n=21) + GSE199939 (n=21) + GSE80178 (n=12)   │
│  Answers questions: Q2, Q3                                          │
│                                                                     │
│  ┌─────────────────────────────────────────────────────────────┐   │
│  │ 2.1 Data Preprocessing                                      │   │
│  │     • Normalization (TPM/RPKM unification)                  │   │
│  │     • Batch effect assessment and correction (ComBat/limma) │   │
│  │     • Sample grouping: DFU vs DM vs Normal                  │   │
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 2.2 Differential Expression Analysis                        │   │
│  │     • DESeq2/limma                                          │   │
│  │     • DFU vs Control                                        │   │
│  │     • Healing vs Non-healing (GSE134431)                    │   │
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 2.3 HIF1A/NFKB1 and Target Gene Expression Validation      │   │
│  │     • HIF1A, RELA, NFKB1 expression levels                 │   │
│  │     • Common target genes: TNF, IL1B, IL6, VEGFA, CXCL8    │   │
│  │     • Expression correlation analysis                       │   │
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 2.4 WGCNA Co-expression Network Analysis                   │   │
│  │     • Module construction and clustering                    │   │
│  │     • Identify modules containing HIF1A/NFKB1               │   │
│  │     • Module association with healing status                │   │
│  │     • Hub gene identification                               │   │
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 2.5 Pathway Enrichment Analysis                             │   │
│  │     • GO/KEGG/Reactome                                      │   │
│  │     • HIF-1 signaling pathway                               │   │
│  │     • NF-κB signaling pathway                               │   │
│  │     • Pathway intersection analysis                         │   │
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 2.6 Immune Infiltration Analysis                            │   │
│  │     • CIBERSORT/xCell/MCP-counter                          │   │
│  │     • Macrophage M1/M2 ratio estimation                     │   │
│  │     • Cross-validation with Phase 1 single-cell results     │   │
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 2.7 Multi-cohort Meta-analysis                              │   │
│  │     • Consistent differentially expressed gene identification│   │
│  │     • Effect size integration                               │   │
│  │     • Heterogeneity assessment                              │   │
│  └─────────────────────────────────────────────────────────────┘   │
│                                                                     │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────┐
│                                                                     │
│  ══════════════════════════════════════════════════════════════    │
│  Phase 3: Temporal Dynamics Analysis                                │
│  ══════════════════════════════════════════════════════════════    │
│                                                                     │
│  Datasets: GSE28914 + GSE50425 + GSE147890                         │
│  Answers question: Q4                                               │
│                                                                     │
│  ┌─────────────────────────────────────────────────────────────┐   │
│  │ 3.1 Human Normal Acute Healing Timeline Construction        │   │
│  │                                                             │   │
│  │     GSE28914: Day 0 → Day 1 → Day 3 → Day 7                │   │
│  │     GSE50425: Day 0 → Day 14 → Day 21                      │   │
│  │                     ↓                                       │   │
│  │     Integrated timeline: D0 → D1 → D3 → D7 → D14 → D21     │   │
│  │                                                             │   │
│  │     Analysis content:                                       │   │
│  │     • HIF1A/NFKB1 expression time curves                    │   │
│  │     • Target gene (TNF/VEGF/IL6) temporal dynamics          │   │
│  │     • Identify "normal healing pattern": Early↑ → Mid peak → Late↓│
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 3.2 Diabetic Wound Early Response Analysis                  │   │
│  │                                                             │   │
│  │     GSE147890: Humanized skin mouse model                   │   │
│  │     Comparison: Control (0h vs 24h) vs Diabetic (0h vs 24h)│   │
│  │                                                             │   │
│  │     Analysis content:                                       │   │
│  │     • Is HIF/NF-κB early response delayed/weakened under    │   │
│  │       diabetic conditions?                                  │   │
│  │     • Differential response gene identification             │   │
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 3.3 DFU-specific Dysregulated Gene Identification           │   │
│  │                                                             │   │
│  │     Strategy:                                               │   │
│  │     (DFU-specific genes) = (DFU DEGs) - (Acute healing common genes)│
│  │                                                             │   │
│  │     Steps:                                                  │   │
│  │     ① GSE28914/50425 identify acute healing response genes │   │
│  │     ② GSE165816/134431 identify DFU differential genes     │   │
│  │     ③ Take set difference → DFU-specific dysregulated genes│   │
│  │     ④ Enrichment analysis to validate gene functions       │   │
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 3.4 Temporal Pattern Comparison                             │   │
│  │                                                             │   │
│  │     Normal healing: HIF/NF-κB → Early activation → Timely resolution → Healing│
│  │     Diabetic:       HIF/NF-κB → Delayed response? Sustained activation? → Non-healing│
│  │                                                             │   │
│  │     Validated through STEM/Mfuzz time series clustering    │   │
│  └─────────────────────────────────────────────────────────────┘   │
│                                                                     │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────┐
│                                                                     │
│  ══════════════════════════════════════════════════════════════    │
│  Phase 4: Macrophage HIF/NF-κB Mechanism In-depth Analysis         │
│  ══════════════════════════════════════════════════════════════    │
│                                                                     │
│  Datasets: GSE15949 + GSE16099 + GSE165816 (macrophage subpopulations)│
│  Answers questions: Q3, Q5                                          │
│                                                                     │
│  ┌─────────────────────────────────────────────────────────────┐   │
│  │ 4.1 Hypoxia Response Gene Signature Construction            │   │
│  │                                                             │   │
│  │     GSE15949 + GSE16099:                                    │   │
│  │     • Human MDMs hypoxia (0.1% O₂, 18h) vs normoxia        │   │
│  │     • Differential gene classification:                     │   │
│  │       - HIF-dependent genes (literature-validated HIF targets)│
│  │       - NF-κB-dependent genes (literature-validated NF-κB targets)│
│  │       - Co-regulated genes (promoter contains HRE+κB sites) │   │
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 4.2 Common Target Gene Network Construction                │   │
│  │                                                             │   │
│  │     Data source integration:                                │   │
│  │     • Literature databases: HIF target genes, NF-κB target gene lists│   │
│  │     • Promoter analysis: Genes containing HRE and κB sites  │   │
│  │     • Expression data validation: GSE15949/16099 differential genes│   │
│  │                                                             │   │
│  │     Core common target genes:                               │   │
│  │     TNF, IL1B, IL6, VEGFA, CXCL8, COX2, iNOS, BCL2         │   │
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 4.3 Hypoxia Signature Application in DFU Single-cell       │   │
│  │                                                             │   │
│  │     Apply Signature constructed in 4.1 to GSE165816:       │   │
│  │     • Calculate "hypoxia response score" for each cell      │   │
│  │     • Compare score differences between healing vs non-healing groups│   │
│  │     • Identify cell subpopulations with abnormal hypoxia response│   │
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 4.4 Macrophage Subpopulation In-depth Analysis             │   │
│  │                                                             │   │
│  │     GSE165816 macrophage subset re-clustering:              │   │
│  │     • Identify specific subpopulations with high HIF/high NF-κB│   │
│  │     • M1 markers: CD86, NOS2, IL1B, TNF                    │   │
│  │     • M2 markers: CD163, MRC1, ARG1, IL10                  │   │
│  │     • Analyze association between subpopulations and healing outcomes│   │
│  ├─────────────────────────────────────────────────────────────┤   │
│  │ 4.5 Macrophage Polarization and Pathway Coupling           │   │
│  │                                                             │   │
│  │     Hypothesis validation:                                  │   │
│  │     • M1 polarization ← HIF1A↑ + NFKB↑ (pro-inflammatory) │   │
│  │     • M2 polarization ← HIF2A↑? (anti-inflammatory repair)│   │
│  │     • Healing group: Moderate M1 activation → Timely transition to M2│   │
│  │     • Non-healing group: Sustained M1 activation → M2 transition impairment│   │
│  └─────────────────────────────────────────────────────────────┘   │
│                                                                     │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```



