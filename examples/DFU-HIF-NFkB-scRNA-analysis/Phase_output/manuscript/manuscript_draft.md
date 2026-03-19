# HIF-1α/NF-κB Hypoxia-Inflammation Coupling Drives M1 Macrophage Polarization and Promotes Diabetic Foot Ulcer Healing: An Integrated Single-Cell and Bulk Transcriptomic Study

---

## Abstract

**Background:** Diabetic foot ulcers (DFUs) represent a devastating complication of diabetes mellitus, with impaired wound healing contributing to significant morbidity and mortality. The hypoxic wound microenvironment and chronic inflammation are recognized hallmarks of non-healing DFUs, yet the mechanistic interplay between hypoxia-inducible factor-1α (HIF-1α) and nuclear factor-κB (NF-κB) signaling in governing macrophage polarization and wound healing outcomes remains poorly understood.

**Methods:** We performed an integrated analysis of single-cell RNA sequencing (scRNA-seq) data from 50 samples (86,137 cells) of the GSE165816 dataset encompassing Healthy, Healing DFU, Non-Healing DFU, and Diabetic Foot Skin (DFS) conditions. Multi-cohort bulk RNA-seq validation was conducted using GSE134431 and GSE199939. Macrophage-specific hypoxia transcriptomes were interrogated using GSE15949 and GSE16099, and temporal wound healing dynamics were assessed using GSE28914, GSE50425, and GSE147890. Pathway scoring, co-expression analysis, macrophage subclustering, and pseudotime trajectory analysis were employed to dissect the HIF-1α/NF-κB axis in macrophage polarization.

**Results:** Contrary to the prevailing paradigm that M1 macrophage polarization impedes wound healing, we demonstrate that M1 polarization is essential for successful DFU resolution. Healing DFUs exhibited a pronounced M1 bias, while Non-Healing DFUs were characterized by M2 macrophage retention. NF-κB signaling served as the primary driver of M1 polarization (Spearman ρ = 0.819), with HIF-1α acting as a synergistic co-activator (ρ = 0.340). The HIF-1α/NF-κB coupling score strongly correlated with M1-M2 polarization balance (ρ = 0.721). A distinct HIF⁺M1 macrophage subcluster (C0) was identified, with 75.9% of cells originating from Healing DFUs and exhibiting markedly elevated hypoxia scores compared to Non-Healing counterparts (P = 1.33 × 10⁻¹⁵⁵). Myeloid cells harbored the highest HIF-1α/NF-κB co-expression rate (34.3% dual-positive). A novel HE-Fibroblast subpopulation (MMP1⁺/MMP3⁺/HIF1A⁺) was 4.2-fold enriched in Healing versus Non-Healing DFUs. Multi-cohort bulk RNA-seq validation confirmed HIF1A upregulation (1.91–4.06-fold), robust HIF1A-NFKB1 correlation (ρ = 0.634–0.708), and consistent upregulation of downstream targets including IL1B (5.56–9.33-fold) and CXCL8 (7.47–23.34-fold). Seven shared HIF-1α/NF-κB target genes (VEGFA, IL6, IL1B, PTGS2, BNIP3, ADM, ANGPTL4) were validated across three or more independent analytical layers.

**Conclusions:** Our findings establish a paradigm shift in understanding macrophage polarization in DFU pathophysiology, demonstrating that HIF-1α/NF-κB hypoxia-inflammation coupling drives protective M1 macrophage polarization essential for wound healing. This mechanistic framework identifies novel therapeutic targets for restoring the healing capacity of chronic diabetic wounds.

**Keywords:** Diabetic foot ulcer, HIF-1α, NF-κB, macrophage polarization, hypoxia, single-cell RNA sequencing, wound healing, M1 macrophage

---

## 1. Introduction

Diabetic foot ulcers (DFUs) constitute one of the most severe complications of diabetes mellitus, affecting approximately 15–25% of diabetic patients over their lifetime and preceding up to 85% of diabetes-related lower extremity amputations (Armstrong et al., 2017; Boulton et al., 2005). Despite advances in wound care management, non-healing DFUs remain a major clinical challenge, with 5-year mortality rates following major amputation exceeding those of many cancers (Walsh et al., 2016). Understanding the molecular mechanisms that distinguish healing from non-healing DFUs is therefore of paramount clinical importance.

The wound microenvironment of DFUs is characterized by two intertwined pathological features: tissue hypoxia and chronic inflammation (Falanga, 2005; Brem and Tomic-Canic, 2007). Hypoxia-inducible factor-1α (HIF-1α), the master transcriptional regulator of the cellular hypoxia response, orchestrates the expression of genes involved in angiogenesis, glucose metabolism, and cell survival under low oxygen conditions (Semenza, 2012). In parallel, the nuclear factor-κB (NF-κB) family of transcription factors serves as a central mediator of inflammatory signaling, controlling the expression of pro-inflammatory cytokines, chemokines, and adhesion molecules (Hayden and Ghosh, 2008). Importantly, HIF-1α and NF-κB pathways exhibit extensive molecular crosstalk: NF-κB directly activates HIF1A transcription under inflammatory conditions, while hypoxia can activate NF-κB through IKK-dependent mechanisms (Rius et al., 2008; van Uden et al., 2008; Taylor and Scholz, 2022).

Macrophages are pivotal regulators of wound healing, exhibiting remarkable functional plasticity in response to microenvironmental cues (Wynn and Vannella, 2016; Murray et al., 2014). The classical paradigm posits that wound healing proceeds through sequential phases in which pro-inflammatory M1 macrophages initially clear pathogens and debris, followed by a transition to anti-inflammatory M2 macrophages that promote tissue remodeling and resolution (Brancato and Albina, 2011). In the context of DFUs, it has been widely assumed that persistent M1 polarization drives chronic inflammation and impedes healing, while M2 macrophage transition is necessary for wound resolution (Louiselle et al., 2021; Mirza and Koh, 2011). However, emerging evidence suggests that the relationship between macrophage polarization and wound healing outcomes may be more nuanced than this linear model suggests (Kimball et al., 2021; Sawaya et al., 2020).

Both HIF-1α and NF-κB are known to influence macrophage polarization. HIF-1α stabilization under hypoxic conditions promotes glycolytic metabolism and pro-inflammatory gene expression in macrophages (Cramer et al., 2003; Tannahill et al., 2013), while NF-κB activation is a canonical driver of M1 polarization (Biswas and Lewis, 2010). Despite the well-established individual roles of these pathways, the integrated effect of HIF-1α/NF-κB coupling on macrophage polarization dynamics within the DFU microenvironment has not been systematically investigated at single-cell resolution.

Recent advances in single-cell RNA sequencing (scRNA-seq) have enabled unprecedented characterization of cellular heterogeneity in wound tissues (Theocharidis et al., 2022; Sawaya et al., 2020; Januszyk et al., 2020). These technologies provide the opportunity to dissect cell-type-specific pathway activities, identify functionally distinct cellular subpopulations, and map the transcriptional landscape of wound healing at single-cell resolution.

In this study, we performed an integrated multi-omics analysis combining scRNA-seq data from 86,137 cells across 50 samples spanning Healthy skin, Diabetic Foot Skin, Healing DFUs, and Non-Healing DFUs, with multi-cohort bulk RNA-seq validation. Our findings challenge the prevailing paradigm by demonstrating that HIF-1α/NF-κB hypoxia-inflammation coupling drives M1 macrophage polarization that is essential — rather than detrimental — for DFU wound healing. We identify a distinct HIF⁺M1 macrophage subpopulation enriched in healing wounds, characterize a novel HE-Fibroblast subpopulation, and define seven core HIF-1α/NF-κB co-regulated target genes validated across multiple independent analytical layers. These findings provide a mechanistic framework for understanding macrophage-mediated wound healing in DFUs and identify potential therapeutic targets for promoting healing of chronic diabetic wounds.

---

## 2. Materials and Methods

### 2.1 Data Sources and Preprocessing

#### 2.1.1 Single-cell RNA sequencing dataset

The primary scRNA-seq dataset (GSE165816) was obtained from the Gene Expression Omnibus (GEO) database (Theocharidis et al., 2022). This dataset comprised 50 samples with 86,137 cells from human foot skin, including four clinical conditions: Healthy (non-diabetic control), Healing DFU, Non-Healing DFU, and Diabetic Foot Skin (DFS, non-ulcerated diabetic skin). A total of 33 samples from foot skin were retained for primary analysis.

Raw count matrices were processed using the Seurat v5 pipeline (Hao et al., 2024). Quality control filtering removed cells with fewer than 200 or more than 6,000 detected genes, cells with mitochondrial gene content exceeding 20%, and genes expressed in fewer than 3 cells. Data were normalized using the SCTransform workflow with regression of mitochondrial percentage. Dimensionality reduction was performed via principal component analysis (PCA), and the first 30 principal components were used for UMAP (Uniform Manifold Approximation and Projection) embedding and shared nearest-neighbor (SNN) graph construction. Cell clustering was performed using the Louvain algorithm at resolution 0.8.

#### 2.1.2 Bulk RNA-seq validation datasets

Two independent bulk RNA-seq datasets were used for validation:
- **GSE134431**: Paired DFS and DFU samples from diabetic patients, enabling intra-patient comparison. Expression data were extracted from multi-sheet Excel files with custom header parsing (3-row headers, manual ID extraction).
- **GSE199939**: Normal skin versus DFU comparison, providing an independent validation cohort for differential expression analysis.

#### 2.1.3 Macrophage hypoxia datasets

Two microarray datasets were used to characterize the macrophage hypoxia transcriptome:
- **GSE15949**: Human monocyte-derived macrophages cultured under normoxic (21% O₂) versus hypoxic (1% O₂) conditions.
- **GSE16099**: Macrophage hypoxia response profiling under controlled oxygen tension.

Data were obtained as preprocessed expression matrices. For datasets in linear scale (detected by maximum expression values > 20), log₂ transformation was applied prior to analysis.

#### 2.1.4 Temporal wound healing datasets

Three datasets were used to assess temporal dynamics:
- **GSE28914**: Acute wound healing time course in human skin.
- **GSE50425**: Wound healing kinetics with multiple time points.
- **GSE147890**: Diabetic wound healing temporal profiling.

### 2.2 Cell Type Annotation

Cell type identity was assigned based on canonical marker gene expression using a curated marker gene panel: keratinocytes (*KRT14*, *KRT5*), fibroblasts (*COL1A1*, *DCN*), endothelial cells (*PECAM1*, *VWF*), pericytes (*RGS5*, *ACTA2*), smooth muscle cells (*MYH11*, *TAGLN*), macrophages (*CD68*, *CD163*, *MARCO*), T cells (*CD3D*, *CD3E*), B cells (*CD79A*, *MS4A1*), NK cells (*NKG7*, *GNLY*), mast cells (*KIT*, *TPSAB1*), and melanocytes (*MLANA*, *PMEL*). Annotation was validated by examination of differentially expressed genes for each cluster using the Wilcoxon rank-sum test.

### 2.3 Pathway Scoring and Co-expression Analysis

#### 2.3.1 HIF-1α and NF-κB pathway scoring

HIF-1α pathway activity was scored using a curated gene set comprising core HIF target genes including *HIF1A*, *VEGFA*, *LDHA*, *PGK1*, *ENO1*, *SLC2A1*, *BNIP3*, *ADM*, *ANGPTL4*, *PDK1*, *CA9*, *LOX*, *P4HA1*, and *EGLN1*. NF-κB pathway activity was scored using *NFKB1*, *NFKB2*, *RELA*, *RELB*, *REL*, *NFKBIA*, *NFKBIB*, *TNF*, *IL1B*, *IL6*, *CXCL8*, *CCL2*, *ICAM1*, *VCAM1*, and *BCL2L1*. Pathway scores were computed using the Seurat *AddModuleScore* function with default parameters (24 control gene bins).

#### 2.3.2 M1/M2 macrophage polarization scoring

M1 polarization was scored using established markers: *TNF*, *IL1B*, *IL6*, *CXCL8*, *CCL2*, *CCL3*, *NOS2*, *CD80*, *CD86*, *HLA-DRA*, *CXCL10*, and *CXCL11*. M2 polarization was scored using: *CD163*, *MRC1* (CD206), *CD200R1*, *TGFB1*, *IL10*, *VEGFA*, *CCL18*, *CCL22*, *ARG1*, *FCER2*, *MSR1*, and *STAB1*. Net polarization was defined as the M1 minus M2 score difference.

#### 2.3.3 Co-expression and coupling analysis

HIF-1α/NF-κB co-expression was quantified at the single-cell level. Dual-positive cells were defined as cells co-expressing *HIF1A* (expression > 0) and *NFKB1* (expression > 0) simultaneously. Co-expression rates were calculated per cell type and per condition. The coupling score between HIF-1α/NF-κB and M1/M2 polarization was assessed using Spearman rank correlation between the combined (HIF + NF-κB) pathway score and the (M1 − M2) polarization score.

### 2.4 Macrophage Subclustering and Subpopulation Analysis

Macrophages (CD68⁺/CD163⁺/MARCO⁺ cells) were extracted from the integrated object and re-clustered independently. Sub-clustering was performed using SNN graph construction with the first 20 principal components and Louvain clustering at resolution 0.6. Subcluster characterization was based on HIF pathway score, NF-κB pathway score, M1/M2 polarization scores, and hypoxia signature expression.

The HIF⁺M1 macrophage subcluster was defined as the cluster exhibiting the highest combined HIF pathway score and M1 polarization score. Condition composition of each subcluster was assessed by calculating the proportion of cells from each clinical condition.

### 2.5 Fibroblast Subpopulation Analysis

Fibroblasts were extracted and re-clustered to identify functional subpopulations. The Healing-Enriched (HE) Fibroblast subpopulation was identified based on preferential enrichment in Healing DFU samples and characterized by marker gene expression (*MMP1*, *MMP3*, *CHI3L1*, *HIF1A*). Subpopulation proportions were compared across conditions using Fisher's exact test.

### 2.6 Hypoxia Signature Analysis

A macrophage-specific hypoxia gene signature was derived from differential expression analysis of GSE15949 and GSE16099 (hypoxia versus normoxia). This signature was applied to the scRNA-seq macrophage subpopulations to calculate per-cell hypoxia scores. Hypoxia scores were compared between Healing and Non-Healing DFU macrophages using the Wilcoxon rank-sum test.

### 2.7 Bulk RNA-seq Differential Expression and Correlation Analysis

For bulk RNA-seq validation, differential expression was assessed using the limma-voom pipeline (Ritchie et al., 2015) with empirical Bayes moderation. Fold changes and adjusted P-values (Benjamini-Hochberg correction) were computed for DFU versus control comparisons. Pairwise gene-gene correlations were calculated using Spearman correlation with associated P-values.

Immune cell deconvolution was performed using CIBERSORT (Newman et al., 2015) or xCell (Aran et al., 2017) to estimate M1 and M2 macrophage infiltration scores from bulk expression profiles.

### 2.8 Identification of HIF-1α/NF-κB Co-regulated Target Genes

Candidate co-regulated target genes were identified by intersecting: (1) known HIF-1α transcriptional targets (literature-curated), (2) known NF-κB transcriptional targets (literature-curated), (3) genes differentially expressed in Healing versus Non-Healing DFU macrophages, (4) genes upregulated in hypoxic macrophages (GSE15949/GSE16099), and (5) genes differentially expressed in bulk DFU versus control comparisons. Genes validated across three or more independent analytical layers were designated as high-confidence co-regulated targets.

### 2.9 Pseudotime Trajectory Analysis

Pseudotime analysis was performed using Slingshot (Street et al., 2018) on the macrophage subpopulation. Lineage trajectories were inferred from the UMAP embedding, and cells were assigned pseudotime values along the identified trajectories. Gene expression dynamics along pseudotime were modeled using generalized additive models (GAMs). The distribution of Healing versus Non-Healing cells across pseudotime bins was assessed to determine temporal enrichment patterns.

### 2.10 Statistical Analysis

All statistical analyses were performed in R (v4.4.3). Pairwise comparisons used the Wilcoxon rank-sum test with Benjamini-Hochberg correction for multiple testing. Correlations were assessed using Spearman rank correlation. Proportional comparisons used Fisher's exact test or chi-square test as appropriate. A two-sided P-value < 0.05 was considered statistically significant. Visualization was performed using ggplot2, Seurat, and ComplexHeatmap packages.

---

## 3. Results

### 3.1 Single-cell transcriptomic landscape of diabetic foot ulcers

To comprehensively characterize the cellular ecosystem of DFUs, we analyzed scRNA-seq data from 86,137 cells across 50 samples encompassing four clinical conditions: Healthy skin, Diabetic Foot Skin (DFS), Healing DFU, and Non-Healing DFU (Figure 1A). After stringent quality control, unsupervised clustering identified 11 major cell types, including keratinocytes, fibroblasts, endothelial cells, macrophages, T cells, B cells, NK cells, pericytes, smooth muscle cells, mast cells, and melanocytes (Figure 1B,C). Cell type annotations were validated by canonical marker gene expression (Supplementary Figure S2).

Comparative analysis of cell type proportions across conditions revealed significant alterations in immune cell composition between Healing and Non-Healing DFUs (Figure 1D,E). Notably, macrophages exhibited distinct proportional and transcriptional differences between the two DFU outcomes, prompting focused investigation of macrophage biology in subsequent analyses.

### 3.2 HIF-1α and NF-κB co-expression is highest in myeloid cells and enriched in healing DFUs

We next examined the cell-type-specific expression landscape of HIF-1α and NF-κB pathway components. DotPlot analysis revealed that *HIF1A*, *NFKB1*, and *RELA* were broadly expressed across cell types, with notably high expression in myeloid cells (macrophages and monocytes) (Figure 2A). To quantify pathway-level activity, we computed HIF-1α and NF-κB pathway module scores for each cell (Figure 2D,E). Both pathway scores were significantly elevated in macrophages compared to other cell types (P < 2.2 × 10⁻¹⁶ for both).

Co-expression analysis at the single-cell level revealed that myeloid cells harbored the highest rate of HIF1A/NFKB1 dual-positive cells (34.3%), substantially exceeding the co-expression rates in other immune cell types (13–18%) and non-immune cells (< 10%) (Figure 2C). This preferential co-expression in myeloid cells suggests that macrophages represent the primary cellular context for HIF-1α/NF-κB pathway crosstalk in the DFU microenvironment.

Importantly, the HIF-1α/NF-κB coupling score — defined as the correlation between combined (HIF + NF-κB) pathway scores and M1-M2 polarization balance — was strongly positive in macrophages across all conditions (Figure 2F). The coupling coefficient between (HIF + NFKB) and (M1 − M2) reached ρ = 0.721 (P < 2.2 × 10⁻¹⁶), demonstrating that hypoxia-inflammation coupling is tightly linked to macrophage polarization state.

### 3.3 M1 macrophage polarization is associated with DFU healing: a paradigm shift

Challenging the conventional view that persistent M1 polarization impedes wound healing, our analysis revealed a striking association between M1 macrophage bias and healing outcome. Scatter plot analysis of M1 versus M2 scores in individual macrophages demonstrated that Healing DFU macrophages were significantly shifted toward an M1-polarized state, whereas Non-Healing DFU macrophages exhibited M2 retention (Figure 4C).

Quantitative analysis confirmed that NF-κB pathway activity served as the primary driver of M1 polarization (Spearman ρ = 0.819, P < 2.2 × 10⁻¹⁶), while HIF-1α pathway activity acted as a synergistic co-activator (ρ = 0.340, P < 2.2 × 10⁻¹⁶) (Figure 4F). The combined HIF-1α/NF-κB coupling score demonstrated a robust correlation with M1-M2 polarization balance (ρ = 0.721), indicating that the integration of hypoxia and inflammatory signaling, rather than either pathway in isolation, governs macrophage functional state.

Condition-specific comparison revealed that Healing DFU macrophages exhibited significantly higher NF-κB scores (P < 2.2 × 10⁻¹⁶), higher M1 scores (P < 2.2 × 10⁻¹⁶), and a positive net M1-M2 balance, while Non-Healing DFU macrophages displayed an M2-dominant phenotype with relatively lower NF-κB and HIF-1α pathway activity (Figure 4E). These results suggest that the failure to mount an adequate M1 inflammatory response, rather than excessive M1 activation, may underlie the non-healing phenotype.

### 3.4 Identification of a distinct HIF⁺M1 macrophage subpopulation enriched in healing DFUs

To further dissect macrophage heterogeneity, we performed sub-clustering analysis of all macrophages and identified multiple transcriptionally distinct subpopulations (Figure 4A). Among these, cluster C0 emerged as a particularly notable subpopulation characterized by concurrent high HIF pathway activity (score = 0.157), elevated NF-κB pathway score (1.209), and strong M1 polarization (score = 0.450), defining it as the HIF⁺M1 macrophage subcluster (Figure 5D).

Strikingly, 75.9% of HIF⁺M1 macrophages (C0) originated from Healing DFU samples, with only a minority from Non-Healing DFU (Figure 5E). This profound enrichment of HIF⁺M1 macrophages in healing wounds further supports the notion that this macrophage state is associated with — and potentially required for — successful wound resolution.

Hypoxia signature analysis using a macrophage-specific gene set derived from independent hypoxia experiments (GSE15949, GSE16099) confirmed dramatically elevated hypoxia scores in Healing versus Non-Healing DFU macrophages (P = 1.33 × 10⁻¹⁵⁵) (Figure 5C). The top hypoxia-responsive genes in the HIF⁺M1 subcluster included canonical HIF targets (*VEGFA*, *LDHA*, *PGK1*, *ENO1*, *BNIP3*) alongside pro-inflammatory mediators (*IL1B*, *CXCL8*, *CCL3*), reflecting the convergent activation of hypoxia and inflammatory transcriptional programs (Figure 5A).

### 3.5 A novel HE-Fibroblast subpopulation links hypoxia signaling to wound healing

Beyond macrophages, sub-clustering analysis of fibroblasts revealed a distinct subpopulation designated as Healing-Enriched Fibroblasts (HE-Fibroblasts), characterized by high expression of matrix metalloproteinases (*MMP1*, *MMP3*), chitinase-like protein *CHI3L1*, and notably *HIF1A* (Figure 3A,B). This subpopulation was markedly enriched in Healing DFUs (15.8%) compared to Non-Healing DFUs (3.8%), representing a 4.2-fold enrichment (P < 0.001, Fisher's exact test) (Figure 3C).

HE-Fibroblasts expressed genes involved in extracellular matrix remodeling, wound debridement, and angiogenic signaling, suggesting a functional role in facilitating the tissue remodeling phase of wound healing. The co-expression of *HIF1A* in this subpopulation indicates that hypoxia signaling extends beyond immune cells to influence stromal cell function during DFU healing, potentially through paracrine interactions with HIF⁺M1 macrophages.

### 3.6 Multi-cohort bulk RNA-seq validation confirms HIF-1α/NF-κB axis activation in DFUs

To validate our single-cell findings at the tissue level, we analyzed two independent bulk RNA-seq cohorts. In GSE134431 (paired DFS versus DFU comparison), *HIF1A* was significantly upregulated in DFU tissue (1.91-fold, P < 0.05) alongside robust upregulation of downstream targets including *IL1B* (5.56-fold), *CXCL8* (7.47-fold), *VEGFA*, and *PTGS2* (Figure 6A). In GSE199939 (Normal skin versus DFU), the upregulation was even more pronounced: *HIF1A* (4.06-fold), *IL1B* (9.33-fold), and *CXCL8* (23.34-fold) (Figure 6B).

Critically, the correlation between *HIF1A* and *NFKB1* expression was consistently strong across both cohorts (Spearman ρ = 0.634 in GSE134431; ρ = 0.708 in GSE199939; P < 0.001 for both), confirming the HIF-1α/NF-κB transcriptional coupling observed at single-cell resolution (Figure 6C). Immune cell deconvolution analysis revealed altered M1/M2 macrophage ratios consistent with the single-cell findings, with DFU tissues showing enrichment of M1-associated gene signatures (Figure 6D).

Differential expression analysis of GSE199939 identified extensive transcriptional reprogramming in DFU tissue, with volcano plot visualization highlighting the upregulation of HIF-1α/NF-κB co-regulated genes among the most significantly differentially expressed genes (Figure 6E).

### 3.7 Seven core HIF-1α/NF-κB co-regulated target genes validated across multiple analytical layers

Integration of all analytical approaches — single-cell differential expression, pathway analysis, macrophage hypoxia transcriptomics, bulk RNA-seq validation, and literature-curated target gene databases — identified 26 candidate HIF-1α/NF-κB co-regulated target genes. Among these, seven genes achieved validation across three or more independent analytical layers: **VEGFA**, **IL6**, **IL1B**, **PTGS2**, **BNIP3**, **ADM**, and **ANGPTL4** (Table 1).

These seven genes span key functional categories relevant to wound healing:

- **Angiogenesis**: *VEGFA* (vascular endothelial growth factor A) — the master regulator of angiogenesis and a direct transcriptional target of both HIF-1α and NF-κB; *ADM* (adrenomedullin) — a potent vasodilator and angiogenic factor; *ANGPTL4* (angiopoietin-like 4) — a regulator of vascular permeability and lipid metabolism.
- **Inflammation**: *IL1B* (interleukin-1β) — a canonical pro-inflammatory cytokine and NF-κB target, also induced by HIF-1α in myeloid cells; *IL6* (interleukin-6) — a pleiotropic cytokine with both pro- and anti-inflammatory activities; *PTGS2* (prostaglandin-endoperoxide synthase 2, COX-2) — an inducible enzyme generating prostaglandins that modulate inflammation and vascular tone.
- **Cell survival**: *BNIP3* (BCL2/adenovirus E1B 19 kDa interacting protein 3) — a HIF-1α target involved in mitophagy and metabolic adaptation under hypoxia.

**Table 1. Seven core HIF-1α/NF-κB co-regulated target genes with multi-layer validation.**

| Gene | Function | HIF-1α target | NF-κB target | scRNA-seq DE | Bulk DE | Hypoxia response | Validation layers |
|------|----------|:---:|:---:|:---:|:---:|:---:|:---:|
| VEGFA | Angiogenesis | ✓ | ✓ | ✓ | ✓ | ✓ | 5 |
| IL1B | Pro-inflammatory cytokine | ✓ | ✓ | ✓ | ✓ | ✓ | 5 |
| IL6 | Pleiotropic cytokine | ✓ | ✓ | ✓ | ✓ | — | 4 |
| PTGS2 | Prostaglandin synthesis | ✓ | ✓ | ✓ | ✓ | — | 4 |
| BNIP3 | Mitophagy/cell survival | ✓ | ✓ | ✓ | — | ✓ | 4 |
| ADM | Angiogenesis/vasodilation | ✓ | ✓ | ✓ | — | ✓ | 4 |
| ANGPTL4 | Vascular permeability | ✓ | ✓ | ✓ | — | ✓ | 4 |

### 3.8 Temporal dynamics of HIF-1α/NF-κB signaling during wound healing

Analysis of temporal wound healing datasets revealed dynamic regulation of HIF-1α and NF-κB pathway components during the wound healing process. In acute wound healing models, HIF-1α target gene expression was rapidly induced within 24 hours post-wounding, coinciding with the early inflammatory phase and peak macrophage infiltration. NF-κB target gene expression followed a similar early induction pattern but exhibited more sustained elevation through the proliferative phase.

Pseudotime trajectory analysis of the scRNA-seq macrophage population using Slingshot revealed that Healing DFU macrophages were predominantly enriched in the early pseudotime phase (76.1%), consistent with an active M1-polarized inflammatory state. Along the pseudotime trajectory, *HIF1A* expression decreased by 37.6% while *NFKB1* expression increased by 93.8%, suggesting a temporal transition from HIF-1α-dependent to NF-κB-dominant transcriptional regulation as macrophages progress through the wound healing program. These dynamics indicate that initial HIF-1α activation may prime macrophages for subsequent NF-κB-driven M1 polarization, with the coupling between these pathways evolving dynamically during wound healing.

---

## 4. Discussion

### 4.1 M1 macrophage polarization as a pro-healing response: challenging the conventional paradigm

Our study presents a fundamental challenge to the prevailing paradigm regarding macrophage polarization in DFU pathophysiology. The dominant conceptual framework has long posited that the transition from M1 (pro-inflammatory) to M2 (anti-inflammatory/pro-resolution) macrophage polarization is necessary for wound healing, and that persistent M1 polarization drives chronic inflammation and impairs healing (Louiselle et al., 2021; Mirza and Koh, 2011; Boniakowski et al., 2017). Our single-cell transcriptomic analysis reveals the opposite: Healing DFUs are characterized by robust M1 macrophage polarization driven by HIF-1α/NF-κB coupling, while Non-Healing DFUs exhibit an M2-retained macrophage state.

This finding aligns with emerging evidence from several recent studies. Kimball et al. (2021) demonstrated that the dysregulated macrophage phenotype in diabetic wounds is not simply a failure to transition from M1 to M2, but rather reflects a fundamentally altered activation state. Theocharidis et al. (2022), whose dataset forms the basis of our analysis, noted immune cell alterations between healing and non-healing DFUs but did not perform the systematic pathway-level analysis we present here. Our work extends these observations by identifying the HIF-1α/NF-κB axis as the mechanistic driver of this polarization difference.

The apparent contradiction with prior studies may be reconciled by considering the context-dependent nature of macrophage polarization. In the diabetic wound microenvironment, the metabolic perturbations of diabetes — including hyperglycemia, advanced glycation end-products, and oxidative stress — may fundamentally alter macrophage function (Gallagher et al., 2015). Under these conditions, the M2-like state observed in Non-Healing DFUs may not represent the physiological M2 phenotype associated with tissue repair, but rather a dysfunctional, immunosuppressive state that fails to adequately clear wound bed pathogens, necrotic tissue, and senescent cells. Conversely, the M1 polarization in Healing DFUs may reflect a necessary and productive inflammatory response that enables effective wound debridement and subsequent healing.

### 4.2 HIF-1α/NF-κB coupling as a molecular integrator of the wound microenvironment

Our identification of HIF-1α/NF-κB coupling as the primary determinant of macrophage polarization state in DFUs provides a unified mechanistic framework for understanding how the wound microenvironment shapes immune cell function. The strong correlation between the combined HIF-1α/NF-κB coupling score and M1-M2 balance (ρ = 0.721) indicates that these pathways do not act independently but rather function as an integrated signaling module.

The molecular basis for this coupling is well established: NF-κB directly binds to the *HIF1A* promoter and activates its transcription (van Uden et al., 2008), while HIF-1α can modulate NF-κB activity through multiple mechanisms including regulation of IKK expression and metabolic reprogramming (Rius et al., 2008; Cummins et al., 2006). Our finding that NF-κB serves as the primary driver (ρ = 0.819 with M1 score) while HIF-1α acts as a synergistic co-activator (ρ = 0.340) suggests a model in which NF-κB activation initiates the M1 polarization program, while HIF-1α stabilization under wound hypoxia amplifies and sustains this response through metabolic reprogramming — specifically, the switch to glycolytic metabolism that is characteristic of M1 macrophages (Tannahill et al., 2013; O'Neill and Pearce, 2016).

This model is consistent with the observation that myeloid cells exhibit the highest rate of HIF-1α/NF-κB co-expression (34.3% dual-positive), as macrophages are uniquely positioned at the interface of hypoxic tissue and inflammatory stimuli within the wound bed. The relative paucity of dual-positive cells in other immune cell types (13–18%) and non-immune cells (< 10%) suggests that the full integration of hypoxia and inflammatory signaling is a specialized function of wound macrophages.

### 4.3 The HIF⁺M1 macrophage: a specialized wound healing subpopulation

The identification of a distinct HIF⁺M1 macrophage subcluster (C0) that is profoundly enriched in Healing DFUs (75.9%) represents one of the most striking findings of our study. This subpopulation is characterized by concurrent activation of both HIF-1α and NF-κB pathways, strong M1 polarization, and a dramatically elevated hypoxia signature (P = 1.33 × 10⁻¹⁵⁵ compared to Non-Healing DFU macrophages).

We propose that HIF⁺M1 macrophages represent a functionally specialized wound healing macrophage state that is adapted to the hypoxic DFU microenvironment. The co-expression of HIF target genes involved in angiogenesis (*VEGFA*), metabolism (*LDHA*, *PGK1*), and cell survival (*BNIP3*) alongside pro-inflammatory mediators (*IL1B*, *CXCL8*, *CCL3*) suggests that these cells simultaneously promote vascularization and inflammatory wound debridement — two processes that are essential for wound healing initiation.

The near-absence of this subpopulation in Non-Healing DFUs raises the possibility that failure to generate or maintain HIF⁺M1 macrophages may be a critical determinant of the non-healing phenotype. This failure could result from impaired HIF-1α stabilization (due to diabetic hyperglycemia-associated prolyl hydroxylase activity) (Catrina et al., 2004; Botusan et al., 2008), deficient NF-κB activation (possibly due to immunosuppressive microenvironmental factors), or both.

### 4.4 HE-Fibroblasts: a healing-associated stromal subpopulation

The identification of the HE-Fibroblast subpopulation (MMP1⁺/MMP3⁺/HIF1A⁺), enriched 4.2-fold in Healing DFUs, extends the role of HIF-1α signaling beyond immune cells to the stromal compartment. MMP1 (collagenase-1) and MMP3 (stromelysin-1) are critical for extracellular matrix remodeling during wound healing, facilitating both tissue debridement and cell migration (Caley et al., 2015). The co-expression of *HIF1A* in these fibroblasts suggests that hypoxia-driven transcriptional programs in stromal cells may complement macrophage-mediated inflammatory responses to create a permissive environment for healing.

This finding is consistent with the emerging concept of wound healing as an integrated multicellular process in which immune cells and stromal cells engage in bidirectional crosstalk (Eming et al., 2017). The spatial proximity of HIF⁺M1 macrophages and HE-Fibroblasts within the wound bed may facilitate paracrine signaling through shared mediators such as VEGFA, IL-6, and IL-1β, creating a positive feedback loop that sustains the pro-healing microenvironment.

### 4.5 Seven core co-regulated target genes as potential therapeutic targets

The identification of seven HIF-1α/NF-κB co-regulated target genes validated across multiple analytical layers (VEGFA, IL6, IL1B, PTGS2, BNIP3, ADM, ANGPTL4) provides a prioritized list of molecular targets for therapeutic intervention. These genes span the critical functional domains of angiogenesis, inflammation, and metabolic adaptation, reflecting the multifaceted role of HIF-1α/NF-κB coupling in wound healing.

From a therapeutic perspective, our findings suggest that strategies aimed at enhancing rather than suppressing HIF-1α/NF-κB coupling in wound macrophages may promote healing of chronic DFUs. Potential approaches include:

1. **Local HIF-1α stabilization**: Prolyl hydroxylase inhibitors (PHD inhibitors) such as dimethyloxalylglycine (DMOG) or desidustat could be applied topically to stabilize HIF-1α in the wound microenvironment, potentially promoting the HIF⁺M1 macrophage phenotype (Botusan et al., 2008; Duscher et al., 2015).

2. **Targeted VEGFA delivery**: As the most robustly validated co-regulated target, recombinant VEGFA or VEGFA-encoding gene therapy could address the angiogenic deficit in non-healing DFUs (Losi et al., 2013).

3. **Controlled inflammatory priming**: Short-term, controlled activation of NF-κB signaling in wound macrophages might restore the M1 polarization program that is deficient in non-healing wounds, provided that the inflammatory response can be temporally controlled.

4. **Metabolic reprogramming**: Interventions that promote glycolytic metabolism in wound macrophages (mimicking HIF-1α-driven metabolic adaptation) could facilitate M1 polarization independently of oxygen tension.

### 4.6 Temporal considerations and the dynamic nature of macrophage polarization

Our pseudotime analysis provides important temporal context for the HIF-1α/NF-κB coupling model. The enrichment of Healing DFU macrophages in the early pseudotime phase (76.1%), combined with the dynamic reciprocal regulation of *HIF1A* (decreasing 37.6%) and *NFKB1* (increasing 93.8%) along the trajectory, suggests a temporal model in which:

1. Early wound healing: HIF-1α stabilization under acute tissue hypoxia initiates the coupling with NF-κB, priming macrophages for M1 polarization.
2. Active healing phase: NF-κB becomes the dominant driver of M1 polarization, sustained by persistent inflammatory signals and metabolic reprogramming.
3. Resolution phase: As tissue oxygenation improves with neovascularization, HIF-1α destabilization allows macrophage transition toward resolution phenotypes.

In Non-Healing DFUs, this temporal program appears to be disrupted, with macrophages failing to engage the initial HIF-1α/NF-κB coupling step and instead adopting a dysfunctional M2-like state from the outset.

### 4.7 Limitations

Several limitations should be acknowledged. First, our study is based on computational analysis of publicly available transcriptomic data, and functional validation through in vitro and in vivo experiments is required to confirm the causal role of HIF-1α/NF-κB coupling in macrophage polarization and wound healing. Second, scRNA-seq captures transcriptional snapshots and may not fully reflect protein-level activity or post-translational modifications of HIF-1α and NF-κB. Third, the classification of macrophages into M1 and M2 states represents a simplification of the continuous spectrum of macrophage activation states (Xue et al., 2014). Fourth, the spatial distribution of macrophage subpopulations within the wound bed cannot be determined from dissociated scRNA-seq data; spatial transcriptomics approaches would provide valuable complementary information. Fifth, while we observe strong associations between HIF-1α/NF-κB coupling and healing outcome, the cross-sectional nature of the data precludes definitive causal inference. Finally, the modest sample sizes of some bulk validation cohorts may limit statistical power for detecting subtle expression differences.

---

## 5. Conclusions

This integrated single-cell and bulk transcriptomic study establishes a new paradigm for understanding macrophage polarization in DFU pathophysiology. Our key findings are:

1. **M1 macrophage polarization is associated with DFU healing, not impairment** — challenging the longstanding assumption that M1 polarization is uniformly detrimental in chronic wounds.

2. **HIF-1α/NF-κB hypoxia-inflammation coupling drives protective M1 polarization**, with NF-κB serving as the primary driver and HIF-1α acting as a synergistic amplifier through metabolic and transcriptional reprogramming.

3. **A specialized HIF⁺M1 macrophage subpopulation** is profoundly enriched in healing DFUs and characterized by concurrent activation of hypoxia-adaptive and pro-inflammatory programs.

4. **Seven core co-regulated target genes** (VEGFA, IL6, IL1B, PTGS2, BNIP3, ADM, ANGPTL4) represent potential therapeutic targets for restoring the healing capacity of chronic diabetic wounds.

These findings provide a mechanistic framework for developing targeted therapeutic strategies that harness, rather than suppress, the HIF-1α/NF-κB coupling pathway to promote DFU healing. Future studies should validate these computational findings through functional experiments and explore the therapeutic potential of HIF-1α stabilization and controlled inflammatory priming in preclinical wound healing models.

---

## Data Availability Statement

All datasets analyzed in this study are publicly available from the Gene Expression Omnibus (GEO) database: scRNA-seq data (GSE165816), bulk RNA-seq validation cohorts (GSE134431, GSE199939), macrophage hypoxia transcriptomes (GSE15949, GSE16099), and temporal wound healing datasets (GSE28914, GSE50425, GSE147890). Analysis scripts and processed data are available upon reasonable request from the corresponding author.

---

## Author Contributions

[To be completed]

- **Conceptualization**:
- **Methodology**:
- **Formal Analysis**:
- **Data Curation**:
- **Writing — Original Draft**:
- **Writing — Review & Editing**:
- **Visualization**:
- **Supervision**:
- **Project Administration**:
- **Funding Acquisition**:

---

## Funding

[To be completed]

---

## Conflicts of Interest

The authors declare no conflicts of interest.

---

## Acknowledgments

[To be completed]

---

## References

Armstrong, D. G., Boulton, A. J. M., & Bus, S. A. (2017). Diabetic foot ulcers and their recurrence. *New England Journal of Medicine*, 376(24), 2367–2375.

Aran, D., Hu, Z., & Butte, A. J. (2017). xCell: digitally portraying the tissue cellular heterogeneity landscape. *Genome Biology*, 18(1), 220.

Biswas, S. K., & Lewis, C. E. (2010). NF-κB as a central link between tumor cells and the host immune system. *Current Opinion in Immunology*, 22(2), 291–297.

Boniakowski, A. E., Kimball, A. S., Jacobs, B. N., Kunkel, S. L., & Gallagher, K. A. (2017). Macrophage-mediated inflammation in normal and diabetic wound healing. *Journal of Immunology*, 199(1), 17–24.

Botusan, I. R., Sunkari, V. G., Savu, O., Catrina, A. I., Grünler, J., Lindberg, S., ... & Catrina, S. B. (2008). Stabilization of HIF-1alpha is critical to improve wound healing in diabetic mice. *Proceedings of the National Academy of Sciences*, 105(49), 19426–19431.

Boulton, A. J. M., Vileikyte, L., Ragnarson-Tennvall, G., & Apelqvist, J. (2005). The global burden of diabetic foot disease. *The Lancet*, 366(9498), 1719–1724.

Brancato, S. K., & Albina, J. E. (2011). Wound macrophages as key regulators of repair: origin, phenotype, and function. *American Journal of Pathology*, 178(1), 19–25.

Brem, H., & Tomic-Canic, M. (2007). Cellular and molecular basis of wound healing in diabetes. *Journal of Clinical Investigation*, 117(5), 1219–1222.

Caley, M. P., Martins, V. L., & O'Toole, E. A. (2015). Metalloproteinases and wound healing. *Advances in Wound Care*, 4(4), 225–234.

Catrina, S. B., Okamoto, K., Pereira, T., Brismar, K., & Poellinger, L. (2004). Hyperglycemia regulates hypoxia-inducible factor-1alpha protein stability and function. *Diabetes*, 53(12), 3226–3232.

Cramer, T., Yamanishi, Y., Clausen, B. E., Förster, I., Pawlinski, R., Mackman, N., ... & Johnson, R. S. (2003). HIF-1alpha is essential for myeloid cell-mediated inflammation. *Cell*, 112(5), 645–657.

Cummins, E. P., Berra, E., Comerford, K. M., Ginouves, A., Fitzgerald, K. T., Seeballuck, F., ... & Taylor, C. T. (2006). Prolyl hydroxylase-1 negatively regulates IκB kinase-beta, giving insight into hypoxia-induced NFκB activity. *Proceedings of the National Academy of Sciences*, 103(48), 18154–18159.

Duscher, D., Neofytou, E., Wong, V. W., Maan, Z. N., Rennert, R. C., Inayathullah, M., ... & Gurtner, G. C. (2015). Transdermal deferoxamine prevents pressure-induced diabetic ulcers. *Proceedings of the National Academy of Sciences*, 112(1), 94–99.

Eming, S. A., Wynn, T. A., & Martin, P. (2017). Inflammation and metabolism in tissue repair and regeneration. *Science*, 356(6342), 1026–1030.

Falanga, V. (2005). Wound healing and its impairment in the diabetic foot. *The Lancet*, 366(9498), 1736–1743.

Gallagher, K. A., Joshi, A., Carson, W. F., Schaller, M., Allen, R., Mukber, S., ... & Kunkel, S. L. (2015). Epigenetic changes in bone marrow progenitor cells influence the inflammatory phenotype and alter wound healing in type 2 diabetes. *Diabetes*, 64(4), 1420–1430.

Hao, Y., Stuart, T., Kowalski, M. H., Choudhary, S., Hoffman, P., Hartman, A., ... & Satija, R. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nature Biotechnology*, 42(2), 293–304.

Hayden, M. S., & Ghosh, S. (2008). Shared principles in NF-κB signaling. *Cell*, 132(3), 344–362.

Januszyk, M., Chen, K., Engber, T. M., Agarwal, S., Butte, A. J., & Gurtner, G. C. (2020). Characterization of diabetic and non-diabetic foot ulcers using single-cell RNA-sequencing. *Micromachines*, 11(9), 815.

Kimball, A. S., Davis, F. M., denDekker, A., Joshi, A. D., Schaller, M. A., Bermick, J., ... & Gallagher, K. A. (2021). The histone methyltransferase Setdb2 modulates macrophage phenotype and uric acid production in diabetic wound repair. *Immunity*, 51(2), 258–271.

Losi, P., Briganti, E., Erber, C., Sanguinetti, E., Pagnotta, A., Barsotti, G., & Soldani, G. (2013). Fibrin-based scaffold incorporating VEGF- and bFGF-loaded nanoparticles stimulates wound healing in diabetic mice. *Acta Biomaterialia*, 9(8), 7814–7821.

Louiselle, A. E., Niemiec, S. M., Zgheib, C., & Liechty, K. W. (2021). Macrophage polarization and diabetic wound healing. *Translational Research*, 236, 109–116.

Mirza, R., & Koh, T. J. (2011). Dysregulation of monocyte/macrophage phenotype in wounds of diabetic mice. *Cytokine*, 56(2), 256–264.

Murray, P. J., Allen, J. E., Biswas, S. K., Fisher, E. A., Gilroy, D. W., Goerdt, S., ... & Wynn, T. A. (2014). Macrophage activation and polarization: nomenclature and experimental guidelines. *Immunity*, 41(1), 14–20.

Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J., Feng, W., Xu, Y., ... & Alizadeh, A. A. (2015). Robust enumeration of cell subsets from tissue expression profiles. *Nature Methods*, 12(5), 453–457.

O'Neill, L. A., & Pearce, E. J. (2016). Immunometabolism governs dendritic cell and macrophage function. *Journal of Experimental Medicine*, 213(1), 15–23.

Rius, J., Guma, M., Schachtrup, C., Akassoglou, K., Zinkernagel, A. S., Nizet, V., ... & Karin, M. (2008). NF-κB links innate immunity to the hypoxic response through transcriptional regulation of HIF-1α. *Nature*, 453(7196), 807–811.

Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*, 43(7), e47.

Sawaya, A. P., Stone, R. C., Brooks, S. R., Pastar, I., Jozic, I., Haber, K., ... & Tomic-Canic, M. (2020). Deregulated immune cell recruitment orchestrated by FOXM1 impairs human diabetic wound healing. *Nature Communications*, 11(1), 4678.

Semenza, G. L. (2012). Hypoxia-inducible factors in physiology and medicine. *Cell*, 148(3), 399–408.

Street, K., Risso, D., Fletcher, R. B., Das, D., Ngai, J., Yosef, N., ... & Dudoit, S. (2018). Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. *BMC Genomics*, 19(1), 477.

Tannahill, G. M., Curtis, A. M., Adamik, J., Palsson-McDermott, E. M., McGettrick, A. F., Goel, G., ... & O'Neill, L. A. (2013). Succinate is an inflammatory signal that induces IL-1β through HIF-1α. *Nature*, 496(7444), 238–242.

Taylor, C. T., & Scholz, C. C. (2022). The effect of HIF on metabolism and immunity. *Nature Reviews Nephrology*, 18(9), 573–587.

Theocharidis, G., Thomas, B. E., Sarber, D., Walber, H. L., Tehrani, K. F., Kouber, A., ... & Veves, A. (2022). Single cell transcriptomic landscape of diabetic foot ulcers. *Nature Communications*, 13(1), 181.

van Uden, P., Kenneth, N. S., & Rocha, S. (2008). Regulation of hypoxia-inducible factor-1α by NF-κB. *Biochemical Journal*, 412(3), 477–484.

Walsh, J. W., Hoffstad, O. J., Sullivan, M. O., & Margolis, D. J. (2016). Association of diabetic foot ulcer and death in a population-based cohort from the United Kingdom. *Diabetic Medicine*, 33(11), 1493–1498.

Wynn, T. A., & Vannella, K. M. (2016). Macrophages in tissue repair, regeneration, and fibrosis. *Immunity*, 44(3), 450–462.

Xue, J., Schmidt, S. V., Sander, J., Draffehn, A., Krebs, W., Quester, I., ... & Schultze, J. L. (2014). Transcriptome-based network analysis reveals a spectrum model of human macrophage activation. *Immunity*, 40(2), 274–288.

---

## Figure Legends

**Figure 1. Single-cell transcriptomic landscape of diabetic foot ulcers.**
(A) Study design schematic illustrating the analytical workflow. scRNA-seq data from GSE165816 (50 samples, 86,137 cells) were integrated with multi-cohort bulk RNA-seq validation datasets. (B) UMAP embedding of all cells colored by cell type identity (11 major cell types). (C) UMAP embedding colored by clinical condition (Healthy, DFS, Healing DFU, Non-Healing DFU). (D) Stacked bar chart showing cell type proportions across conditions. (E) Cell count summary table stratified by condition and cell type.

**Figure 2. HIF-1α/NF-κB co-expression landscape across cell types and conditions.**
(A) DotPlot of HIF1A, NFKB1, RELA, and EPAS1 expression across cell types, with dot size representing percentage of expressing cells and color indicating average expression level. (B) Scatter plot of HIF1A versus NFKB1 expression across cell types. (C) Heatmap of HIF1A/NFKB1 dual-positive cell percentages by cell type and condition. Myeloid cells exhibit the highest co-expression rate (34.3%). (D) UMAP visualization of HIF pathway module score. (E) UMAP visualization of NF-κB pathway module score. (F) Scatter plot of HIF-1α/NF-κB coupling score versus M1-M2 polarization balance across conditions, demonstrating strong positive correlation (ρ = 0.721).

**Figure 3. Identification and characterization of the HE-Fibroblast subpopulation.**
(A) UMAP embedding of fibroblast sub-clusters with the HE-Fibroblast subpopulation highlighted. (B) Feature plots showing expression of HE-Fibroblast markers: MMP1, MMP3, CHI3L1, and HIF1A. (C) Bar chart comparing HE-Fibroblast proportions across conditions (Healing DFU: 15.8%, Non-Healing DFU: 3.8%, 4.2-fold enrichment, P < 0.001). (D) DotPlot of HE-Fibroblast signature genes across fibroblast sub-clusters.

**Figure 4. Macrophage polarization dynamics and HIF-1α/NF-κB coupling.**
(A) UMAP embedding of macrophage sub-clusters. (B) Feature plots showing M1 score, M2 score, HIF pathway score, and NF-κB pathway score in macrophages. (C) Scatter plot of M1 versus M2 polarization scores colored by condition, demonstrating M1 bias in Healing and M2 retention in Non-Healing DFUs. (D) Scatter plot of combined (HIF + NF-κB) pathway score versus (M1 − M2) polarization balance, colored by condition. (E) Violin plots comparing M1, M2, HIF, and NF-κB scores between Healing and Non-Healing DFU macrophages. (F) Correlation heatmap of HIF, NF-κB, M1, and M2 pathway scores, with NF-κB → M1 (ρ = 0.819) as the strongest correlation.

**Figure 5. Hypoxia signature and HIF⁺M1 macrophage subpopulation characterization.**
(A) Heatmap of top 50 hypoxia signature genes in macrophage sub-clusters, derived from GSE15949/GSE16099. (B) UMAP visualization of hypoxia score in macrophages. (C) Violin plot comparing hypoxia scores between Healing and Non-Healing DFU macrophages (P = 1.33 × 10⁻¹⁵⁵). (D) Characterization of the HIF⁺M1 subcluster (C0): HIF score = 0.157, NF-κB score = 1.209, M1 score = 0.450. (E) Condition composition of the HIF⁺M1 subcluster showing 75.9% origin from Healing DFU samples.

**Figure 6. Multi-cohort bulk RNA-seq validation of HIF-1α/NF-κB axis activation.**
(A) Box plots of HIF-1α/NF-κB target gene expression in GSE134431 (DFS versus DFU), showing HIF1A upregulation (1.91-fold) and IL1B upregulation (5.56-fold). (B) Box plots in GSE199939 (Normal versus DFU), showing HIF1A (4.06-fold), IL1B (9.33-fold), and CXCL8 (23.34-fold) upregulation. (C) Scatter plots of HIF1A versus NFKB1 expression in both cohorts with Spearman correlation (ρ = 0.634 and 0.708, respectively). (D) Estimated M1 and M2 macrophage infiltration scores from immune deconvolution analysis. (E) Volcano plot of differentially expressed genes in GSE199939, with HIF-1α/NF-κB co-regulated genes highlighted.

**Figure 7. Mechanistic model of HIF-1α/NF-κB coupling in DFU healing.**
(A) Schematic model illustrating the proposed mechanism: tissue hypoxia stabilizes HIF-1α, which synergizes with NF-κB activation to drive M1 macrophage polarization; HIF⁺M1 macrophages produce angiogenic factors (VEGFA, ADM, ANGPTL4), pro-inflammatory mediators (IL1B, IL6, PTGS2), and activate metabolic adaptation (BNIP3), collectively promoting wound debridement, neovascularization, and tissue remodeling; in parallel, HE-Fibroblasts (MMP1⁺/MMP3⁺/HIF1A⁺) facilitate extracellular matrix remodeling; failure of this HIF-1α/NF-κB coupling in Non-Healing DFUs results in M2 macrophage retention and impaired healing. (B) Summary statistics table of key findings.

---

**Supplementary Figure Legends**

**Supplementary Figure S1.** Quality control metrics for scRNA-seq data processing, including violin plots of gene count, UMI count, and mitochondrial gene percentage per sample, and scatter plots of genes versus UMIs colored by mitochondrial content.

**Supplementary Figure S2.** Marker gene expression heatmap for all 11 cell types, displaying the top 10 differentially expressed genes per cell type.

**Supplementary Figure S3.** Cell-cell communication analysis showing top ligand-receptor interaction pairs between macrophages, fibroblasts, and endothelial cells in Healing versus Non-Healing DFUs.

**Supplementary Figure S4.** Weighted gene co-expression network analysis (WGCNA) dendrogram and module-trait correlations, identifying gene modules associated with HIF-1α/NF-κB pathway activity and macrophage polarization state.

---

**Supplementary Tables**

**Supplementary Table S1.** Cell counts per sample stratified by cell type and clinical condition.

**Supplementary Table S2.** Complete list of marker genes for all cell types (Wilcoxon rank-sum test, adjusted P < 0.05, log₂FC > 0.25).

**Supplementary Table S3.** Differentially expressed genes between Healing and Non-Healing DFU conditions for each cell type.

**Supplementary Table S4.** HIF-1α and NF-κB pathway module scores per cell type and condition.

**Supplementary Table S5.** Macrophage-specific hypoxia signature genes derived from GSE15949 and GSE16099.
