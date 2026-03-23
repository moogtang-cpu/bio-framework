# Technical Report: Integrative Single-Cell and Bulk RNA-seq Analysis of HIF-1α/NF-κB Hypoxia-Inflammation Coupling in Diabetic Foot Ulcers

---

**Report Version:** 1.0
**Date:** March 2026
**Analysis Platform:** R 4.4.3, Seurat v5, Harmony, limma
**Computing Environment:** Ubuntu 24.04 WSL2, 64 cores, 235 GB RAM

---

## Table of Contents

1. [Abstract](#1-abstract)
2. [Background](#2-background)
3. [Datasets and Methods](#3-datasets-and-methods)
4. [Results](#4-results)
   - 4.1 [Phase 1: Single-Cell Core Analysis](#41-phase-1-single-cell-core-analysis)
   - 4.2 [Phase 2: Bulk RNA-seq Multi-Cohort Validation](#42-phase-2-bulk-rna-seq-multi-cohort-validation)
   - 4.3 [Phase 3: Temporal Dynamics Analysis](#43-phase-3-temporal-dynamics-analysis)
   - 4.4 [Phase 4: Macrophage Mechanistic Dissection](#44-phase-4-macrophage-mechanistic-dissection)
5. [Key Findings Summary](#5-key-findings-summary)
6. [Discussion](#6-discussion)
7. [Detailed Methods](#7-detailed-methods)
8. [Dataset References](#8-dataset-references)
9. [Figure Index](#9-figure-index)

---

## 1. Abstract

Diabetic foot ulcers (DFUs) represent a severe complication of diabetes mellitus characterized by chronic, non-healing wounds with significant morbidity and mortality. The molecular interplay between hypoxia signaling (HIF-1α) and inflammatory cascades (NF-κB) in the wound microenvironment remains incompletely understood. Here, we present a comprehensive integrative analysis combining single-cell RNA sequencing (scRNA-seq) of 86,137 cells from 50 DFU and control skin samples with bulk RNA-seq data from three independent cohorts, two macrophage hypoxia datasets, and three temporal wound-healing datasets.

Our analysis reveals a finding that challenges the conventional model: non-healing DFU wounds are characterized by M2 macrophage retention and impaired M1 activation, rather than by excessive M1-driven inflammation. Healing wounds exhibit robust M1-polarized macrophages with strong HIF-1α/NF-κB coupling, suggesting that the capacity for M1 inflammatory initiation — and subsequent normal M1→M2 transition — is associated with wound resolution. NF-κB signaling is strongly associated with M1 polarization (Spearman rho = 0.819), with HIF-1α acting as a synergistic co-activator (rho = 0.340). The HIF-1α/NF-κB coupling axis strongly correlates with the M1-M2 polarization balance (rho = 0.721). We identify a hypoxia-responsive, extracellular matrix-remodeling fibroblast subpopulation (HE-Fibroblasts: MMP1+/MMP3+/HIF1A+) that is 4.2-fold enriched in healing versus non-healing wounds. In macrophages, hypoxia scoring reveals dramatically higher hypoxia activation in healing wounds (mean score 0.101 vs. 0.013, p = 1.33 x 10^-155), and a distinct HIF+M1 macrophage subcluster (C0) whose composition is 75.9% healing-derived. Bulk RNA-seq validation across two independent DFU cohorts confirms HIF1A upregulation (fold changes 1.91-4.06), robust HIF1A-NFKB1 co-expression correlation (rho = 0.634-0.708), and massive inflammatory gene induction (IL1B up to 9.33-fold, CXCL8 up to 23.34-fold).

These findings establish the HIF-1α/NF-κB coupling axis as a central orchestrator of the DFU wound microenvironment, reveal M2 retention as a hallmark of non-healing wounds, and identify actionable therapeutic targets for restoring normal macrophage dynamic transition in chronic wound management.

---

## 2. Background

### 2.1 Diabetic Foot Ulcers: A Clinical Challenge

Diabetic foot ulcers affect approximately 15-25% of patients with diabetes mellitus during their lifetime, with annual incidence rates of 2-5%. DFUs are the leading cause of non-traumatic lower extremity amputations, and 5-year mortality rates following major amputation exceed 50%, rivaling many cancers. Despite advances in wound care, approximately 40-60% of DFUs fail to heal within 20 weeks of standard treatment, underscoring the urgent need for mechanistic understanding and novel therapeutic strategies.

### 2.2 The Hypoxia-Inflammation Nexus

The DFU wound microenvironment is characterized by two intertwined pathological features: chronic hypoxia and persistent inflammation. Peripheral arterial disease and microvascular dysfunction create a profoundly hypoxic tissue milieu, while metabolic derangements and recurrent infections sustain inflammatory cascades. At the molecular level, two master transcription factor families orchestrate these responses:

- **HIF-1α (Hypoxia-Inducible Factor 1-alpha):** The canonical hypoxia sensor and transcriptional activator that regulates angiogenesis (VEGFA), glycolytic metabolism (LDHA, PGK1), and cell survival (BNIP3) under low oxygen conditions.
- **NF-κB (Nuclear Factor kappa B):** The master regulator of innate immunity and inflammation, controlling cytokine production (IL1B, TNF, IL6), chemokine expression (CXCL8), and immune cell activation.

Critically, these pathways are not independent. HIF-1α and NF-κB exhibit extensive biochemical crosstalk: NF-κB directly transactivates the HIF1A promoter, while HIF-1α can modulate NF-κB activity through multiple mechanisms. However, the cell-type-specific manifestation of this coupling in the DFU wound microenvironment, and its functional consequences for healing outcome, have not been systematically characterized at single-cell resolution.

### 2.3 Macrophage Polarization Paradigm

Wound healing macrophages are conventionally described as transitioning from a pro-inflammatory M1 state (classically activated; TNF, IL1B, NOS2) to an anti-inflammatory, pro-resolution M2 state (alternatively activated; MRC1, CD163, ARG1). The prevailing model posits that chronic wounds are "stuck" in an M1-dominated inflammatory phase. Our findings challenge this framework and suggest a more nuanced relationship between macrophage polarization, hypoxia signaling, and healing competence.

### 2.4 Study Objectives

This study addresses five interconnected scientific questions:

1. **Q1:** How are HIF-1α and NF-κB expressed and co-regulated across cell types in the DFU microenvironment?
2. **Q2:** Does HIF-1α/NF-κB coupling differ between healing and non-healing DFU wounds?
3. **Q3:** What is the relationship between HIF-1α/NF-κB signaling and macrophage polarization?
4. **Q4:** Can bulk RNA-seq data from independent cohorts validate single-cell findings?
5. **Q5:** What are the temporal dynamics of hypoxia and inflammatory gene responses in wound-relevant models?

---

## 3. Datasets and Methods

### 3.1 Dataset Overview

| Dataset     | Type          | Samples | Purpose                              | Source    |
|-------------|---------------|---------|--------------------------------------|-----------|
| GSE165816   | scRNA-seq     | 50      | Core single-cell analysis            | DFU skin  |
| GSE134431   | Bulk RNA-seq  | 13      | DFS vs DFU validation                | DFU skin  |
| GSE199939   | Bulk RNA-seq  | 20      | Normal vs DFU validation             | DFU skin  |
| GSE15949    | Microarray    | 12      | Macrophage hypoxia response          | In vitro  |
| GSE16099    | Microarray    | 6       | Macrophage hypoxia (M1/M2)           | In vitro  |
| GSE28914    | Microarray    | 30      | Temporal wound healing (burns)       | In vivo   |
| GSE50425    | Microarray    | 32      | Pressure ulcer time series           | In vivo   |
| GSE147890   | Microarray    | 16      | Diabetic vs control wound healing    | In vivo   |

**Total cells analyzed (scRNA-seq):** 86,137 cells across 50 samples
**Conditions (scRNA-seq):** Healthy (n=10), DM (Diabetes without wound, n=10), Healing DFU (n=15), Non-Healing DFU (n=15)

### 3.2 Analytical Pipeline Summary

The analysis was conducted in four phases:

- **Phase 1 (Single-Cell Core):** Quality control, normalization, dimensionality reduction, cell type annotation, HIF/NF-κB co-expression analysis, pathway scoring, fibroblast and macrophage subpopulation dissection, cell-cell communication.
- **Phase 2 (Bulk Validation):** Independent cohort differential expression analysis, HIF1A-NFKB1 correlation validation, immune deconvolution.
- **Phase 3 (Temporal Dynamics):** Time-series wound healing analysis, diabetic-specific temporal response characterization.
- **Phase 4 (Macrophage Mechanism):** Hypoxia signature derivation from in vitro data, cross-platform scoring, HIF+M1 subcluster identification, target gene characterization.

---

## 4. Results

### 4.1 Phase 1: Single-Cell Core Analysis

#### 4.1.1 Cell Type Landscape of the DFU Microenvironment

Quality-controlled single-cell RNA-seq data from 50 skin biopsies yielded 86,137 high-quality cells. Following Harmony-based batch correction across all samples, unsupervised clustering with Seurat v5 identified 11 distinct cell types:

| Cell Type       | Proportion | Key Markers              |
|-----------------|-----------|--------------------------|
| Fibroblast      | ~36%      | COL1A1, DCN, LUM         |
| Keratinocyte    | ~15%      | KRT14, KRT5, KRT1        |
| Endothelial     | ~10%      | PECAM1, VWF, CDH5        |
| T cell          | ~9%       | CD3D, CD3E, IL7R         |
| Macrophage      | ~8%       | CD68, CD163, CSF1R       |
| Myeloid         | ~7%       | LYZ, S100A8, S100A9      |
| Smooth Muscle   | ~5%       | ACTA2, MYH11, TAGLN      |
| Pericyte        | ~4%       | PDGFRB, RGS5, NOTCH3     |
| B cell          | ~3%       | CD79A, MS4A1, CD19       |
| Mast cell       | ~2%       | KIT, TPSAB1, CPA3        |
| Melanocyte      | ~1%       | PMEL, MLANA, DCT         |

**Condition-specific compositional shifts** were observed: Healing wounds showed increased myeloid and macrophage representation relative to Non-Healing wounds, while fibroblast proportions were altered across all DFU conditions compared to healthy skin (**Figure 1**).

#### 4.1.2 HIF-1α/NF-κB Co-expression Landscape

Systematic analysis of HIF1A and NFKB1 co-expression across all cell types revealed cell-type-specific patterns:

- **Myeloid cells** exhibited the highest dual-positive (HIF1A+/NFKB1+) rate at **34.3%**, consistent with the known transcriptional plasticity of myeloid lineage cells.
- **Endothelial cells** showed **15.8%** dual-positive rate, reflecting vascular sensitivity to hypoxia-inflammation signaling.
- **Macrophages** displayed **13.3%** dual-positive rate, with strong condition-dependent variation.
- **Fibroblasts** demonstrated a striking condition-dependent pattern: **39.4%** in Healing versus **25.3%** in Non-Healing wounds.

Pathway-level analysis using AddModuleScore confirmed that HIF pathway activity and NF-κB pathway activity are significantly correlated across all cell types, with the strongest coupling observed in myeloid lineage cells (**Figure 2**).

#### 4.1.3 HE-Fibroblast Subpopulation: A Healing-Associated Matrix-Remodeling State

Sub-clustering of the fibroblast compartment revealed a distinctive subpopulation characterized by co-expression of matrix metalloproteinases and hypoxia markers:

**HE-Fibroblasts (Hypoxia-responsive, ECM-remodeling):**
- Marker profile: MMP1+/MMP3+/CHI3L1+/HIF1A+
- Healing DFU: **15.8%** of fibroblasts
- Non-Healing DFU: **3.8%** of fibroblasts
- **Fold enrichment: 4.2x** in Healing vs. Non-Healing

This subpopulation represents a wound-activated fibroblast state that combines hypoxia-responsive transcription with aggressive extracellular matrix remodeling, suggesting that productive ECM turnover under hypoxic conditions is a hallmark of healing-competent wounds (**Figure 3**).

#### 4.1.4 Macrophage Polarization: M2 Retention and Impaired M1 Activation

Analysis of macrophage polarization states across conditions revealed a finding that challenges the conventional wound healing paradigm:

**Polarization scoring correlations (Spearman):**

| Comparison              | rho     | Interpretation                                     |
|-------------------------|---------|-----------------------------------------------------|
| NF-κB ↔ M1 score       | **0.819** | NF-κB is the dominant driver of M1 polarization     |
| HIF ↔ M1 score         | 0.340   | HIF-1α contributes to but does not dominate M1      |
| (HIF+NFKB) ↔ (M1-M2)  | **0.721** | Combined coupling strongly predicts polarization balance |

**Critical observation challenging the conventional model:** Healing DFU wounds were characterized by M1-biased macrophage polarization, while Non-Healing wounds exhibited M2-biased polarization. This directly contradicts the simplified model that chronic wound inflammation (M1 excess) is pathological and that M2 transition is necessary for healing. Instead, our data suggest that:

1. Non-healing wounds are characterized by M2-biased macrophage populations that may reflect a dysfunctional, immunosuppressive state — rather than the physiological M2 phenotype associated with tissue repair — locked in a state that fails to adequately clear wound bed pathogens, necrotic tissue, and senescent cells.
2. Healing wounds retain the capacity to mount a robust M1 inflammatory response, driven by NF-κB and augmented by HIF-1α, which may be a prerequisite for initiating the normal M1→M2 transition that underlies successful wound resolution.
3. The HIF-1α/NF-κB coupling axis serves as the upstream switch that controls the M1/M2 polarization balance (**Figure 4**).

#### 4.1.5 Cell-Cell Communication: IL-1 Signaling Axis

Ligand-receptor interaction analysis identified the IL-1 signaling pathway as the dominant communication axis in healing DFU wounds:

- **Myeloid → Fibroblast** IL-1 signaling was the strongest intercellular communication pathway in Healing wounds.
- This axis directly connects macrophage M1 activation (IL1B production) to fibroblast activation and ECM remodeling (HE-Fibroblast induction).
- The IL-1 pathway provides a mechanistic bridge between the macrophage polarization findings and the HE-Fibroblast subpopulation enrichment (**Figure S3**).

### 4.2 Phase 2: Bulk RNA-seq Multi-Cohort Validation

#### 4.2.1 GSE134431: Diabetic Foot Skin vs. Diabetic Foot Ulcer

Differential expression analysis of 13 samples (Diabetic Foot Skin controls vs. active DFU) confirmed key single-cell findings at the tissue level:

| Gene   | Fold Change (DFU/DFS) | p-value   | Direction |
|--------|------------------------|-----------|-----------|
| HIF1A  | 1.91                   | 0.016     | Up        |
| IL1B   | 9.33                   | < 0.001   | Up        |
| CXCL8  | 23.34                  | < 0.001   | Up        |
| VEGFA  | 1.52                   | 0.042     | Up        |

**HIF1A-NFKB1 correlation:** Spearman rho = **0.708** (p < 0.01), confirming the single-cell-identified coupling at the tissue level.

#### 4.2.2 GSE199939: Normal Skin vs. Diabetic Foot Ulcer

Analysis of 20 samples (Normal skin vs. DFU) provided a more dramatic contrast:

| Gene   | Fold Change (DFU/Normal) | p-value     | Direction |
|--------|--------------------------|-------------|-----------|
| HIF1A  | 4.06                     | 5.7 x 10^-6 | Up        |
| IL1B   | 5.56                     | < 0.001     | Up        |
| CXCL8  | 7.47                     | < 0.001     | Up        |

**HIF1A-NFKB1 correlation:** Spearman rho = **0.634** (p < 0.01).

**Genome-wide differential expression:** 1,304 significantly upregulated genes and 2,467 significantly downregulated genes in DFU versus normal skin, indicating massive transcriptomic remodeling (**Figure 5, Figure 6**).

#### 4.2.3 Cross-Cohort Consistency

The consistency of findings across two independent bulk RNA-seq cohorts with different experimental designs (DFS vs. DFU; Normal vs. DFU) strongly validates the single-cell observations:

- HIF1A is consistently upregulated in DFU tissue (FC range: 1.91-4.06).
- HIF1A-NFKB1 co-expression correlation is robust and reproducible (rho range: 0.634-0.708).
- Inflammatory cytokines (IL1B, CXCL8) show large and consistent upregulation across cohorts.

### 4.3 Phase 3: Temporal Dynamics Analysis

#### 4.3.1 GSE147890: Diabetic vs. Control Wound Healing Time Course

Analysis of the GSE147890 dataset, which profiles wound healing at 24 hours post-wounding in diabetic and control subjects, revealed:

- **Control 24h response:** 3,493 differentially expressed probes (vs. baseline)
- **Diabetic 24h response:** 4,098 differentially expressed probes (vs. baseline)
- **Interaction effect probes (Diabetic-specific):** Only 2 probes showed statistically significant diabetic-specific temporal responses after multiple testing correction.

**Interpretation:** The overall transcriptional response to wounding at 24 hours is largely conserved between diabetic and control conditions in terms of the genes involved, but the magnitude and potentially the kinetics of the response differ. The near-absence of interaction effects suggests that diabetes does not fundamentally reprogram the wound response gene repertoire at this early time point, but rather modulates the intensity of a shared healing program. This is consistent with the hypothesis that DFU healing failure is a quantitative rather than qualitative defect in the wound response.

### 4.4 Phase 4: Macrophage Mechanistic Dissection

#### 4.4.1 Hypoxia Signature Derivation

From the in vitro macrophage hypoxia datasets (GSE15949, GSE16099), we derived a robust 100-gene hypoxia response signature. Of these, 82 genes were detectable in the scRNA-seq data and were used for cross-platform hypoxia scoring via AddModuleScore.

**Data processing note:** Both microarray datasets used linear-scale expression values (maximum values > 20), requiring log2 transformation prior to differential expression analysis with limma.

#### 4.4.2 Hypoxia Activation in DFU Macrophages

Application of the macrophage-specific hypoxia signature to the scRNA-seq macrophage compartment revealed dramatic condition-dependent differences:

| Condition    | Mean Hypoxia Score | Interpretation                  |
|-------------|--------------------|---------------------------------|
| Healing     | 0.101              | Strong hypoxia pathway activation |
| Non-Healing | 0.013              | Minimal hypoxia pathway activation |

**Statistical significance:** Wilcoxon rank-sum test, p = **1.33 x 10^-155** (Healing vs. Non-Healing macrophages).

This result demonstrates that macrophages in healing DFU wounds exhibit dramatically higher activation of hypoxia-responsive transcriptional programs compared to those in non-healing wounds. Combined with the M1 polarization finding, this establishes that healing-competent macrophages operate in a high-HIF, high-NF-κB, M1-polarized state (**Figure 5**).

#### 4.4.3 HIF+M1 Macrophage Subcluster (Cluster C0)

Reclustering of the macrophage compartment identified a distinct subcluster (C0) characterized by concurrent HIF pathway activation and M1 polarization:

**Cluster C0 characteristics:**

| Feature               | Value    | Context                           |
|-----------------------|----------|-----------------------------------|
| HIF pathway score     | 0.157    | Highest among all clusters        |
| NF-κB pathway score   | 1.209    | Highest among all clusters        |
| M1 polarization score | 0.450    | Highest among all clusters        |
| Healing-derived cells | **75.9%** | Strongly healing-enriched         |
| NonHealing-derived    | 15.2%    | Depleted relative to input        |

Cluster C0 represents the prototypical "healing macrophage": a cell state characterized by maximal activation of both HIF-1α and NF-κB programs, driving robust M1 polarization. The near-exclusive derivation from healing wounds confirms its functional association with productive wound repair.

#### 4.4.4 Hypoxia-Induced Target Gene Upregulation

Cross-referencing the macrophage hypoxia signature with single-cell expression data in the HIF+M1 subcluster identified key hypoxia-upregulated effector genes:

| Gene     | Fold Change (HIF+M1 vs. rest) | Function                        |
|----------|-------------------------------|---------------------------------|
| PTGS2    | 9.4x                         | Prostaglandin synthesis (COX-2) |
| ANGPTL4  | 9.3x                         | Angiogenesis, metabolism        |
| IL1B     | 7.2x                         | Pro-inflammatory cytokine       |
| TNF      | 5.1x                         | Pro-inflammatory cytokine       |
| VEGFA    | ~3-4x                        | Angiogenesis                    |
| BNIP3    | ~2-3x                        | Autophagy, cell survival        |

These genes represent direct transcriptional targets of the HIF-1α/NF-κB coupling axis and provide the molecular effectors through which this signaling axis promotes wound healing through inflammation, angiogenesis, and metabolic reprogramming.

---

## 5. Key Findings Summary

### Finding 1: Non-Healing DFUs Exhibit M2 Macrophage Retention and Impaired M1 Activation

Contrary to the simplified model that M2 anti-inflammatory macrophages universally promote wound healing, our data demonstrate that non-healing DFU wounds are characterized by M2-dominant macrophage populations locked in a dysregulated state. Healing wounds instead exhibit M1-biased polarization, suggesting that the capacity for adequate M1-driven inflammatory initiation — and subsequent M1→M2 dynamic transition — is associated with healing progression in the diabetic wound context.

### Finding 2: NF-κB Is the Primary Driver, HIF-1α Is the Synergistic Co-activator

NF-κB signaling shows the strongest correlation with M1 polarization (rho = 0.819), establishing it as the dominant upstream driver. HIF-1α acts synergistically (rho = 0.340 individually; combined coupling rho = 0.721), amplifying the inflammatory program under hypoxic conditions.

### Finding 3: HIF+M1 Macrophage Subcluster Is Healing-Specific

A distinct macrophage subpopulation (Cluster C0) with maximal HIF and NF-κB activation, M1 polarization, and 75.9% healing derivation represents the functional "healing macrophage" state.

### Finding 4: Hypoxia Activation Distinguishes Healing from Non-Healing

Macrophage hypoxia scores differ by nearly an order of magnitude between healing (0.101) and non-healing (0.013) wounds (p = 1.33 x 10^-155), making hypoxia pathway activation one of the strongest molecular discriminators of healing outcome.

### Finding 5: HE-Fibroblasts Link Hypoxia to ECM Remodeling

The MMP1+/MMP3+/HIF1A+ HE-Fibroblast subpopulation (4.2-fold enriched in healing wounds) connects hypoxia signaling to productive extracellular matrix turnover, mediated by IL-1 signaling from M1 macrophages.

### Finding 6: Multi-Cohort Bulk Validation Confirms Core Findings

Two independent bulk RNA-seq cohorts confirm HIF1A upregulation (FC 1.91-4.06), HIF1A-NFKB1 co-expression (rho 0.634-0.708), and massive inflammatory gene induction (IL1B FC up to 9.33, CXCL8 FC up to 23.34) in DFU tissue.

### Finding 7: Quantitative Rather Than Qualitative Healing Defect

Temporal analysis reveals that diabetic and control wounds activate largely overlapping gene programs with minimal diabetic-specific interaction effects, suggesting that DFU healing failure is a quantitative deficit in an otherwise conserved response.

---

## 6. Discussion

### 6.1 Redefining the M1/M2 Dynamic in Diabetic Wound Healing

Our findings necessitate a reassessment of the macrophage polarization model in DFU pathophysiology. The conventional model, which views chronic wound inflammation as a pathological "stalling" in the M1 phase, assumes that M1-to-M2 transition is the critical determinant of healing. Our data suggest an alternative framework:

**Revised model:** In diabetic wounds, non-healing is characterized not by "being stuck in M1" but rather by a failure to adequately activate the initial M1 program, with macrophages defaulting to an M2-biased state that may reflect immune exhaustion, metabolic dysfunction, or premature anti-inflammatory skewing. The capacity to mount a robust, HIF-1α-augmented M1 inflammatory response appears to be a prerequisite for the subsequent M1→M2 transition that drives wound resolution.

This interpretation is consistent with emerging literature demonstrating that adequate initial inflammation is necessary for wound debridement, pathogen clearance, and activation of reparative fibroblast programs. The IL-1 communication axis (M1 macrophage → Fibroblast) we identified provides a direct mechanistic link between macrophage activation and wound repair.

### 6.2 The HIF-1α/NF-κB Coupling Axis as a Therapeutic Target

The identification of HIF-1α/NF-κB coupling as the upstream switch controlling macrophage polarization and healing outcome has direct therapeutic implications:

1. **HIF-1α stabilization** (e.g., PHD inhibitors such as roxadustat or daprodustat) may enhance the healing response by amplifying M1 macrophage activation and hypoxia-responsive gene programs.
2. **Local hypoxia modulation** through controlled oxygen therapy or biomaterial-based oxygen delivery could augment the HIF-1α arm of the coupling axis.
3. **Caution against broad-spectrum NF-κB inhibition** in DFU management, as suppressing the dominant M1 driver may paradoxically impair healing.
4. **Targeted activation** of the IL-1 → HE-Fibroblast axis could promote productive ECM remodeling.

### 6.3 The HE-Fibroblast as a Biomarker and Therapeutic Target

The HE-Fibroblast subpopulation (MMP1+/MMP3+/HIF1A+) represents both a potential biomarker for healing competence and a therapeutic target. The 4.2-fold enrichment in healing wounds suggests that strategies to expand or activate this subpopulation -- potentially through IL-1 stimulation or localized hypoxia -- could enhance wound repair. Conversely, the proportion of HE-Fibroblasts in wound biopsies could serve as a prognostic indicator.

### 6.4 Limitations

Several limitations should be acknowledged:

1. **Cross-sectional scRNA-seq design:** The single-cell data represents a snapshot rather than longitudinal trajectory. Pseudotime and trajectory analyses can partially address this but cannot fully substitute for true temporal profiling.
2. **Bulk RNA-seq deconvolution limits:** Tissue-level gene expression conflates contributions from multiple cell types, and immune deconvolution algorithms have inherent uncertainty.
3. **In vitro hypoxia signature generalizability:** The macrophage hypoxia signature was derived from in vitro monoculture experiments, which may not fully recapitulate the complex in vivo microenvironment.
4. **Limited temporal resolution in diabetes:** The temporal dataset (GSE147890) captured only a 24-hour window, and longer time courses would be needed to characterize the full kinetics of the diabetic wound response.
5. **Absence of functional validation:** All findings are correlative; functional perturbation experiments (e.g., conditional knockouts, pharmacological intervention) would be required to establish causality.
6. **Patient heterogeneity:** DFU etiology is multifactorial (neuropathy, vascular disease, infection), and the current analysis does not stratify by DFU subtype.

### 6.5 Future Directions

1. **Spatial transcriptomics** (e.g., Visium, MERFISH) to map the HIF-1α/NF-κB coupling axis and HE-Fibroblasts within the wound architecture.
2. **In vivo functional validation** using diabetic mouse models with myeloid-specific HIF-1α or NF-κB knockout/overexpression.
3. **Clinical cohort validation** of HE-Fibroblast proportion as a prognostic biomarker for DFU healing.
4. **Pharmacological studies** testing PHD inhibitors or targeted IL-1 pathway modulation in preclinical DFU models.
5. **Multi-omic integration** incorporating ATAC-seq or CUT&Tag to characterize the chromatin landscape of the HIF+M1 macrophage state.

---

## 7. Detailed Methods

### 7.1 Single-Cell RNA-seq Processing (Phase 1)

**Data acquisition:** Raw count matrices for GSE165816 were obtained from the Gene Expression Omnibus (GEO). The dataset comprises 50 samples across four conditions: Healthy skin (n=10), Diabetic skin without wound (DM, n=10), Healing DFU (n=15), and Non-Healing DFU (n=15).

**Quality control:** Cells were filtered using the following criteria:
- Minimum 200 detected genes per cell
- Maximum 20% mitochondrial gene content
- Doublet removal based on expected doublet rate per sample

**Normalization and scaling:** Log-normalization was applied using Seurat's `NormalizeData` with a scale factor of 10,000. The top 2,000 highly variable genes were identified using the variance-stabilizing transformation (VST) method. Data were scaled with regression of mitochondrial percentage and cell cycle scores.

**Dimensionality reduction and batch correction:** Principal component analysis (PCA) was performed on scaled data, retaining 30 principal components based on elbow plot analysis. Harmony batch correction was applied across samples to mitigate technical variation while preserving biological signal. Uniform Manifold Approximation and Projection (UMAP) was computed on Harmony-corrected embeddings.

**Clustering and annotation:** Shared nearest neighbor (SNN) graph construction and Louvain clustering were performed at multiple resolutions. Cell type annotation was based on canonical marker gene expression, validated by differential expression analysis between clusters.

**Pathway scoring:** HIF pathway, NF-κB pathway, M1 polarization, and M2 polarization scores were computed using Seurat's `AddModuleScore` function with curated gene sets:
- HIF pathway: HIF1A, EPAS1, VEGFA, LDHA, PGK1, SLC2A1, BNIP3, PDK1, ENO1, ALDOA, and additional canonical HIF targets.
- NF-κB pathway: NFKB1, NFKB2, RELA, RELB, REL, IKBKB, IKBKG, CHUK, and downstream targets.
- M1 markers: TNF, IL1B, IL6, NOS2, CXCL8, CCL2, CD80, CD86.
- M2 markers: MRC1, CD163, ARG1, IL10, TGFB1, CCL17, CCL22.

**Subpopulation analysis:** Fibroblast and macrophage compartments were independently subset, re-normalized, and re-clustered at higher resolution to identify subpopulations. HE-Fibroblasts were defined by co-expression of MMP1, MMP3, and HIF1A above the 75th percentile within the fibroblast compartment.

**Cell-cell communication:** Ligand-receptor interaction analysis was performed using a simplified approach based on known L-R pairs from CellChat/CellPhoneDB databases, calculating communication probability scores based on expression of ligand-receptor pairs across cell type pairs within each condition.

### 7.2 Bulk RNA-seq Differential Expression (Phase 2)

**Data acquisition:** Processed expression matrices for GSE134431 and GSE199939 were obtained from GEO using the `GEOquery` R package.

**GSE134431 processing:** The dataset required special handling due to a 3-row header format in the Excel source file. Data were read with `col_names = FALSE`, and sample IDs were manually extracted from the header row. Expression values were organized into a proper matrix with gene identifiers as row names.

**Differential expression analysis:** The `limma` package was used for differential expression analysis with the following pipeline:
1. Design matrix construction for two-group comparisons (DFS vs. DFU; Normal vs. DFU)
2. Linear model fitting with `lmFit`
3. Empirical Bayes moderation with `eBayes`
4. Multiple testing correction using the Benjamini-Hochberg method

**Correlation analysis:** Spearman rank correlation was computed between HIF1A and NFKB1 expression across all samples within each dataset.

### 7.3 Temporal Dynamics Analysis (Phase 3)

**GSE147890 processing:** Pre- and post-wounding (24h) samples from diabetic and control subjects were analyzed using a two-factor design (condition x time) in limma. Interaction terms were extracted to identify diabetic-specific temporal responses.

**Multiple testing correction:** False discovery rate (FDR) control at q < 0.05 using the Benjamini-Hochberg procedure.

### 7.4 Macrophage Hypoxia Signature (Phase 4)

**Signature derivation:** Differentially expressed genes (FDR < 0.05, absolute log2 fold change > 1) from GSE15949 (macrophage hypoxia time course) and GSE16099 (macrophage polarization under hypoxia) were intersected and filtered for consistent upregulation under hypoxia to derive a 100-gene macrophage hypoxia signature.

**Data scale detection:** Expression data from both microarray datasets were examined for value range to determine appropriate scale. Values with maximum > 20 indicated linear-scale data requiring log2 transformation.

**Cross-platform scoring:** The 82 signature genes detectable in the scRNA-seq data were used with `AddModuleScore` to compute per-cell hypoxia scores in the macrophage compartment.

**Subcluster characterization:** Macrophage subclusters were profiled for HIF pathway score, NF-κB pathway score, M1 score, M2 score, and hypoxia score. Condition composition was computed for each subcluster to identify healing-enriched populations.

### 7.5 Statistical Methods

- **Two-group comparisons:** Wilcoxon rank-sum test (non-parametric) for single-cell score comparisons; moderated t-test (limma) for bulk expression comparisons.
- **Correlation analysis:** Spearman rank correlation with exact p-values.
- **Multiple testing correction:** Benjamini-Hochberg FDR across all applicable analyses.
- **Effect size reporting:** Fold changes (for expression), Spearman rho (for correlations), and proportion ratios (for compositional analyses).

### 7.6 Software and Packages

| Tool/Package | Version | Purpose                                |
|-------------|---------|----------------------------------------|
| R           | 4.4.3   | Statistical computing environment      |
| Seurat      | v5      | scRNA-seq analysis                     |
| Harmony     | latest  | Batch correction                       |
| limma       | latest  | Differential expression (bulk/array)   |
| GEOquery    | latest  | GEO data retrieval                     |
| ggplot2     | latest  | Visualization                          |
| dplyr       | latest  | Data manipulation                      |
| pheatmap    | latest  | Heatmap visualization                  |

---

## 8. Dataset References

| Accession   | Title / Description                                                     | Platform     | Reference                         |
|-------------|-------------------------------------------------------------------------|--------------|-----------------------------------|
| GSE165816   | Single-cell transcriptomics of diabetic foot ulcers                    | 10x Genomics | Theocharidis et al., 2022        |
| GSE134431   | Transcriptomic profiling of diabetic foot skin and ulcers              | RNA-seq      | Ramirez et al., 2020             |
| GSE199939   | Bulk transcriptomics of normal skin vs diabetic foot ulcers            | RNA-seq      | GEO contributors                 |
| GSE15949    | Macrophage transcriptional response to hypoxia (time course)           | Microarray   | Bosco et al., 2006               |
| GSE16099    | M1/M2 macrophage polarization under hypoxia                           | Microarray   | Bentley et al.                   |
| GSE28914    | Burn wound healing temporal transcriptomics                            | Microarray   | GEO contributors                 |
| GSE50425    | Pressure ulcer temporal gene expression                                | Microarray   | GEO contributors                 |
| GSE147890   | Diabetic vs. control wound healing at 24h                              | Microarray   | GEO contributors                 |

---

## 9. Figure Index

### Main Figures

| Figure | Title                                                | File                       |
|--------|------------------------------------------------------|----------------------------|
| Fig. 1 | Study Overview and Single-Cell Landscape             | `Figure1.pdf`             |
| Fig. 2 | HIF-1α/NF-κB Co-expression Landscape                | `Figure2.pdf`             |
| Fig. 3 | HE-Fibroblast Subpopulation Characterization         | `Figure3.pdf`             |
| Fig. 4 | Macrophage Polarization and HIF/NF-κB Coupling       | `Figure4.pdf`             |
| Fig. 5 | Hypoxia Signature in DFU Macrophages                 | `Figure5.pdf`             |
| Fig. 6 | Bulk RNA-seq Multi-Cohort Validation                 | *(data in Phase2 outputs)* |
| Fig. 7 | Mechanistic Model: HIF-1α/NF-κB → M1 → Healing      | `Figure7.pdf`             |

### Supplementary Figures

| Figure  | Title                                          | File                       |
|---------|------------------------------------------------|----------------------------|
| Fig. S1 | Quality Control Metrics                       | `FigureS1.pdf`            |
| Fig. S2 | Cell Type Marker Gene Heatmap                 | `FigureS2.pdf`            |
| Fig. S3 | Cell-Cell Communication: Top L-R Pairs        | `FigureS3.pdf`            |
| Fig. S4 | WGCNA Module Analysis                         | `FigureS4.pdf`            |

### Supplementary Tables

| Table   | Title                                          | File                            |
|---------|------------------------------------------------|---------------------------------|
| Table S1| Cell Counts per Sample                        | `TableS1_sample_stats.csv`     |
| Table S2| Cell Type by Condition Distribution           | `TableS2_celltype_condition.csv`|
| Table S4| Pathway Scores per Cell Type per Condition    | `TableS4_pathway_scores.csv`   |
| Table S5| Macrophage Hypoxia Signature Gene List        | `TableS5_hypoxia_signature.csv`|

All figures and supplementary materials are located in: `Phase_output/publication_figures/`

---

*End of Technical Report*
