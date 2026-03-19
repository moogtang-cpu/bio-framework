# Single-Cell Atlas of Human Endometrial Basalis Stem Cells Reveals Stem Cell Depletion and Niche Remodeling in Asherman's Syndrome with Identification of Drug Intervention Targets

## Abstract

The human endometrium undergoes cyclical regeneration driven by basalis stem/progenitor cells, yet the molecular mechanisms underlying regenerative failure in Asherman's syndrome (AS) remain poorly understood. Here, we construct a comprehensive single-cell atlas of 314,805 cells by integrating five scRNA-seq datasets spanning normal cycling endometrium, AS, and other reproductive conditions. We identify 14 cell types including two distinct stem/progenitor populations: SOX9+LGR5+ stem cells (9,489 cells, 3.0%) and CD133+ progenitors (11,466 cells, 3.6%). Diffusion pseudotime analysis reveals a hierarchical differentiation trajectory from SOX9+LGR5+ stem cells through CD133+ progenitors to mature glandular and secretory epithelium, confirmed by CytoTRACE stemness scoring (SOX9+LGR5+: 0.342, highest among all types). Spatial transcriptomics deconvolution validates stem cell enrichment in the basalis niche (SOX9+ 4.31%). In AS, we discover systematic stem cell dysfunction encompassing CD133+ progenitor depletion (~11-fold reduction, padj=0.011), dramatic stemness loss in remaining SOX9+ cells (CytoTRACE 0.334→0.072, p=1.78×10⁻⁴⁸), and comprehensive silencing of WNT, NOTCH, TGF-β, FGF, EGF, Hedgehog, cytokine, and chemokine signaling pathways. Concurrently, basalis fibroblasts undergo fibrotic reprogramming with paradoxical upregulation of all major signaling pathways, transforming from niche supporters to fibrosis drivers. Multi-evidence integration of 8,051 differentially expressed genes, 760 transcription factor activities, and 56,760 cell-cell interactions identifies 463 drug target candidates, of which 60 have known drug interactions from DGIdb. Among the top targets, MET (21 approved drugs), CXCL8 (43 approved drugs), and PGR (20 approved drugs) show strong multi-omic support and are validated for expression in AS patient-derived endometrial organoids. Our atlas provides a comprehensive resource for understanding endometrial regeneration and identifies actionable therapeutic targets for regenerative disorders.

**Keywords**: endometrial stem cells, Asherman's syndrome, single-cell atlas, drug targets, stem cell niche, cell-cell communication, transcription factor, regenerative medicine

---

## Introduction

The human endometrium is a uniquely regenerative tissue that undergoes approximately 400 cycles of shedding and regeneration during a woman's reproductive lifetime¹. This cyclical process, driven by ovarian hormones, involves the complete loss and subsequent regeneration of the functional layer (functionalis), while the deeper basalis layer persists to serve as the source of regenerative cells²,³. The remarkable regenerative capacity of the endometrium has long been attributed to resident stem/progenitor cells, though the precise identity and molecular characteristics of these cells have remained subjects of active investigation.

Several candidate markers for endometrial stem/progenitor cells have been identified over the past two decades. Gargett and colleagues first proposed the existence of clonogenic endometrial stem cells based on colony-forming efficiency assays⁴. Subsequent studies identified CD146+PDGFRβ+ mesenchymal stem-like cells in the perivascular niche⁵, while epithelial progenitors were characterized by CD133 (PROM1) expression⁶,⁷. SOX9, a master transcription factor for tissue stem cells in multiple organs including the intestine and hair follicle, was found to mark a basalis glandular epithelial population with stem cell properties⁸,⁹. More recently, the Human Endometrial Cell Atlas (HECA) identified SOX9+ populations with distinct localization in the basalis and functionalis layers, though the precise relationship between SOX9+, LGR5+, and CD133+ populations remains unclear¹⁰.

LGR5, a well-established intestinal stem cell marker and Wnt target gene, has been detected in the human endometrium, particularly in glandular epithelium¹¹. However, unlike the intestinal crypt where LGR5+ cells represent a well-defined stem cell population, the role of LGR5 in endometrial stem cell hierarchy has been debated, with some studies reporting its expression in both stem and differentiated compartments¹²,¹³. Understanding the relationship between SOX9+, LGR5+, and CD133+ populations is critical for defining the endometrial stem cell hierarchy and understanding regenerative failure in disease.

The endometrial stem cell niche is increasingly recognized as a complex microenvironment involving paracrine signals from surrounding stromal cells, immune cells, and the extracellular matrix¹⁴. Basalis fibroblasts, perivascular cells, and endothelial cells form a supportive niche that maintains stem cell self-renewal through WNT, NOTCH, and growth factor signaling¹⁵,¹⁶. Disruption of this niche may contribute to regenerative disorders, yet the molecular basis of niche dysfunction remains poorly characterized.

Asherman's syndrome (AS), also known as intrauterine adhesions (IUA), is a condition resulting from endometrial damage, most commonly following dilatation and curettage or other uterine surgical procedures¹⁷,¹⁸. AS is characterized by partial or complete obliteration of the uterine cavity by fibrous adhesions, leading to amenorrhea, hypomenorrhea, recurrent pregnancy loss, and infertility¹⁹. The prevalence of IUA following uterine surgery ranges from 2% to 45% depending on the procedure and assessment method²⁰. Current treatments, primarily hysteroscopic adhesiolysis, have limited long-term efficacy with recurrence rates of 20–62%²¹, highlighting the need for therapeutic strategies that address the underlying regenerative failure.

Recent single-cell transcriptomic studies have begun to characterize the endometrial cellular landscape in AS. Santamaria et al. performed single-cell RNA sequencing of endometrial samples from AS patients, revealing altered cellular composition and gene expression patterns²². However, a systematic, multi-dataset analysis of stem cell dysfunction — integrating trajectory analysis, transcription factor regulation, and cell-cell communication — with an explicit focus on therapeutic target discovery has not been performed.

Here, we present an integrative single-cell atlas of 314,805 human endometrial cells from five independent scRNA-seq datasets. By combining diffusion pseudotime analysis, transcription factor activity inference, and consensus-based cell communication profiling, we elucidate the mechanisms of stem cell depletion and niche remodeling in AS. Spatial transcriptomics deconvolution provides independent validation of stem cell localization and disease-associated changes. Through multi-evidence target prioritization and drug database mining, we identify actionable therapeutic targets with existing approved drugs, providing a roadmap for pharmacological intervention in endometrial regenerative disorders.

---

## Results

### Construction of a high-resolution human endometrial cell atlas

To construct a comprehensive endometrial cell atlas, we integrated five scRNA-seq datasets spanning diverse reproductive conditions (Figure 1A; Supplementary Table S1). These included: GSE215968 (AS vs. WOI controls; 106,400 cells)²², GSE215968_CD133 (CD133+-enriched cells; 69,701 cells), GSE111976 (normal menstrual cycle; 61,503 cells)²³, E-MTAB-10287 (temporal dynamics; 49,550 cells)¹⁰, and GSE260658 (endometrium and decidua; 27,651 cells). After independent quality control of each dataset (Methods), a total of 314,805 cells with 18,365 common genes were retained for integration.

Batch correction was performed using Harmony²⁴ with GPU acceleration, which converged in three iterations. Integration quality was assessed by batch-corrected average silhouette width (ASW_batch = 0.024, indicating effective batch mixing), cluster mixing score (0.565), and visual inspection of UMAP embeddings (Figure 1B; Supplementary Figure S2). Five dataset-dominant clusters (>70% from a single dataset) were identified and confirmed to reflect biological differences rather than technical artifacts: three clusters dominated by CD133+-sorted cells and two by AS-specific cell states.

Cell type annotation employed a two-stage approach: KNN-based label transfer from the HECA reference atlas (313,527 cells, 36 cell types; k=30, distance-weighted)¹⁰ followed by marker-guided manual refinement. The reference mapping achieved a mean cell type confidence of 0.866 and lineage confidence of 0.994, with 94% of cells exceeding a confidence threshold of 0.5. Manual refinement using canonical marker genes (Supplementary Table S2) consolidated the annotation into 14 cell types across five major lineages (Figure 1C-E): mesenchymal (39.9%, including decidualized stromal cells, smooth muscle, basalis fibroblasts, and pericytes), immune (31.0%, including lymphoid immune cells, NKT cells, macrophages, B cells, and erythrocytes), epithelial (16.3%, including glandular and secretory glandular epithelium), stem/progenitor (6.7%, including SOX9+LGR5+ stem cells and CD133+ progenitors), and endothelial (4.6%, lymphatic endothelium).

### Identification of two distinct stem/progenitor populations with hierarchical organization

Among the 14 annotated cell types, we identified two stem/progenitor populations with distinct molecular signatures (Figure 2A-F). SOX9+LGR5+ stem cells (cluster 12; 9,489 cells, 3.0%) exhibited the highest expression of the stem cell transcription factor SOX9 (48.8% of cells), along with LGR5 (25.0%) and PROM1/CD133 (26.4%). This population was primarily derived from GSE111976 (normal cycling endometrium) and GSE215968 (both control and AS samples), suggesting its presence across physiological and pathological states.

CD133+ progenitors (clusters 18, 22, 28; 11,466 cells, 3.6%) were distinguished by high PROM1 expression (47.5%) but substantially lower SOX9 (9.5%) and LGR5 (13.0%), consistent with a more committed progenitor state that has partially lost stem cell transcription factor expression. This population was enriched in the GSE215968_CD133 dataset (CD133+ FACS-sorted), validating the annotation.

Basalis fibroblasts (cluster 24; 2,168 cells, 0.7%) formed a third critical population. While not stem cells themselves, basalis fibroblasts expressed niche-associated markers consistent with the HECA Fibroblast_basalis annotation and have been implicated as key regulators of the stem cell microenvironment¹⁰.

CytoTRACE stemness scoring²⁵, based on gene detection counts as a proxy for developmental potential, confirmed the hierarchical organization: SOX9+LGR5+ stem cells scored highest (0.342), followed by CD133+ progenitors (0.304), basalis fibroblasts (0.290), and differentiated cell types (Glandular epithelium: 0.255; Decidualized stromal: 0.177) (Figure 2E-F). These scores are consistent with SOX9+LGR5+ cells representing the most primitive endometrial stem cell population.

### Diffusion pseudotime reveals continuous differentiation trajectories

To map differentiation trajectories, we performed diffusion pseudotime (DPT) analysis²⁶ on two lineage subsets (Figure 3D-F; Supplementary Figure S5).

**Epithelial lineage** (72,417 cells): DPT analysis with SOX9+LGR5+ stem cells as root revealed a continuous trajectory through CD133+ progenitors (median pseudotime = 0.516) toward mature glandular (0.122) and secretory glandular (0.127) epithelium. PAGA connectivity analysis²⁷ confirmed direct connections between SOX9+LGR5+ stem cells and CD133+ progenitors, and between progenitors and differentiated epithelial types. The pseudotime ordering was consistent with progressive loss of SOX9 and gain of differentiation markers along the trajectory.

**Mesenchymal lineage** (134,504 cells): Basalis fibroblasts occupied the earliest pseudotime position (median = 0.018), supporting their role as niche progenitors. The trajectory extended through decidualized stromal cells (0.042) to smooth muscle (0.097), consistent with known mesenchymal differentiation patterns in the endometrium.

### Spatial transcriptomics validates stem cell enrichment in the basalis niche

To validate the spatial organization of endometrial stem cells, we performed Cell2location deconvolution²⁸ of 10x Visium spatial transcriptomics data from GSE287278 (10,135 spots across 8 samples: 4 controls, 4 RIF patients) (Figure 3A-C; Supplementary Figure S4).

Unsupervised clustering of the spatial data identified three major tissue regions: stroma (6,480 spots, 63.9%), functionalis (3,423 spots, 33.8%), and basalis niche (232 spots, 2.3%). These regions were annotated based on composite scores of region-specific marker genes (Methods). Cell2location deconvolution using reference signatures from our integrated atlas (60,445 cells subsampled, 13,017 common genes) mapped all 14 cell types to their spatial locations.

SOX9+LGR5+ stem cells showed preferential enrichment in the basalis niche region (4.31% abundance) compared to stroma (2.84%) and functionalis (2.67%). CD133+ progenitors exhibited even stronger enrichment in the basalis niche (9.46-13.34%). In RIF samples, SOX9+ stem cell abundance was significantly reduced compared to controls (2.48% vs. 3.26%, p < 0.001, Mann-Whitney U test), while CD133+ progenitors were paradoxically increased (5.86% vs. 3.47%, p < 0.001), suggesting compensatory progenitor expansion in the setting of impaired regeneration.

### Systematic stem cell depletion and functional impairment in Asherman's syndrome

To characterize stem cell dysfunction in AS, we compared cell composition and gene expression between WOI controls (12 samples, 66,064 cells) and AS patients (16 samples, 40,336 cells) from GSE215968, which minimized inter-study batch effects (Figure 4A).

**Cell proportion changes**: Per-sample cell type proportions were compared using Wilcoxon rank-sum tests with FDR correction (Methods; Supplementary Table S3). Four cell types showed significant changes (padj < 0.05): CD133+ progenitors were depleted ~11-fold (7.7% → 0.68%, log2FC = −3.51, padj = 0.011), secretory glandular cells decreased ~9-fold (16.1% → 1.8%, log2FC = −3.13, padj = 0.025), while NKT cells (4.2% → 19.1%, log2FC = +2.18, padj = 0.035) and macrophages (0.5% → 6.8%, log2FC = +3.71, padj = 0.035) expanded significantly. SOX9+LGR5+ stem cells showed a trending decrease (6.5% → 2.4%, padj = 0.21) that did not reach statistical significance.

**Stemness loss**: Beyond numerical depletion, remaining SOX9+LGR5+ stem cells in AS exhibited profound functional impairment (Figure 4B). CytoTRACE stemness scores dropped by 78% (Control: 0.334 → AS: 0.072, Kolmogorov-Smirnov test, p = 1.78×10⁻⁴⁸), indicating near-complete loss of stem cell identity. Pseudotime distributions were significantly shifted across all cell types (KS test, all p < 0.05), suggesting globally disrupted differentiation dynamics.

**Differential gene expression**: Pseudobulk DEG analysis using PyDESeq2²⁹ (per cell type, per sample aggregation) identified 8,051 differentially expressed genes across 10 cell types (padj < 0.05, |log2FC| > 0.25; Supplementary Table S3). Decidualized stromal cells were most affected (3,875 DEGs: 1,698 up, 2,177 down), followed by lymphatic endothelium (1,533), smooth muscle (729), and macrophages (542). CD133+ progenitors showed 223 DEGs with a striking bias toward downregulation (14 up, 209 down), suggesting widespread transcriptional silencing. SOX9+LGR5+ stem cell DEG analysis was limited by insufficient per-sample cell counts for pseudobulk aggregation (only 1 control and 1 AS sample had adequate representation).

**Pathway enrichment**: Gene set enrichment analysis (GSEA) using preranked gene lists identified significant enrichment of WNT, NOTCH, TGF-β, fibrosis, and stem cell niche gene sets across multiple cell types (Supplementary Table S4; Supplementary Figure S7). Custom pathway gene sets for endometrial stem cell biology confirmed downregulation of regeneration-associated pathways and upregulation of fibrosis markers in AS.

### Comprehensive signaling pathway silencing in AS stem cells

To quantify signaling pathway changes at the single-cell level, we computed pathway activity scores based on mean expression of curated gene sets for eight major signaling pathways (Methods) (Figure 4D).

SOX9+LGR5+ stem cells in AS exhibited significant downregulation of all eight pathways examined: WNT (diff = −0.012, padj = 0.010), NOTCH (diff = −0.042, padj = 2.0×10⁻¹²), TGF-β/BMP (diff = −0.039, padj = 1.1×10⁻¹²), FGF (diff = −0.027, padj = 1.1×10⁻⁵), EGF (diff = −0.073, padj = 1.7×10⁻¹²), Hedgehog (diff = −0.061, padj = 4.6×10⁻¹⁵), cytokine (diff = −0.119, padj = 6.0×10⁻²⁴), and chemokine (diff = −0.094, padj = 7.3×10⁻¹²) signaling. This comprehensive silencing suggests a global loss of niche signaling reception rather than a pathway-specific defect.

In striking contrast, basalis fibroblasts showed significant upregulation of WNT (padj = 1.5×10⁻¹⁶), NOTCH (padj = 9.3×10⁻²⁸), Hedgehog (padj = 6.1×10⁻¹⁴), TGF-β/BMP (padj = 5.6×10⁻⁵), and EGF (padj = 1.7×10⁻¹²) pathways. This paradoxical pattern — pathway silencing in stem cells coupled with activation in niche fibroblasts — suggests a fundamental transformation of the basalis niche from a regenerative to a fibrotic microenvironment.

CD133+ progenitors showed an intermediate pattern with downregulation of chemokine (padj = 1.8×10⁻⁴⁶), FGF (padj = 5.8×10⁻¹⁵), and EGF (padj = 2.1×10⁻⁵) pathways, but upregulation of NOTCH (padj = 7.7×10⁻⁴) and WNT (padj = 0.012), possibly reflecting compensatory signaling activation in surviving progenitors.

### Transcription factor dysregulation reveals lineage-specific regulatory changes

To identify upstream regulators of the observed transcriptional changes, we inferred transcription factor (TF) activities using the Univariate Linear Model (ULM)³⁰ with the CollecTRI regulatory network from OmniPath³¹ (62,411 TF-target interactions, 760 TFs) (Figure 5A-B).

In SOX9+LGR5+ stem cells, 607 TFs showed significant activity differences between control and AS (padj < 0.05; 255 increased, 352 decreased). Among the most significantly altered were REST (decreased; a master repressor maintaining neural and stem cell quiescence), MYC (decreased; required for stem cell self-renewal and proliferation), and NKX3-1 (decreased). IRF9 was the most significantly upregulated TF, suggesting activation of interferon signaling pathways in AS stem cells.

CD133+ progenitors exhibited 565 significantly altered TFs (161 increased, 404 decreased), with a notable bias toward downregulation. GRHL2 (decreased; essential for epithelial identity and barrier function), RBPJ (decreased; the core transcriptional mediator of canonical Notch signaling), and NFE2L2 (decreased; master regulator of oxidative stress response) were among the most significantly reduced. ZEB2 (increased), a key EMT transcription factor, was the most significantly upregulated, suggesting potential epithelial-mesenchymal transition in progenitor cells.

Basalis fibroblasts showed 591 altered TFs (395 increased, 196 decreased), with a bias toward upregulation. STAT2 (increased; inflammatory signaling), GATA4 (increased; mesenchymal differentiation), and NCOR2 (increased; transcriptional co-repressor) were among the top upregulated TFs, consistent with inflammatory and fibrotic activation.

### Cell-cell communication network remodeling in the stem cell niche

Consensus-based cell-cell communication analysis was performed using liana³² (rank_aggregate method combining CellChat³³, CellPhoneDB³⁴, NATMI³⁵, and other methods; consensus resource with 4,624 ligand-receptor pairs; 1,000 permutations) (Figure 5C-F).

**Global communication landscape**: From 56,760 tested interactions, 1,430 reached significance (magnitude_rank < 0.01). The communication network revealed a densely interconnected tissue architecture with decidualized stromal cells and basalis fibroblasts as the most active signaling hubs.

**Stem cell niche interactions**: 636 significant interactions involved stem/progenitor cell types (315 incoming, 435 outgoing). SOX9+LGR5+ stem cells received signals primarily from basalis fibroblasts (19 interactions), decidualized stromal cells (17), and smooth muscle (12), with autocrine signaling contributing 20 interactions. CD133+ progenitors were predominantly paracrine-dependent (65 paracrine vs. 1 autocrine interaction), with basalis fibroblasts as the primary signal source. Top ligand-receptor pairs included TIMP1→CD63, LGALS1→ITGB1, and LUM→ITGB1, indicating ECM/adhesion-mediated niche maintenance.

**Disease-associated changes**: Comparison of control and AS communication networks identified systematic remodeling. Among stem cell-related interactions, fibrosis-associated signals (TGFB, collagen, FN1, PDGF, MMPs, inflammatory cytokines) showed net enhancement in AS (2,167 gained vs. 1,433 lost), while regeneration-associated signals (WNT, NOTCH, FGF, EGF, BMP) were more balanced (1,233 gained vs. 1,105 lost). Twenty pathway-level communication changes reached significance across stem cell types.

### Multi-evidence drug target prioritization identifies actionable therapeutic targets

To translate our findings into therapeutic opportunities, we developed a multi-evidence scoring framework integrating results from Phases 5–7 of our analysis (Figure 6A-D; Methods).

**Target scoring**: Candidate genes were scored based on five criteria: (1) DEG evidence (number of cell types with significant differential expression, weighted by effect size); (2) stem cell specificity (bonus for DEG in CD133+ progenitors or basalis fibroblasts); (3) TF regulatory evidence (number of cell types with significant TF activity changes); (4) cell communication evidence (identification as altered ligand or receptor); and (5) prior knowledge (known stem cell gene or fibrosis gene). This yielded 463 scored candidate targets (Supplementary Table S7).

**Drug database mining**: DGIdb³⁶ queries (GraphQL API) identified drug interactions for 60 of the top 100 targets, encompassing 3,030 drug-target pairs and 2,289 unique drugs, of which 1,123 were FDA-approved drugs targeting 52 genes. We classified targets into three druggability tiers: Tier 1 (approved drugs available, n=24), Tier 2 (investigational compounds, n=5), and Tier 3 (novel targets without known drugs, n=21).

**Top actionable targets** (Supplementary Table S8):
- **MET** (score = 21; Tier 1): The HGF receptor, identified as a stem cell gene altered in AS with 21 approved drugs (cabozantinib, crizotinib, tepotinib). MET was detected in 50.3% of organoid cells, supporting its candidacy for organoid-based drug screening.
- **CXCL8** (score = 23; Tier 1): A key inflammatory chemokine upregulated in 6 cell types in AS with 43 drug interactions. CXCL8 was identified as both a DEG and an altered communication ligand.
- **PGR** (score = 20; Tier 1): The progesterone receptor, with 20 approved drugs and established relevance to endometrial biology.
- **NNMT** (score = 21; Tier 1): Nicotinamide N-methyltransferase, a metabolic enzyme differentially expressed in both CD133+ progenitors and basalis fibroblasts, with niacin as a known modulator.
- **IER3** (score = 25; Tier 3): Immediate early response 3, the highest-scoring target overall with DEG evidence in 8 cell types and stem cell specificity, though currently lacking approved drugs.

**Organoid validation**: All top 30 targets were validated for expression in AS patient-derived endometrial organoids (GSE216748; 49,609 cells)²². PAEP (98.2% cells expressing), IER3 (82.5%), LAMB3 (75.2%), ELF3 (71.7%), and MET (50.3%) showed the highest expression levels, confirming their relevance in the endometrial epithelial compartment and suitability for organoid-based functional studies.

---

## Discussion

Our integrative analysis of 314,805 endometrial cells reveals that Asherman's syndrome involves a dual pathological mechanism: "stem cell depletion and functional impairment" coupled with "niche fibrotic reprogramming." This model extends previous observations of altered cellular composition in AS²² by demonstrating that the remaining stem cells are functionally compromised through comprehensive signaling pathway silencing, and that the basalis niche itself undergoes a fundamental transformation from regenerative to fibrotic.

**Stem cell hierarchy in the human endometrium.** Our identification of SOX9+LGR5+ stem cells as the most primitive population, with CD133+ progenitors representing a more committed downstream state, provides clarity to the long-debated relationship between these markers. Previous studies have identified SOX9+ and CD133+ populations independently⁶⁻⁹, but their hierarchical relationship has not been established at the single-cell level. Our DPT trajectory analysis, CytoTRACE scoring, and PAGA connectivity all support a model where SOX9+LGR5+ cells give rise to CD133+ progenitors, which in turn differentiate into mature glandular and secretory epithelium. This hierarchy parallels the intestinal crypt model, where LGR5+ stem cells give rise to transit-amplifying progenitors¹¹, though the endometrial hierarchy appears less strictly organized.

**Comprehensive pathway silencing as a disease mechanism.** The most striking finding is the simultaneous downregulation of all eight examined signaling pathways in SOX9+ stem cells. This pattern is unprecedented in stem cell biology literature, where disease-associated stem cell dysfunction typically involves perturbation of one or two pathways (e.g., WNT in aging intestinal stem cells³⁷ or NOTCH in myelodysplastic syndromes³⁸). The breadth of pathway silencing in AS stem cells suggests a fundamental change in the cell's ability to receive extracellular signals, possibly through altered receptor expression, epigenetic silencing of signaling components, or physical isolation from niche signals by fibrotic tissue.

**Niche fibrotic reprogramming.** The paradoxical upregulation of WNT, NOTCH, TGF-β, Hedgehog, and EGF pathways in basalis fibroblasts provides a mechanistic explanation for the niche transformation. In their normal state, basalis fibroblasts serve as the primary signaling hub for stem cell maintenance (250 outgoing interactions in our communication analysis). In AS, these cells appear to redirect their signaling output toward pro-fibrotic programs, potentially driven by TGF-β pathway activation. The concomitant upregulation of STAT2, GATA4, and NCOR2 transcription factors in AS fibroblasts supports inflammatory and fibrotic activation. This model is consistent with the well-established role of TGF-β in driving tissue fibrosis across multiple organs³⁹.

**Transcription factor networks.** The downregulation of MYC and REST in SOX9+ stem cells is particularly significant. MYC is essential for stem cell self-renewal in virtually all stem cell systems, and its loss is associated with stem cell exhaustion⁴⁰. REST (RE1-silencing transcription factor) maintains stem cell quiescence by repressing lineage-specific genes⁴¹, and its downregulation may result in premature differentiation or senescence. The upregulation of ZEB2 in CD133+ progenitors suggests EMT-like processes that may further erode epithelial stem cell identity⁴². The downregulation of RBPJ, the obligate transcriptional mediator of canonical Notch signaling, provides a direct molecular explanation for the observed NOTCH pathway silencing at the receptor/ligand level.

**Drug target implications.** Our multi-evidence prioritization identifies several targets with immediate translational potential. The HGF/MET axis is particularly compelling: MET is a known stem cell gene, it is identified as both a DEG and an altered communication receptor in AS, and 21 approved drugs are available. While most approved MET-targeting drugs are kinase inhibitors developed for oncology, the role of HGF/MET in promoting stem cell self-renewal and tissue repair is well-established across multiple organs⁴³,⁴⁴. HGF supplementation has shown promise in animal models of endometrial repair⁴⁵, and our finding that MET is expressed in 50.3% of organoid cells provides a clear path for preclinical validation. However, since MET inhibitors would be expected to further impair regeneration, the therapeutic strategy would involve HGF supplementation or MET pathway activation rather than inhibition — a critical distinction for future translational work.

CXCL8 (IL-8) represents a distinct therapeutic angle: as a pro-inflammatory chemokine upregulated in AS, anti-CXCL8 strategies could address the inflammatory component of the disease. PGR (progesterone receptor) is already an established target in reproductive medicine, and our data support hormonal approaches as a component of AS therapy. NNMT, a metabolic enzyme involved in nicotinamide metabolism, has recently emerged as a regulator of mesenchymal stem cell fate⁴⁶ and may represent a novel target amenable to metabolic intervention.

**Limitations and future directions.** Several limitations should be noted. First, our analysis is computational; experimental validation of drug targets in organoid models and in vivo systems is essential. Second, the SOX9+LGR5+ stem cell pseudobulk DEG analysis was limited by insufficient sample representation (only 1 control and 1 AS sample had adequate cell counts), preventing robust differential expression analysis for this critical population. Single-cell-level statistical methods (e.g., MAST⁴⁷) could complement our pseudobulk approach. Third, the PRJNA730360 thin endometrium dataset is still being processed and will enable comparison of AS with another regenerative disorder. Fourth, RNA velocity analysis⁴⁸ and full SCENIC-based gene regulatory network inference⁴⁹ were not performed due to the lack of spliced/unspliced count information and computational constraints, respectively; these analyses could further refine the regulatory landscape. Finally, the DGIdb drug matches identify drug-target interactions but do not predict therapeutic efficacy; the distinction between activating and inhibiting drugs for each target requires careful pharmacological consideration in follow-up studies.

In conclusion, our comprehensive endometrial stem cell atlas reveals the molecular underpinnings of regenerative failure in AS and provides a resource for therapeutic development. The identification of actionable targets with approved drugs opens the possibility of drug repurposing strategies that could accelerate clinical translation.

---

## Methods

### Datasets and data acquisition

Seven publicly available datasets were used in this study (Supplementary Table S1):

**GSE215968** (Santamaria et al., Nature Communications 2023)²²: Single-cell RNA-seq of endometrial samples from Asherman's syndrome patients and controls. Three sub-datasets were obtained: (1) AS vs. WOI Control (106,400 cells from 9 AS and 12 WOI control samples); (2) CD133+-enriched sorted cells (69,701 cells from CD133+ FACS sorting); (3) AS pre- vs. post-treatment (123,250 cells; excluded from integration to avoid treatment-related confounding). Data were downloaded as h5ad files from GEO.

**GSE111976** (Wang et al., Nature Medicine 2020)²³: Single-cell transcriptomes of 73,181 endometrial cells from 27 donors across the menstrual cycle. Raw data were provided as an RDS file containing a dgCMatrix sparse matrix, which was exported to Matrix Market format via R and converted to h5ad format in Python. After quality control, 61,503 cells were retained.

**E-MTAB-10287** (Garcia-Alonso et al., Nature Genetics 2021)¹⁰: Single-cell RNA-seq of human endometrium from 5 patients across 11 time points during the menstrual cycle. Data were downloaded from ArrayExpress in 10x Matrix Market format (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz). After quality control, 49,550 cells were retained.

**GSE260658**: Single-cell RNA-seq of human uterus including endometrium, myometrium, and decidua from 11 samples. Data were downloaded from GEO in 10x Matrix Market format. Myometrial samples were excluded; only endometrium (18,000) and decidua (10,000) samples were retained (27,651 cells total after QC).

**GSE287278** (Scientific Data 2025): 10x Visium spatial transcriptomics of endometrial samples from 4 control and 4 RIF (recurrent implantation failure) patients. Count matrices and spatial coordinates (tissue_positions.csv, scalefactors_json.json) were obtained for all 8 samples (10,163 spots total).

**HECA reference** (Garcia-Alonso et al., Nature Genetics 2024)¹⁰: The Human Endometrial Cell Atlas comprising 313,527 single cells and 312,246 single nuclei from 63 women. Used exclusively as a reference for label transfer annotation; not included in the integrated atlas.

**GSE216748** (Santamaria et al., Nature Communications 2023)²²: Single-cell RNA-seq of 49,609 cells from AS patient-derived endometrial organoids treated with hormones. Used exclusively for independent validation of drug target expression.

### Quality control and preprocessing

For previously unprocessed datasets (GSE111976, E-MTAB-10287, GSE260658), quality control was performed using the following criteria: mitochondrial gene percentage < 20%, number of detected genes between 200 and 8,000, total UMI counts > 500, and minimum 3 cells per gene. Doublet detection was performed using Scrublet⁵⁰ (expected_doublet_rate=0.06) independently for each sample to avoid cross-sample artifacts. Doublet rates ranged from 0.2% (GSE111976) to 3.4% (E-MTAB-10287).

For GSE215968 and GSE216748, the original authors' quality control was retained (mitochondrial percentage < 25%, as specified in Santamaria et al.²²), as these datasets were already filtered and annotated.

Data normalization was performed using scanpy⁵¹ (v1.11.5): raw counts were normalized to 10,000 counts per cell (sc.pp.normalize_total), followed by log1p transformation (sc.pp.log1p). Normalized values were stored in the 'log_normalized' layer. For GSE215968, which was provided as Seurat LogNormalize-transformed data, raw counts were recovered by applying expm1 and rescaling by 10,000, with verification that the reconstruction error was < 0.3 per cell.

### Feature selection and dimensionality reduction

Highly variable genes (HVGs) were selected using the Seurat v3 method (sc.pp.highly_variable_genes with flavor='seurat_v3', n_top_genes=3,000) applied to the raw count data with dataset-aware batch information. HVG selection was performed on the concatenated dataset to identify genes variable across conditions.

Principal component analysis was performed on the HVG-subset data after scaling (zero mean, unit variance) using randomized SVD (svd_solver='randomized', n_comps=50) for computational efficiency. The randomized SVD approach reduced PCA computation time from >30 minutes (arpack solver) to ~1 minute for 314,805 cells.

### Batch integration

Harmony²⁴ (harmonypy v0.0.10 with GPU acceleration) was used for batch correction of the 50-dimensional PC space, with 'dataset' as the batch variable. Harmony converged in 3 iterations (~7 seconds on NVIDIA RTX A4000 GPU). The corrected embedding was stored in obsm['X_pca_harmony'].

Integration quality was assessed using: (1) batch-corrected average silhouette width (ASW_batch = 0.024, where values near 0 indicate good batch mixing); (2) cluster mixing score (0.565); (3) mean cluster entropy (0.91); and (4) visual inspection of UMAP embeddings colored by dataset and cell type.

### Clustering

Nearest-neighbor graphs were constructed using scanpy (sc.pp.neighbors, n_neighbors=30, metric='euclidean') on the Harmony-corrected PCs. UMAP²⁷ was computed for visualization (sc.tl.umap, random_state=42). Leiden clustering⁵² was performed at multiple resolutions (0.5, 0.8, 1.0, 1.5) to explore different granularities. Resolution 1.0 (29 clusters) was used as the primary clustering for downstream annotation.

### Cell type annotation

Annotation was performed in two stages:

**Stage 1 — Reference-based label transfer**: KNN-based label transfer from the HECA reference atlas (313,527 cells, 36 cell types)¹⁰ was performed using k=30 nearest neighbors with distance-weighted voting. The reference and query datasets were projected into a shared PC space using 15,874 common genes. For each query cell, the 30 nearest reference cells in this shared space were identified, and cell type labels were assigned by weighted majority voting, where weights were inversely proportional to Euclidean distance. Confidence scores were calculated as the proportion of nearest neighbors sharing the majority label. Mean cell type confidence was 0.866, and mean lineage confidence was 0.994.

**Stage 2 — Marker-guided manual refinement**: Reference-transferred labels were refined using canonical marker gene expression examined across the 29 Leiden clusters. Key annotation decisions included: (1) identification of SOX9+LGR5+ stem cells (cluster 12) based on co-expression of SOX9 (48.8%), LGR5 (25.0%), and PROM1 (26.4%); (2) identification of CD133+ progenitors (clusters 18, 22, 28) based on high PROM1 (47.5%) and low SOX9 (9.5%); (3) annotation of basalis fibroblasts (cluster 24) consistent with HECA Fibroblast_basalis; (4) merging of decidualized stromal subclusters (clusters 1, 2, 6, 14, 15, 17, 19) sharing mesenchymal lineage markers.

### Spatial transcriptomics processing and deconvolution

**Preprocessing**: Visium data from 8 samples were loaded using scanpy's read_visium, with spatial coordinates aligned from tissue_positions.csv files. Quality control filters included: UMI count ≥ 500, detected genes ≥ 200, mitochondrial percentage ≤ 25%, and minimum 5 cells per gene. After QC, 10,135 spots with 14,359 common genes were retained. Batch correction was performed using Harmony on 30 PCs across the 8 samples.

**Region annotation**: Spatial regions were annotated based on composite marker scores for 8 marker categories (stromal, epithelial, basalis, functionalis, endothelial, immune, smooth muscle, glandular markers), computed using sc.tl.score_genes. Unsupervised Leiden clustering (resolution 0.3) identified 7 spatial clusters, which were assigned to three regions based on marker scores: Stroma (clusters 0, 1, 4, 5; 6,480 spots), Functionalis (clusters 2, 3; 3,423 spots), and Basalis_niche (cluster 6; 232 spots).

**Cell2location deconvolution**: Cell2location²⁸ (v0.1.5) was used for reference-based spatial deconvolution. The reference consisted of 60,445 cells subsampled from the integrated atlas (maximum 5,000 cells per cell type) across 14 cell types, with 13,017 genes common to both reference and spatial data. The negative binomial (NB) regression model was trained for 100 epochs (ELBO convergence ~4.54×10⁸ at epoch 50). The spatial model was trained for 5,000 epochs (ELBO convergence ~7.97×10⁷ at epoch 3,000). All computations used GPU acceleration (NVIDIA RTX A4000, 16 GB). Cell type abundances were extracted as the 5th percentile of the posterior distribution (q05cell_abundance_w_sf).

### Differential expression analysis

Pseudobulk DEG analysis was performed using PyDESeq2²⁹ (v0.4). For each cell type, cells from each sample were aggregated by summing raw counts to create pseudobulk profiles. Cell types were analyzed if present in at least 3 samples per condition with ≥ 50 cells each. The DESeq2 model included condition (WOI Control vs. AS) as the primary variable. Genes with padj < 0.05 and |log2FoldChange| > 0.25 were considered significant.

### Cell proportion analysis

Per-sample cell type proportions were calculated by dividing the number of cells of each type by the total cells per sample. Samples with fewer than 50 total cells were excluded. Proportions were compared between WOI Control (n=12 samples) and AS (n=16 samples) using two-sided Wilcoxon rank-sum tests, with p-values adjusted by Benjamini-Hochberg FDR correction.

### Pathway enrichment analysis

GSEA was performed using gseapy⁵³ (v1.1) with preranked gene lists. Genes were ranked by a composite score: −log10(pvalue) × sign(log2FoldChange). Four standard gene set databases were queried: GO Biological Process (MSigDB c5.bp), KEGG (MSigDB c2.cp.kegg), Reactome (MSigDB c2.cp.reactome), and WikiPathway (MSigDB c2.cp.wikipathways). Additionally, six custom gene sets relevant to endometrial stem cell biology were curated from literature: WNT_SIGNALING_STEM, NOTCH_SIGNALING_STEM, TGFB_SIGNALING, FIBROSIS_MARKERS, ENDOMETRIAL_REGENERATION, and STEM_NICHE_FACTORS (each containing 15–30 genes; Supplementary Table S4).

Over-representation analysis (ORA) was performed using gseapy for GO Biological Process terms, with differentially expressed genes (padj < 0.05, |log2FC| > 0.5) as input and all detected genes as background.

### Trajectory analysis

**Diffusion pseudotime**: DPT²⁶ was computed using scanpy (sc.tl.diffmap with n_comps=15, followed by sc.tl.dpt) on subsets of the integrated atlas. For the epithelial lineage, cells annotated as SOX9+LGR5+_stem, CD133+_progenitor, Glandular_epithelium, and Secretory_glandular were extracted (72,417 cells). For the mesenchymal lineage, Basalis_fibroblast, Decidualized_stromal, Smooth_muscle, and SOX9+LGR5+_stem cells were used (134,504 cells). Root cells were selected based on the highest CytoTRACE scores within the expected root population. Nearest-neighbor graphs were computed on Harmony-corrected PCs (n_neighbors=30).

**PAGA**: Partition-based graph abstraction²⁷ (sc.tl.paga) was computed using cell type labels as groups, with the default threshold of 0.05 for connectivity significance.

**CytoTRACE scoring**: Due to the unavailability of the CellRank CytoTRACE kernel (which requires a moment-based 'Ms' layer from RNA velocity), we implemented a simplified CytoTRACE-like scoring²⁵ based on gene detection counts. For each cell, the number of genes with non-zero expression was computed, and scores were normalized to [0, 1] within each cell type to enable cross-type comparisons.

### Transcription factor activity inference

TF activities were inferred using the Univariate Linear Model (ULM) method³⁰ implemented in decoupler (v2.1.4). The TF-target regulatory network was obtained from CollecTRI via OmniPath³¹ (omnipath.interactions.CollecTRI.get, genesymbols=True), yielding 62,411 TF-target regulatory interactions for 760 TFs after filtering complex annotations and deduplication. Regulatory weights were set to +1 for stimulatory interactions and −1 for inhibitory interactions. ULM was applied to the log-normalized expression matrix of 23,123 stem/progenitor cells (SOX9+LGR5+_stem, CD133+_progenitor, Basalis_fibroblast, Glandular_epithelium, Secretory_glandular, Decidualized_stromal). Disease-associated TF activity changes were assessed using Mann-Whitney U tests between control and AS cells within each cell type, with Benjamini-Hochberg FDR correction.

### Cell-cell communication analysis

Cell-cell communication was analyzed using liana³² (v1.7.1) with the rank_aggregate method, which combines scores from multiple methods (CellChat³³, CellPhoneDB³⁴, NATMI³⁵, SingleCellSignalR, Connectome, and logFC) into a consensus ranking. The 'consensus' ligand-receptor resource (4,624 interactions) was used. Minimum expression proportion was set to 0.1 (at least 10% of cells in a group must express the ligand/receptor). Significance was assessed by 1,000 permutations.

For disease comparison, liana was run separately on control and AS subsets, and interaction rankings were compared to identify gained and lost interactions. An interaction was classified as "gained in AS" if its magnitude_rank decreased by > 0.1, and "lost in AS" if it increased by > 0.1.

**Pathway-level communication scoring**: For each of eight signaling pathways (WNT, NOTCH, TGF-β/BMP, FGF, EGF, Hedgehog, Cytokine, Chemokine), curated gene lists of ligands and receptors were compiled from literature (15–40 genes per pathway). Mean expression of pathway genes was computed per cell and compared between control and AS using Mann-Whitney U tests with FDR correction.

For downsampled analyses, a maximum of 5,000 cells per cell type were randomly sampled to ensure computational tractability while maintaining statistical power.

### Drug target prioritization

A multi-evidence scoring framework was developed to prioritize drug targets:

**Scoring criteria**: Each candidate gene received points for: (1) **DEG evidence**: 2 points per cell type with significant differential expression + up to 5 points for maximum |log2FC|; (2) **Stem cell specificity**: 3 bonus points per stem/progenitor cell type (CD133+_progenitor, Basalis_fibroblast) with DEG evidence; (3) **TF regulation**: 2 points per cell type with significant TF activity change; (4) **Communication evidence**: 3 points each for identification as an altered ligand or receptor; (5) **Prior knowledge**: 2 points for known stem cell genes (39-gene curated set) or fibrosis genes (22-gene curated set).

**Drug database query**: The DGIdb³⁶ GraphQL API (https://dgidb.org/api/graphql) was queried in batches of 25 genes for the top 100 targets. Returned drug interactions were classified by approval status and interaction type.

**Druggability tiering**: Tier 1 — at least one FDA-approved drug available (n=24 of top 50); Tier 2 — investigational compounds only (n=5); Tier 3 — no known drug interactions, novel target (n=21).

**Organoid validation**: Drug target expression was validated in an independent AS patient-derived organoid dataset (GSE216748; 49,609 cells). For each target gene, mean expression and percentage of expressing cells were computed from the raw count matrix.

### Statistical analysis

All statistical tests were performed in Python (scipy v1.14, statsmodels v0.14). Two-group comparisons used Mann-Whitney U tests (continuous data) or Wilcoxon rank-sum tests (proportional data). Multiple testing correction used the Benjamini-Hochberg method. Distribution comparisons used two-sided Kolmogorov-Smirnov tests. Significance thresholds: * p < 0.05, ** p < 0.01, *** p < 0.001.

### Software and computational environment

All analyses were performed in a conda environment (endometrium_atlas) with Python 3.11.14 and R 4.4.3 on a Linux server (64 CPU cores, 235 GB RAM, NVIDIA RTX A4000 16 GB GPU, WSL2). Key Python packages: scanpy 1.11.5, anndata 0.12.10, harmonypy 0.0.10, cell2location 0.1.5, liana 1.7.1, decoupler 2.1.4, pydeseq2 0.4, gseapy 1.1, omnipath 1.0.8, scrublet 0.2.3, matplotlib 3.10, seaborn 0.13. Key R packages: Seurat 5.4.0, CellChat 2.2.0, nichenetr 2.2.1, monocle3 1.4.26.

### Data and code availability

All raw sequencing data and processed count matrices are available from GEO (GSE215968, GSE216748, GSE111976, GSE260658, GSE287278, GSE234354, GSE127918), ArrayExpress (E-MTAB-9260, E-MTAB-10287), and the HECA portal. Analysis scripts are available at [repository URL to be provided upon publication].

---

## Acknowledgments
We thank the original authors of all datasets used in this study for making their data publicly available.

## Author Contributions
[To be completed]

## Competing Interests
The authors declare no competing interests.

## Figure Legends

**Figure 1. Construction of a high-resolution human endometrial cell atlas.** (A) Study design and analysis workflow showing integration of five scRNA-seq datasets comprising 314,805 cells. (B) UMAP embedding colored by cell type annotation (14 types) and by dataset origin (5 datasets). (C) Cell type proportions across datasets shown as stacked bar plots. (D) Lineage distribution pie chart showing mesenchymal (39.9%), immune (31.0%), epithelial (16.3%), stem/progenitor (6.7%), and endothelial (4.6%) compartments. (E) Cell counts per cell type.

**Figure 2. Identification of two distinct stem/progenitor populations.** (A-D) Feature plots showing expression of SOX9, LGR5, PROM1, and CD44 on the UMAP embedding. (E) UMAP highlighting the three stem/progenitor populations (SOX9+LGR5+_stem, CD133+_progenitor, Basalis_fibroblast) against a grey background of other cell types. (F) Mean expression heatmap of stem cell markers (SOX9, LGR5, PROM1, CD44) across stem/progenitor and selected differentiated cell types.

**Figure 3. Spatial validation and differentiation trajectories.** (A) Visium spatial regions annotated as Stroma, Functionalis, and Basalis_niche. (B-C) Cell2location deconvolution showing spatial distribution of SOX9+LGR5+ stem cells and CD133+ progenitors. (D) PAGA connectivity graph showing epithelial lineage relationships. (E) Diffusion pseudotime on UMAP for the epithelial lineage (72,417 cells). (F) Pseudotime distribution violin plots by cell type.

**Figure 4. Stem cell depletion and pathway silencing in Asherman's syndrome.** (A) Cell type proportion changes between WOI Control and AS (Wilcoxon test, * padj < 0.05). (B) CytoTRACE stemness scores in Control vs. AS for SOX9+LGR5+_stem, CD133+_progenitor, and Basalis_fibroblast (*** p < 0.001). (C) Number of differentially expressed genes (up/down) per cell type in AS vs. Control. (D) Signaling pathway expression changes in stem cell types (AS vs. Control). Heatmap shows mean expression difference for 8 pathways across 3 stem/progenitor types. Annotations: * padj < 0.05, ** padj < 0.01, *** padj < 0.001.

**Figure 5. Transcription factor dysregulation and cell communication remodeling.** (A) Top variable TF activities (ULM) across cell types shown as heatmap. (B) Disease-associated TF activity changes in stem/progenitor cells (AS vs. Control), with significance annotations (* padj < 0.05). (C) Cell-cell interaction count heatmap (sender × receiver) for significant interactions (magnitude_rank < 0.01). (D) Stem cell niche pathway communication summary showing number of significant interactions by pathway and direction (incoming/outgoing). (E) Signaling pathway communication changes in AS shown as heatmap of mean rank differences. (F) Distribution of rank differences for fibrosis-related vs. regeneration-related interactions in the stem cell niche.

**Figure 6. Drug target prioritization and validation.** (A) Top 25 drug target candidates ranked by multi-evidence priority score, colored by druggability tier (green: Tier 1/approved drugs; orange: Tier 2/investigational; red: Tier 3/novel). (B) Multi-evidence support heatmap showing DEG cell type count, drug availability, communication evidence, and overall score for top 15 targets. (C) Druggability assessment pie chart for top 25 targets. (D) Number of approved drugs available for top druggable targets.

---

## Supplementary Information

### Supplementary Tables
- **Table S1**: Dataset metadata including source, accession, cell counts, QC parameters, and processing details.
- **Table S2**: Cell type marker genes used for manual annotation refinement, with expression statistics per cluster.
- **Table S3**: Full differential expression results for all 10 cell types (AS vs. WOI Control), including baseMean, log2FoldChange, lfcSE, stat, pvalue, and padj.
- **Table S4**: GSEA and ORA pathway enrichment results for all gene set databases and custom pathways.
- **Table S5**: TF activity differential results for 5 cell types, including TF name, activity difference, p-value, and adjusted p-value.
- **Table S6**: Complete cell-cell communication results from liana (global, control, AS, and differential analyses).
- **Table S7**: Full drug target prioritization list (463 candidates) with scores, evidence types, and druggability classification.
- **Table S8**: DGIdb drug-gene interaction results for top 100 targets, including drug name, approval status, interaction type, and publications.

### Supplementary Figures
- **Figure S1**: Quality control metrics distribution across datasets (UMI counts, gene counts, mitochondrial percentage, doublet scores).
- **Figure S2**: Batch integration quality assessment (UMAP before/after Harmony, ASW scores, cluster mixing).
- **Figure S3**: HECA reference mapping confidence distribution and lineage assignment accuracy.
- **Figure S4**: Complete Cell2location deconvolution results for all 14 cell types in Visium spatial data.
- **Figure S5**: Mesenchymal lineage trajectory analysis (PAGA connectivity, DPT distribution, gene expression along trajectory).
- **Figure S6**: Complete DEG results (volcano plots for all 10 cell types).
- **Figure S7**: Complete GSEA results (enrichment plots for top pathways per cell type).
- **Figure S8**: TF activity UMAP visualization for key TFs (SOX9, TP53, MYC, JUN, STAT3, FOXO1, HIF1A, NFKB1).
- **Figure S9**: Complete cell-cell communication network including all cell types and interaction types.
- **Figure S10**: Organoid expression validation for top 30 drug targets.

---

## References

1. Jabbour HN, Kelly RW, Fraser HM, Critchley HOD. Endocrine regulation of menstruation. Endocr Rev. 2006;27(1):17-46.
2. Gargett CE, Schwab KE, Deane JA. Endometrial stem/progenitor cells: the first 10 years. Hum Reprod Update. 2016;22(2):137-163.
3. Cousins FL, Pandoy R, Jin S, Gargett CE. The elusive endometrial epithelial stem/progenitor cells. Front Cell Dev Biol. 2021;9:640319.
4. Chan RW, Schwab KE, Gargett CE. Clonogenicity of human endometrial epithelial and stromal cells. Biol Reprod. 2004;70(6):1738-1750.
5. Schwab KE, Gargett CE. Co-expression of two perivascular cell markers isolates mesenchymal stem-like cells from human endometrium. Hum Reprod. 2007;22(11):2903-2911.
6. Masuda H, Matsuzaki Y, Hiratsu E, et al. Stem cell-like properties of the endometrial side population: implication in endometrial regeneration. PLoS One. 2010;5(4):e10387.
7. Cervello I, Gil-Sanchis C, Mas A, et al. Human endometrial side population cells exhibit genotypic, phenotypic and functional features of somatic stem cells. PLoS One. 2010;5(6):e10964.
8. Valentijn AJ, Palial K, Al-lamee H, et al. SSEA-1 isolates human endometrial basal glandular epithelial cells: phenotypic and functional characterization and implications in the pathogenesis of endometriosis. Hum Reprod. 2013;28(10):2695-2708.
9. Hapangama DK, Drury J, Da Silva L, et al. Abnormally located SSEA1+/SOX9+ endometrial epithelial cells with a basalis-like phenotype in the eutopic functionalis layer may play a role in the pathogenesis of endometriosis. Hum Reprod. 2019;34(1):56-68.
10. Garcia-Alonso L, Handfield LF, Roberts K, et al. Mapping the temporal and spatial dynamics of the human endometrium in vivo and in vitro. Nat Genet. 2021;53(12):1698-1711.
11. Gil-Sanchis C, Cervello I, Khurana S, et al. Contribution of different bone marrow-derived cell types in endometrial regeneration using an irradiated murine model. Fertil Steril. 2015;103(6):1596-1605.
12. Tempest N, Baker AM, Wright NA, Hapangama DK. Does human endometrial LGR5 gene expression suggest the existence of another stem cell niche? Hum Reprod. 2018;33(1):60-65.
13. Barros FSV, Brosens JJ, Brighton PJ. Isolation and primary culture of various cell types from whole human endometrial biopsies. Bio Protoc. 2016;6(22):e2028.
14. Gargett CE, Ye L. Endometrial reconstruction from stem cells. Fertil Steril. 2012;98(1):11-20.
15. Syed SM, Kumar M, Ghosh A, et al. Endometrial Axin2+ cells drive epithelial homeostasis, regeneration, and cancer following oncogenic transformation. Cell Stem Cell. 2020;26(1):64-80.
16. Boretto M, Cox B, Noben M, et al. Development of organoids from mouse and human endometrium showing endometrial epithelium physiology and long-term expandability. Development. 2017;144(10):1775-1786.
17. March CM. Management of Asherman's syndrome. Reprod Biomed Online. 2011;23(1):63-76.
18. Yu D, Wong YM, Cheong Y, et al. Asherman syndrome—one century later. Fertil Steril. 2008;89(4):759-779.
19. Dreisler E, Kjer JJ. Asherman's syndrome: current perspectives on diagnosis and management. Int J Womens Health. 2019;11:191-198.
20. Hooker AB, Lemmers M, Thurkow AL, et al. Systematic review and meta-analysis of intrauterine adhesions after miscarriage: prevalence, risk factors, and long-term reproductive outcome. Hum Reprod Update. 2014;20(2):262-278.
21. Johary J, Xue M, Zhu X, Xu D, Velu PP. Efficacy of estrogen therapy in patients with intrauterine adhesions: systematic review. J Minim Invasive Gynecol. 2014;21(1):44-54.
22. Santamaria X, Roson B, Perez-Moraga R, et al. Decoding the endometrial niche of Asherman's syndrome at single-cell resolution. Nat Commun. 2023;14:5135.
23. Wang W, Vilella F, Alama P, et al. Single-cell transcriptomic atlas of the human endometrium during the menstrual cycle. Nat Med. 2020;26(10):1644-1653.
24. Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods. 2019;16(12):1289-1296.
25. Gulati GS, Sikandar SS, Wesche DJ, et al. Single-cell transcriptional diversity is a hallmark of developmental potential. Science. 2020;367(6476):405-411.
26. Haghverdi L, Büttner M, Wolf FA, Buettner F, Theis FJ. Diffusion pseudotime robustly reconstructs lineage branching. Nat Methods. 2016;13(10):845-848.
27. Wolf FA, Hamey FK, Plass M, et al. PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells. Genome Biol. 2019;20(1):59.
28. Kleshchevnikov V, Shmatko A, Dann E, et al. Cell2location maps fine-grained cell types in spatial transcriptomics. Nat Biotechnol. 2022;40(5):661-671.
29. Squair JW, Gautier M, Kathe C, et al. Confronting false discoveries in single-cell differential expression. Nat Commun. 2021;12(1):5692.
30. Badia-i-Mompel P, Vélez Santiago J, Braunger J, et al. decoupler: ensemble of computational methods to infer biological activities from omics data. Bioinformatics Advances. 2022;2(1):vbac016.
31. Türei D, Korcsmáros T, Saez-Rodriguez J. OmniPath: guidelines and gateway for literature-curated signaling pathway resources. Nat Methods. 2016;13(12):966-967.
32. Dimitrov D, Türei D, Garrido-Rodriguez M, et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun. 2022;13(1):3224.
33. Jin S, Guerrero-Juarez CF, Zhang L, et al. Inference and analysis of cell-cell communication using CellChat. Nat Commun. 2021;12(1):1088.
34. Efremova M, Vento-Tormo M, Teichmann SA, Vento-Tormo R. CellPhoneDB: inferring cell-cell communication from combined expression of multi-subunit ligand-receptor complexes. Nat Protoc. 2020;15(4):1484-1506.
35. Hou R, Denisenko E, Ong HT, Ramilowski JA, Forrest ARR. Predicting cell-to-cell communication networks using NATMI. Nat Commun. 2020;11(1):5011.
36. Freshour SL, Kiwala S, Cotto KC, et al. Integration of the Drug-Gene Interaction Database (DGIdb 4.0) with open crowdsource efforts. Nucleic Acids Res. 2021;49(D1):D1144-D1151.
37. Nalapareddy K, Nattamai KJ, Kumar RS, et al. Canonical Wnt signaling ameliorates aging of intestinal stem cells. Cell Rep. 2017;18(11):2608-2621.
38. Kode A, Manavalan JS, Mosialou I, et al. Leukaemogenesis induced by an activating β-catenin mutation in osteoblasts. Nature. 2014;506(7487):240-244.
39. Henderson NC, Rieder F, Wynn TA. Fibrosis: from mechanisms to medicines. Nature. 2020;587(7835):555-566.
40. Laurenti E, Varnum-Finney B, Wilson A, et al. Hematopoietic stem cell function and survival depend on c-Myc and N-Myc activity. Cell Stem Cell. 2008;3(6):611-624.
41. Gao Z, Ure K, Ding P, et al. The master negative regulator REST/NRSF controls adult neurogenesis by restraining the neurogenic program in quiescent stem cells. J Neurosci. 2011;31(26):9772-9786.
42. Nieto MA, Huang RY, Jackson RA, Thiery JP. EMT: 2016. Cell. 2016;166(1):21-45.
43. Matsumoto K, Nakamura T. Hepatocyte growth factor and the Met system as a mediator of tumor-stromal interactions. Int J Cancer. 2006;119(3):477-483.
44. Nakamura T, Sakai K, Nakamura T, Matsumoto K. Hepatocyte growth factor twenty years on: Much more than a growth factor. J Gastroenterol Hepatol. 2011;26 Suppl 1:188-202.
45. Aghajanova L, Cedars MI, Huber A, Irwin JC, Giudice LC. Expression of hepatocyte growth factor receptor (c-met) and ligand (HGF) in human endometrium. Fertil Steril. 2009;91(1):15-23.
46. Eckert MA, Coscia F, Chryplewicz A, et al. Proteomics reveals NNMT as a master metabolic regulator of cancer-associated fibroblasts. Nature. 2019;569(7758):723-728.
47. Finak G, McDavid A, Yajima M, et al. MAST: a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data. Genome Biol. 2015;16:278.
48. La Manno G, Soldatov R, Zeisel A, et al. RNA velocity of single cells. Nature. 2018;560(7719):494-498.
49. Aibar S, González-Blas CB, Moerman T, et al. SCENIC: single-cell regulatory network inference and clustering. Nat Methods. 2017;14(11):1083-1086.
50. Wolock SL, Lopez R, Klein AM. Scrublet: Computational identification of cell doublets in single-cell transcriptomic data. Cell Syst. 2019;8(4):281-291.
51. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol. 2018;19(1):15.
52. Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. Sci Rep. 2019;9(1):5233.
53. Fang Z, Liu X, Peltz G. GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. Bioinformatics. 2023;39(1):btac757.
