# ============================================================
# Phase 1.6: Fibroblast subpopulation analysis
# Phase 1.7: Macrophage analysis
# Phase 1.8: Cell communication analysis
# ============================================================

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

cat("=== Phase 1.6: Fibroblast subpopulation analysis ===\n")

outdir_fib <- "/home/moog/test/skin/Phase_output/Phase1/fibroblast"
outdir_mac <- "/home/moog/test/skin/Phase_output/Phase1/macrophage"

seu <- readRDS("/home/moog/test/skin/Phase_output/Phase1/clustering/seurat_annotated.rds")

# ============ 1.6 Fibroblast subpopulations ============

fib <- subset(seu, celltype == "Fibroblast")
cat(sprintf("Fibroblast cell count: %d\n", ncol(fib)))

# Re-cluster
fib <- NormalizeData(fib)
fib <- FindVariableFeatures(fib, nfeatures = 2000)
fib <- ScaleData(fib, vars.to.regress = "percent.mt")
fib <- RunPCA(fib, npcs = 30, verbose = FALSE)

# Harmony
library(harmony)
fib <- RunHarmony(fib, group.by.vars = "sample", dims.use = 1:20)
fib <- RunUMAP(fib, reduction = "harmony", dims = 1:20)
fib <- FindNeighbors(fib, reduction = "harmony", dims = 1:20)
fib <- FindClusters(fib, resolution = 0.6)

cat(sprintf("Number of subpopulations: %d\n", length(levels(Idents(fib)))))

# HE-Fibro marker genes
he_fibro_markers <- c("MMP1", "MMP3", "CHI3L1", "TNFAIP6", "IL11", "MMP13",
                      "HIF1A", "NFKB1", "CXCL6", "CXCL8")
existing_he <- he_fibro_markers[he_fibro_markers %in% rownames(fib)]

# DotPlot
p_fib_dot <- DotPlot(fib, features = existing_he) +
  coord_flip() + ggtitle("HE-Fibro Markers in Fibroblast Subclusters")
ggsave(file.path(outdir_fib, "HE_fibro_markers_dotplot.png"), p_fib_dot,
       width = 12, height = 6, dpi = 150)

# Identify HE-Fibro subpopulation (high MMP1/MMP3/HIF1A expression)
# Calculate HE-Fibro module score
fib <- AddModuleScore(fib, features = list(existing_he), name = "HEFibro_score")

# Find cluster with highest score
he_scores <- FetchData(fib, vars = c("HEFibro_score1", "seurat_clusters")) %>%
  group_by(seurat_clusters) %>%
  summarise(mean_score = mean(HEFibro_score1), .groups = "drop") %>%
  arrange(desc(mean_score))

cat("\nHE-Fibro score ranking:\n")
print(as.data.frame(he_scores))

# Label HE-Fibro (take the cluster with highest score)
he_cluster <- as.character(he_scores$seurat_clusters[1])
cat(sprintf("HE-Fibro cluster: %s\n", he_cluster))

fib@meta.data$fib_subtype <- ifelse(fib@meta.data$seurat_clusters == he_cluster,
                                     "HE-Fibro", "Other-Fibro")

# Proportion of HE-Fibro in each condition
he_prop <- fib@meta.data %>%
  group_by(condition) %>%
  summarise(
    total = n(),
    HE_fibro = sum(fib_subtype == "HE-Fibro"),
    HE_pct = HE_fibro / total * 100,
    .groups = "drop"
  )
cat("\nHE-Fibro proportion:\n")
print(as.data.frame(he_prop))
write.csv(he_prop, file.path(outdir_fib, "HE_fibro_proportion.csv"), row.names = FALSE)

# UMAP
p_fib1 <- DimPlot(fib, reduction = "umap", label = TRUE, pt.size = 0.2) +
  ggtitle("Fibroblast Subclusters") + NoLegend()
p_fib2 <- DimPlot(fib, reduction = "umap", group.by = "fib_subtype", pt.size = 0.2,
                  cols = c("HE-Fibro" = "#E41A1C", "Other-Fibro" = "#999999")) +
  ggtitle("HE-Fibro Identification")
p_fib3 <- DimPlot(fib, reduction = "umap", group.by = "condition", pt.size = 0.2,
                  cols = c("Healthy"="#4DAF4A","DM"="#377EB8",
                           "Healing"="#FF7F00","NonHealing"="#E41A1C")) +
  ggtitle("Condition")
ggsave(file.path(outdir_fib, "fibroblast_UMAP.png"), p_fib1 + p_fib2 + p_fib3,
       width = 20, height = 6, dpi = 150)

# HIF1A FeaturePlot
p_hif <- FeaturePlot(fib, features = c("HIF1A", "MMP1", "MMP3", "CHI3L1"),
                     reduction = "umap", pt.size = 0.2, ncol = 2)
ggsave(file.path(outdir_fib, "HE_fibro_featureplot.png"), p_hif,
       width = 12, height = 10, dpi = 150)

saveRDS(fib, file.path(outdir_fib, "fibroblast_subclustered.rds"))

# ============ 1.7 Macrophage analysis ============
cat("\n=== Phase 1.7: Macrophage analysis ===\n")

# Merge Macrophage + Myeloid
mac <- subset(seu, celltype %in% c("Macrophage", "Myeloid"))
cat(sprintf("Macrophage+Myeloid cell count: %d\n", ncol(mac)))

mac <- NormalizeData(mac)
mac <- FindVariableFeatures(mac, nfeatures = 2000)
mac <- ScaleData(mac, vars.to.regress = "percent.mt")
mac <- RunPCA(mac, npcs = 30, verbose = FALSE)
mac <- RunHarmony(mac, group.by.vars = "sample", dims.use = 1:20)
mac <- RunUMAP(mac, reduction = "harmony", dims = 1:20)
mac <- FindNeighbors(mac, reduction = "harmony", dims = 1:20)
mac <- FindClusters(mac, resolution = 0.6)

cat(sprintf("Number of subpopulations: %d\n", length(levels(Idents(mac)))))

# M1/M2 marker genes
m1_markers <- c("CD86", "NOS2", "IL1B", "TNF", "IL6", "IRF5", "STAT1", "PTGS2")
m2_markers <- c("CD163", "MRC1", "ARG1", "IL10", "TGFB1", "FOLR2", "CD209", "STAB1")
hif_nfkb <- c("HIF1A", "NFKB1", "RELA", "EPAS1")

all_mac_markers <- c(m1_markers, m2_markers, hif_nfkb)
existing_mac <- all_mac_markers[all_mac_markers %in% rownames(mac)]

p_mac_dot <- DotPlot(mac, features = existing_mac) +
  coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("M1/M2 + HIF/NF-kB Markers")
ggsave(file.path(outdir_mac, "macrophage_markers_dotplot.png"), p_mac_dot,
       width = 12, height = 8, dpi = 150)

# M1/M2 scores
m1_use <- m1_markers[m1_markers %in% rownames(mac)]
m2_use <- m2_markers[m2_markers %in% rownames(mac)]
hif_use <- c("HIF1A", "VEGFA", "SLC2A1", "LDHA", "PGK1", "ENO1", "PDK1",
             "BNIP3", "CA9", "ADM", "ANGPTL4")
hif_use <- hif_use[hif_use %in% rownames(mac)]
nfkb_use <- c("NFKB1", "RELA", "TNF", "IL1B", "IL6", "CXCL8", "NFKBIA", "PTGS2")
nfkb_use <- nfkb_use[nfkb_use %in% rownames(mac)]

mac <- AddModuleScore(mac, features = list(m1_use), name = "M1_score")
mac <- AddModuleScore(mac, features = list(m2_use), name = "M2_score")
mac <- AddModuleScore(mac, features = list(hif_use), name = "HIF_score")
mac <- AddModuleScore(mac, features = list(nfkb_use), name = "NFKB_score")

# M1-M2 differential score
mac@meta.data$M1_M2_diff <- mac@meta.data$M1_score1 - mac@meta.data$M2_score1

# M1/M2 scores per condition
mac_polar <- mac@meta.data %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    M1_mean = mean(M1_score1),
    M2_mean = mean(M2_score1),
    M1_M2_ratio = mean(M1_score1) / mean(abs(M2_score1)),
    HIF_mean = mean(HIF_score1),
    NFKB_mean = mean(NFKB_score1),
    .groups = "drop"
  )
cat("\nMacrophage polarization status:\n")
print(as.data.frame(mac_polar))
write.csv(mac_polar, file.path(outdir_mac, "macrophage_polarization.csv"), row.names = FALSE)

# Correlation between HIF/NF-kB and M1/M2
mac_data <- FetchData(mac, vars = c("HIF_score1", "NFKB_score1", "M1_score1", "M2_score1",
                                    "condition", "seurat_clusters"))

# Correlation
cat("\nHIF/NF-kB and M1/M2 correlation (all cells):\n")
cat(sprintf("  HIF vs M1: rho=%.3f\n", cor(mac_data$HIF_score1, mac_data$M1_score1, method="spearman")))
cat(sprintf("  HIF vs M2: rho=%.3f\n", cor(mac_data$HIF_score1, mac_data$M2_score1, method="spearman")))
cat(sprintf("  NFKB vs M1: rho=%.3f\n", cor(mac_data$NFKB_score1, mac_data$M1_score1, method="spearman")))
cat(sprintf("  NFKB vs M2: rho=%.3f\n", cor(mac_data$NFKB_score1, mac_data$M2_score1, method="spearman")))
cat(sprintf("  HIF vs NFKB: rho=%.3f\n", cor(mac_data$HIF_score1, mac_data$NFKB_score1, method="spearman")))

# UMAP visualization
p_mac1 <- DimPlot(mac, reduction = "umap", label = TRUE, pt.size = 0.3) +
  ggtitle("Macrophage Subclusters") + NoLegend()
p_mac2 <- DimPlot(mac, reduction = "umap", group.by = "condition", pt.size = 0.3,
                  cols = c("Healthy"="#4DAF4A","DM"="#377EB8",
                           "Healing"="#FF7F00","NonHealing"="#E41A1C")) +
  ggtitle("Condition")

p_mac3 <- FeaturePlot(mac, features = c("M1_score1", "M2_score1", "HIF_score1", "NFKB_score1"),
                      reduction = "umap", pt.size = 0.3, ncol = 2) &
  scale_color_viridis_c()
ggsave(file.path(outdir_mac, "macrophage_UMAP.png"), p_mac1 + p_mac2, width = 16, height = 6, dpi = 150)
ggsave(file.path(outdir_mac, "macrophage_scores_UMAP.png"), p_mac3, width = 12, height = 10, dpi = 150)

# M1 vs M2 scatter plot (colored by condition)
p_m1m2 <- ggplot(mac_data, aes(x = M1_score1, y = M2_score1, color = condition)) +
  geom_point(size = 0.3, alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Healthy"="#4DAF4A","DM"="#377EB8",
                                "Healing"="#FF7F00","NonHealing"="#E41A1C")) +
  theme_minimal() + ggtitle("M1 vs M2 Polarization by Condition")
ggsave(file.path(outdir_mac, "M1_vs_M2_scatter.png"), p_m1m2, width = 8, height = 6, dpi = 150)

# Identify HIF+NF-kB high-expression subpopulations
mac_cluster_scores <- mac_data %>%
  group_by(seurat_clusters) %>%
  summarise(
    HIF = mean(HIF_score1), NFKB = mean(NFKB_score1),
    M1 = mean(M1_score1), M2 = mean(M2_score1), n = n(),
    .groups = "drop"
  )
cat("\nMacrophage subcluster scores:\n")
print(as.data.frame(mac_cluster_scores))
write.csv(mac_cluster_scores, file.path(outdir_mac, "macrophage_subcluster_scores.csv"),
          row.names = FALSE)

saveRDS(mac, file.path(outdir_mac, "macrophage_subclustered.rds"))

# ============ 1.8 Cell communication (simplified version) ============
cat("\n=== Phase 1.8: Cell communication analysis (ligand-receptor analysis) ===\n")

# Use simplified expression-based analysis instead of full CellChat
# Focus on HIF/NF-kB related ligand-receptor pairs
lr_pairs <- data.frame(
  ligand = c("VEGFA","TNF","IL1B","IL6","CXCL8","CCL2","TGFB1","HGF","FGF2","PDGFA"),
  receptor = c("KDR","TNFRSF1A","IL1R1","IL6R","CXCR1","CCR2","TGFBR1","MET","FGFR1","PDGFRA"),
  pathway = c("VEGF","TNF","IL1","IL6","CXCL","CCL","TGFb","HGF","FGF","PDGF")
)

# Calculate mean ligand/receptor expression per cell type
ct_levels <- unique(seu$celltype)
lr_results <- list()

for (i in 1:nrow(lr_pairs)) {
  lig <- lr_pairs$ligand[i]
  rec <- lr_pairs$receptor[i]
  if (!(lig %in% rownames(seu) & rec %in% rownames(seu))) next

  lig_expr <- FetchData(seu, vars = c(lig, "celltype", "condition"))
  rec_expr <- FetchData(seu, vars = c(rec, "celltype", "condition"))

  for (cond in c("Healing", "NonHealing")) {
    lig_sub <- lig_expr[lig_expr$condition == cond, ]
    rec_sub <- rec_expr[rec_expr$condition == cond, ]

    lig_mean <- tapply(lig_sub[[lig]], lig_sub$celltype, mean)
    rec_mean <- tapply(rec_sub[[rec]], rec_sub$celltype, mean)

    for (ct_l in names(lig_mean)) {
      for (ct_r in names(rec_mean)) {
        score <- lig_mean[ct_l] * rec_mean[ct_r]
        if (score > 0.01) {
          lr_results[[length(lr_results) + 1]] <- data.frame(
            ligand = lig, receptor = rec, pathway = lr_pairs$pathway[i],
            sender = ct_l, receiver = ct_r, condition = cond,
            lig_expr = lig_mean[ct_l], rec_expr = rec_mean[ct_r],
            lr_score = score
          )
        }
      }
    }
  }
}

lr_df <- do.call(rbind, lr_results)
write.csv(lr_df, "/home/moog/test/skin/Phase_output/Phase1/cellchat/lr_analysis.csv",
          row.names = FALSE)

# Differential communication: Healing vs NonHealing
lr_diff <- lr_df %>%
  group_by(ligand, receptor, pathway, sender, receiver) %>%
  summarise(
    score_healing = sum(lr_score[condition == "Healing"]),
    score_nonhealing = sum(lr_score[condition == "NonHealing"]),
    .groups = "drop"
  ) %>%
  mutate(delta = score_healing - score_nonhealing) %>%
  arrange(desc(abs(delta)))

write.csv(lr_diff, "/home/moog/test/skin/Phase_output/Phase1/cellchat/lr_diff_H_vs_NH.csv",
          row.names = FALSE)

cat("\nTop differential ligand-receptor interactions (Healing vs NonHealing):\n")
print(head(as.data.frame(lr_diff[, c("pathway", "sender", "receiver", "delta")]), 15))

cat("\n=== Phase 1.6-1.8 all complete ===\n")
