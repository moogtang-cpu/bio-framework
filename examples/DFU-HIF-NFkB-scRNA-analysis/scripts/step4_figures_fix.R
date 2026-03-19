# ============================================================
# Step 4 Fix: Generate missing Figure 3, 5, 7 + supplementary figures
# ============================================================

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(pheatmap)
  library(tidyr)
})

cat("=== Step 4 Fix: Generating missing figures ===\n")

figdir <- "/home/moog/test/skin/Phase_output/publication_figures"
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

ct_colors <- c(
  "Fibroblast"="#E41A1C","Keratinocyte"="#377EB8","Endothelial"="#4DAF4A",
  "Macrophage"="#984EA3","T_cell"="#FF7F00","Myeloid"="#A65628",
  "Pericyte"="#F781BF","Smooth_Muscle"="#999999","Mast_cell"="#66C2A5",
  "Melanocyte"="#FC8D62","Lymphatic_Endo"="#8DA0CB"
)
cond_colors <- c("Healthy"="#4DAF4A","DM"="#377EB8","Healing"="#FF7F00","NonHealing"="#E41A1C")

# ========== Figure 3: HE-Fibroblast ==========
cat("Figure 3: HE-Fibroblast...\n")

fib <- readRDS("/home/moog/test/skin/Phase_output/Phase1/fibroblast/fibroblast_subclustered.rds")

# 3a: Fibroblast UMAP, highlight HE-Fibro
fib_meta <- fib@meta.data
if (!"HE_fibro" %in% colnames(fib_meta)) {
  # Recalculate HE-Fibro markers
  he_genes <- c("MMP1","MMP3","CHI3L1","HIF1A","IL6","CXCL8")
  he_exist <- he_genes[he_genes %in% rownames(fib)]
  fib <- AddModuleScore(fib, features = list(he_exist), name = "HE_fibro_score")
  fib@meta.data$HE_fibro <- ifelse(fib@meta.data$HE_fibro_score1 > 0.5, "HE-Fibro", "Other")
}

p3a <- DimPlot(fib, reduction = "umap", group.by = "HE_fibro", pt.size = 0.1,
               cols = c("HE-Fibro" = "#E41A1C", "Other" = "grey85")) +
  ggtitle("Fibroblast Subclusters") +
  theme(plot.title = element_text(size = 12, face = "bold"))

# 3b: FeaturePlot key genes
he_markers <- c("MMP1", "MMP3", "CHI3L1", "HIF1A")
he_markers_exist <- he_markers[he_markers %in% rownames(fib)]
p3b_list <- lapply(he_markers_exist, function(g) {
  FeaturePlot(fib, features = g, pt.size = 0.05, reduction = "umap") +
    scale_color_viridis_c() +
    theme(plot.title = element_text(size = 9),
          legend.key.size = unit(0.3, "cm"))
})
p3b <- wrap_plots(p3b_list, ncol = 2)

# 3c: HE-Fibro proportion bar chart
he_prop <- read.csv("/home/moog/test/skin/Phase_output/Phase1/fibroblast/HE_fibro_proportion.csv")
he_prop$condition <- factor(he_prop$condition, levels = c("Healthy","DM","Healing","NonHealing"))

p3c <- ggplot(he_prop, aes(x = condition, y = HE_pct, fill = condition)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = cond_colors) +
  theme_classic(base_size = 10) +
  labs(y = "HE-Fibro (%)", x = "", title = "HE-Fibroblast Proportion by Condition") +
  theme(legend.position = "none") +
  geom_text(aes(label = sprintf("%.1f%%", HE_pct)), vjust = -0.5, size = 3)

# 3d: HE-Fibro marker DotPlot
he_dot_genes <- c("MMP1","MMP3","MMP9","CHI3L1","HIF1A","NFKB1","IL6","CXCL8","COL1A1","COL3A1")
he_dot_exist <- he_dot_genes[he_dot_genes %in% rownames(fib)]

# Group by HE_fibro for DotPlot
p3d <- DotPlot(fib, features = he_dot_exist, group.by = "HE_fibro") +
  coord_flip() + theme_classic(base_size = 9) +
  ggtitle("HE-Fibro vs Other Fibroblast Markers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig3 <- (p3a | p3c) / (p3b) / p3d +
  plot_layout(heights = c(1, 1.2, 0.8)) +
  plot_annotation(tag_levels = "A")
ggsave(file.path(figdir, "Figure3.pdf"), fig3, width = 14, height = 16, dpi = 300)
ggsave(file.path(figdir, "Figure3.png"), fig3, width = 14, height = 16, dpi = 150)
cat("Figure 3 done\n")

rm(fib); gc()

# ========== Figure 5: Hypoxia Signature (fix) ==========
cat("Figure 5: Hypoxia Signature (fix)...\n")

mac <- readRDS("/home/moog/test/skin/Phase_output/Phase1/macrophage/macrophage_subclustered.rds")

# Recalculate Hypoxia_sig1 (from Phase4-saved gene list)
hyp_genes <- read.csv("/home/moog/test/skin/Phase_output/Phase4/hypoxia_signature_genes.csv")$x
hyp_in_mac <- hyp_genes[hyp_genes %in% rownames(mac)]
cat(sprintf("Hypoxia signature genes: %d/%d available in scRNA\n", length(hyp_in_mac), length(hyp_genes)))

mac <- AddModuleScore(mac, features = list(hyp_in_mac), name = "Hypoxia_sig")

# 5a: Hypoxia signature gene heatmap (top genes by subcluster)
# Select top 30 genes, plot heatmap of average expression per subcluster
top_hyp <- hyp_in_mac[1:min(30, length(hyp_in_mac))]
avg_expr <- AverageExpression(mac, features = top_hyp, group.by = "seurat_clusters")
mat_hm <- as.matrix(avg_expr$RNA)

# Clean NA/NaN/Inf values
mat_hm[is.na(mat_hm)] <- 0
mat_hm[is.infinite(mat_hm)] <- 0
# Remove all-zero rows
mat_hm <- mat_hm[rowSums(abs(mat_hm)) > 0, , drop = FALSE]

# Save heatmap as separate file
pdf(file.path(figdir, "Figure5a_heatmap.pdf"), width = 8, height = 10)
pheatmap(mat_hm, scale = "row", color = colorRampPalette(c("#377EB8","white","#E41A1C"))(100),
         fontsize_row = 8, fontsize_col = 10,
         main = "Top Hypoxia Response Genes by Macrophage Subcluster",
         cluster_cols = TRUE, cluster_rows = TRUE)
dev.off()

png(file.path(figdir, "Figure5a_heatmap.png"), width = 8, height = 10, units = "in", res = 150)
pheatmap(mat_hm, scale = "row", color = colorRampPalette(c("#377EB8","white","#E41A1C"))(100),
         fontsize_row = 8, fontsize_col = 10,
         main = "Top Hypoxia Response Genes by Macrophage Subcluster",
         cluster_cols = TRUE, cluster_rows = TRUE)
dev.off()

# 5b: Hypoxia score UMAP
p5b <- FeaturePlot(mac, features = "Hypoxia_sig1", pt.size = 0.3) +
  scale_color_viridis_c() + ggtitle("Hypoxia Signature Score") +
  theme(plot.title = element_text(size = 10))

# 5c: Hypoxia score violin by condition
p5c <- VlnPlot(mac, features = "Hypoxia_sig1", group.by = "condition", pt.size = 0,
               cols = cond_colors) +
  ggtitle("Hypoxia Score by Condition") +
  theme(legend.position = "none")

# 5d: C0 characterization (HIF+M1 subcluster)
# Read cluster profile
cluster_prof <- read.csv("/home/moog/test/skin/Phase_output/Phase4/macrophage_cluster_profile.csv")
top_cluster <- cluster_prof$seurat_clusters[1]  # subcluster with highest HIF+NFKB

c0_pct <- mac@meta.data %>%
  filter(seurat_clusters == top_cluster) %>%
  group_by(condition) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(pct = n / sum(n) * 100)

p5d <- ggplot(c0_pct, aes(x = condition, y = pct, fill = condition)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = cond_colors) +
  theme_classic(base_size = 10) +
  labs(y = paste0("% of Cluster ", top_cluster), x = "",
       title = paste0("HIF+M1 Cluster (C", top_cluster, ") Composition")) +
  theme(legend.position = "none") +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.5, size = 3)

# 5e: Subcluster score radar/bar chart
cluster_scores <- cluster_prof %>%
  select(seurat_clusters, HIF, NFKB, M1, M2, Hypoxia) %>%
  tidyr::pivot_longer(cols = -seurat_clusters, names_to = "score_type", values_to = "value") %>%
  mutate(cluster = factor(seurat_clusters))

p5e <- ggplot(cluster_scores, aes(x = cluster, y = value, fill = score_type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 9) +
  labs(x = "Subcluster", y = "Mean Score", fill = "Score",
       title = "Macrophage Subcluster Score Profiles") +
  theme(legend.text = element_text(size = 7))

fig5 <- (p5b | p5c) / (p5d | p5e) + plot_annotation(tag_levels = "A")
ggsave(file.path(figdir, "Figure5.pdf"), fig5, width = 14, height = 10, dpi = 300)
ggsave(file.path(figdir, "Figure5.png"), fig5, width = 14, height = 10, dpi = 150)
cat("Figure 5 done\n")

rm(mac); gc()

# ========== Figure 7: Mechanistic Summary ==========
cat("Figure 7: Mechanistic Summary...\n")

# Summary statistics table
summary_stats <- data.frame(
  Finding = c(
    "Cell types identified",
    "Total cells analyzed",
    "HIF1A-NFKB1 dual+ (Myeloid)",
    "HE-Fibro enrichment (Healing vs NonHealing)",
    "NF-κB ↔ M1 correlation",
    "HIF ↔ M1 correlation",
    "(HIF+NFKB) ↔ (M1-M2) coupling",
    "HIF+M1 subtype (Healing vs NonHealing)",
    "Hypoxia score Healing vs NonHealing (Mac)",
    "HIF1A bulk validation (GSE134431 FC)",
    "HIF1A bulk validation (GSE199939 FC)",
    "HIF1A-NFKB1 bulk correlation",
    "Validated co-regulated targets (≥3 layers)"
  ),
  Value = c(
    "11", "86,137", "34.3%", "4.2x (15.8% vs 3.8%)",
    "ρ = 0.819", "ρ = 0.340", "ρ = 0.721",
    "42.1% vs 9.9%", "p = 2.9e-194",
    "1.91x (p=0.016)", "4.06x (p=5.7e-6)",
    "ρ = 0.634–0.708",
    "7 (VEGFA, IL6, IL1B, PTGS2, BNIP3, ADM, ANGPTL4)"
  ),
  stringsAsFactors = FALSE
)

# Plot summary statistics table
p7 <- ggplot(summary_stats, aes(x = 1, y = rev(seq_len(nrow(summary_stats))))) +
  geom_text(aes(label = Finding, x = 0.3), hjust = 0, size = 3.2, fontface = "bold") +
  geom_text(aes(label = Value, x = 0.75), hjust = 0, size = 3.2, color = "#E41A1C") +
  xlim(0.25, 1.1) +
  theme_void() +
  ggtitle("Key Findings Summary") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.margin = margin(20, 20, 20, 20))

ggsave(file.path(figdir, "Figure7.pdf"), p7, width = 12, height = 8, dpi = 300)
ggsave(file.path(figdir, "Figure7.png"), p7, width = 12, height = 8, dpi = 150)
cat("Figure 7 done\n")

# ========== Supplementary figures ==========
cat("Supplementary figures...\n")

seu <- readRDS("/home/moog/test/skin/Phase_output/Phase1/clustering/seurat_annotated.rds")

# FigS1: QC metrics
p_s1a <- VlnPlot(seu, features = "nFeature_RNA", group.by = "sample", pt.size = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
        legend.position = "none") +
  ggtitle("Genes per Cell by Sample")

p_s1b <- VlnPlot(seu, features = "nCount_RNA", group.by = "sample", pt.size = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
        legend.position = "none") +
  ggtitle("UMI per Cell by Sample")

p_s1c <- VlnPlot(seu, features = "percent.mt", group.by = "sample", pt.size = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
        legend.position = "none") +
  ggtitle("Mitochondrial % by Sample")

figs1 <- p_s1a / p_s1b / p_s1c + plot_annotation(title = "Figure S1: QC Metrics", tag_levels = "A")
ggsave(file.path(figdir, "FigureS1.pdf"), figs1, width = 16, height = 14, dpi = 300)
ggsave(file.path(figdir, "FigureS1.png"), figs1, width = 16, height = 14, dpi = 150)
cat("FigS1 done\n")

# FigS2: Marker gene heatmap
marker_genes <- c(
  "COL1A1","COL3A1","DCN",        # Fibroblast
  "KRT14","KRT5","KRT1",          # Keratinocyte
  "PECAM1","VWF","CDH5",          # Endothelial
  "CD68","CD163","CSF1R",         # Macrophage
  "CD3D","CD3E","CD2",            # T_cell
  "S100A8","S100A9","LYZ",        # Myeloid
  "RGS5","PDGFRB","ACTA2",       # Pericyte/SM
  "KIT","TPSAB1","CPA3",         # Mast_cell
  "MLANA","PMEL","TYRP1",        # Melanocyte
  "PROX1","LYVE1","FLT4"         # Lymphatic
)
marker_exist <- marker_genes[marker_genes %in% rownames(seu)]

p_s2 <- DotPlot(seu, features = marker_exist, group.by = "celltype") +
  coord_flip() + theme_classic(base_size = 8) +
  ggtitle("Canonical Markers by Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(figdir, "FigureS2.pdf"), p_s2, width = 10, height = 12, dpi = 300)
ggsave(file.path(figdir, "FigureS2.png"), p_s2, width = 10, height = 12, dpi = 150)
cat("FigS2 done\n")

# FigS3: Co-expression statistics
coexpr <- read.csv("/home/moog/test/skin/Phase_output/Phase1/expression/coexpression_stats.csv")
p_s3 <- ggplot(coexpr, aes(x = reorder(celltype, -dual_pos), y = dual_pos, fill = cor_spearman)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C", midpoint = 0.13,
                       name = "Spearman ρ") +
  geom_text(aes(label = sprintf("%.1f%%", dual_pos)), vjust = -0.3, size = 3) +
  theme_classic(base_size = 10) +
  labs(x = "", y = "HIF1A+/NFKB1+ Dual Positive (%)",
       title = "HIF1A/NFKB1 Co-expression by Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(figdir, "FigureS3.pdf"), p_s3, width = 10, height = 6, dpi = 300)
ggsave(file.path(figdir, "FigureS3.png"), p_s3, width = 10, height = 6, dpi = 150)
cat("FigS3 done\n")

# FigS4: UMAP faceted by condition
p_s4 <- DimPlot(seu, reduction = "umap", group.by = "celltype", split.by = "condition",
                pt.size = 0.01, cols = ct_colors, ncol = 2) +
  ggtitle("Cell Types by Condition") +
  theme(plot.title = element_text(size = 12, face = "bold"))

ggsave(file.path(figdir, "FigureS4.pdf"), p_s4, width = 14, height = 12, dpi = 300)
ggsave(file.path(figdir, "FigureS4.png"), p_s4, width = 14, height = 12, dpi = 150)
cat("FigS4 done\n")

# ========== Supplementary tables ==========
cat("Supplementary tables...\n")

# TableS1: Sample statistics
ts1 <- seu@meta.data %>%
  group_by(sample, condition) %>%
  summarise(n_cells = n(), median_nFeature = median(nFeature_RNA),
            median_nCount = median(nCount_RNA), median_mt = median(percent.mt), .groups = "drop")
write.csv(ts1, file.path(figdir, "TableS1_sample_stats.csv"), row.names = FALSE)

# TableS2: Cell type x condition statistics
ts2 <- seu@meta.data %>%
  group_by(condition, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(pct = n / sum(n) * 100)
write.csv(ts2, file.path(figdir, "TableS2_celltype_condition.csv"), row.names = FALSE)

# TableS4: Pathway scores
if (file.exists("/home/moog/test/skin/Phase_output/Phase1/pathway/pathway_scores_summary.csv")) {
  ts4 <- read.csv("/home/moog/test/skin/Phase_output/Phase1/pathway/pathway_scores_summary.csv")
  write.csv(ts4, file.path(figdir, "TableS4_pathway_scores.csv"), row.names = FALSE)
}

# TableS5: Hypoxia signature genes
ts5 <- read.csv("/home/moog/test/skin/Phase_output/Phase4/hypoxia_signature_genes.csv")
write.csv(ts5, file.path(figdir, "TableS5_hypoxia_signature.csv"), row.names = FALSE)

rm(seu); gc()

cat("\n=== Step 4 Fix complete: All figures generated ===\n")
