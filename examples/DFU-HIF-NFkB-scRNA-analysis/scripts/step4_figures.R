# ============================================================
# Step 4: Publication-quality figures
# Figure 1-7 key panels
# ============================================================

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(pheatmap)
})

cat("=== Step 4: Generating publication-quality figures ===\n")

figdir <- "/home/moog/test/skin/Phase_output/publication_figures"
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

seu <- readRDS("/home/moog/test/skin/Phase_output/Phase1/clustering/seurat_annotated.rds")

ct_colors <- c(
  "Fibroblast"="#E41A1C","Keratinocyte"="#377EB8","Endothelial"="#4DAF4A",
  "Macrophage"="#984EA3","T_cell"="#FF7F00","Myeloid"="#A65628",
  "Pericyte"="#F781BF","Smooth_Muscle"="#999999","Mast_cell"="#66C2A5",
  "Melanocyte"="#FC8D62","Lymphatic_Endo"="#8DA0CB"
)
cond_colors <- c("Healthy"="#4DAF4A","DM"="#377EB8","Healing"="#FF7F00","NonHealing"="#E41A1C")

# ========== Figure 1: Overview ==========
cat("Figure 1...\n")

p1b <- DimPlot(seu, reduction="umap", group.by="celltype", label=TRUE, repel=TRUE,
               pt.size=0.05, cols=ct_colors, label.size=3) +
  ggtitle("Cell Types") + NoLegend() +
  theme(plot.title=element_text(size=12, face="bold"))

p1c <- DimPlot(seu, reduction="umap", group.by="condition", pt.size=0.05,
               cols=cond_colors) +
  ggtitle("Condition") +
  theme(plot.title=element_text(size=12, face="bold"))

# Cell type proportions
ct_prop <- seu@meta.data %>%
  group_by(condition, celltype) %>% summarise(n=n(), .groups="drop") %>%
  group_by(condition) %>% mutate(pct=n/sum(n)*100)

p1d <- ggplot(ct_prop, aes(x=condition, y=pct, fill=celltype)) +
  geom_bar(stat="identity", position="stack", width=0.7) +
  scale_fill_manual(values=ct_colors) +
  theme_classic(base_size=10) +
  labs(y="Percentage (%)", x="", fill="Cell Type") +
  theme(legend.position="right", legend.text=element_text(size=7),
        legend.key.size=unit(0.3, "cm"))

fig1 <- (p1b | p1c) / p1d + plot_annotation(tag_levels="A")
ggsave(file.path(figdir, "Figure1.pdf"), fig1, width=14, height=12, dpi=300)
ggsave(file.path(figdir, "Figure1.png"), fig1, width=14, height=12, dpi=150)

# ========== Figure 2: HIF/NFKB Expression ==========
cat("Figure 2...\n")

genes_plot <- c("HIF1A","NFKB1","RELA","EPAS1","HIF1AN")
genes_exist <- genes_plot[genes_plot %in% rownames(seu)]

p2a <- DotPlot(seu, features=genes_exist, group.by="celltype") +
  coord_flip() + theme_classic(base_size=9) +
  ggtitle("HIF/NF-κB Gene Expression") +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# Co-expression proportion heatmap
coexpr <- read.csv("/home/moog/test/skin/Phase_output/Phase1/expression/coexpression_stats.csv")
p2c <- ggplot(coexpr, aes(x=reorder(celltype, -dual_pos), y=1, fill=dual_pos)) +
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradient2(low="white", high="red", mid="orange", midpoint=15,
                       name="Dual+ (%)") +
  geom_text(aes(label=sprintf("%.1f%%", dual_pos)), size=2.5) +
  theme_minimal(base_size=9) +
  theme(axis.text.x=element_text(angle=45, hjust=1), axis.text.y=element_blank()) +
  labs(x="", y="", title="HIF1A+/NFKB1+ Co-expressing Cells")

# HIF/NFKB score UMAPs
p2d <- FeaturePlot(seu, features="HIF_score1", pt.size=0.05) +
  scale_color_viridis_c() + ggtitle("HIF Pathway Score") +
  theme(plot.title=element_text(size=10))
p2e <- FeaturePlot(seu, features="NFKB_score1", pt.size=0.05) +
  scale_color_viridis_c() + ggtitle("NF-κB Pathway Score") +
  theme(plot.title=element_text(size=10))

fig2 <- (p2a | p2c) / (p2d | p2e) + plot_annotation(tag_levels="A")
ggsave(file.path(figdir, "Figure2.pdf"), fig2, width=14, height=10, dpi=300)
ggsave(file.path(figdir, "Figure2.png"), fig2, width=14, height=10, dpi=150)

# ========== Figure 4: Macrophage Polarization ==========
cat("Figure 4...\n")

mac <- readRDS("/home/moog/test/skin/Phase_output/Phase1/macrophage/macrophage_subclustered.rds")

p4a <- DimPlot(mac, reduction="umap", label=TRUE, repel=TRUE, pt.size=0.3) +
  ggtitle("Macrophage Subclusters") + NoLegend()

p4b <- DimPlot(mac, reduction="umap", group.by="condition", pt.size=0.3, cols=cond_colors) +
  ggtitle("Condition")

mac_data <- FetchData(mac, vars=c("M1_score1","M2_score1","HIF_score1","NFKB_score1","condition"))

p4c <- ggplot(mac_data, aes(x=M1_score1, y=M2_score1, color=condition)) +
  geom_point(size=0.2, alpha=0.3) +
  geom_smooth(method="lm", se=FALSE, linewidth=0.8) +
  scale_color_manual(values=cond_colors) +
  theme_classic(base_size=10) +
  labs(x="M1 Score", y="M2 Score", title="M1 vs M2 Polarization")

p4d <- ggplot(mac_data, aes(x=HIF_score1+NFKB_score1, y=M1_score1-M2_score1, color=condition)) +
  geom_point(size=0.2, alpha=0.3) +
  geom_smooth(method="lm", se=FALSE, linewidth=0.8) +
  scale_color_manual(values=cond_colors) +
  theme_classic(base_size=10) +
  labs(x="HIF + NF-κB Score", y="M1 - M2 Bias", title="Pathway-Polarization Coupling (ρ=0.721)")

# Polarization boxplot
mac_polar_long <- mac_data %>%
  tidyr::pivot_longer(cols=c(M1_score1, M2_score1), names_to="polarization", values_to="score") %>%
  mutate(polarization = gsub("_score1", "", polarization))

p4e <- ggplot(mac_polar_long %>% filter(condition %in% c("Healing","NonHealing")),
              aes(x=condition, y=score, fill=polarization)) +
  geom_boxplot(outlier.size=0.3, width=0.6) +
  scale_fill_manual(values=c("M1"="#E41A1C","M2"="#377EB8")) +
  theme_classic(base_size=10) +
  labs(y="Score", x="", fill="Polarization", title="M1/M2 in Healing vs NonHealing")

fig4 <- (p4a | p4b) / (p4c | p4d | p4e) + plot_annotation(tag_levels="A")
ggsave(file.path(figdir, "Figure4.pdf"), fig4, width=16, height=10, dpi=300)
ggsave(file.path(figdir, "Figure4.png"), fig4, width=16, height=10, dpi=150)

# ========== Figure 5: Hypoxia Signature ==========
cat("Figure 5...\n")

p5b <- FeaturePlot(mac, features="Hypoxia_sig1", pt.size=0.3) +
  scale_color_viridis_c() + ggtitle("Hypoxia Signature Score")

p5c <- VlnPlot(mac, features="Hypoxia_sig1", group.by="condition", pt.size=0,
               cols=cond_colors) +
  ggtitle("Hypoxia Score by Condition") +
  theme(legend.position="none")

# C0 characterization
c0_pct <- mac@meta.data %>%
  filter(seurat_clusters == 0) %>%
  group_by(condition) %>%
  summarise(n=n(), .groups="drop") %>%
  mutate(pct=n/sum(n)*100)

p5e <- ggplot(c0_pct, aes(x=condition, y=pct, fill=condition)) +
  geom_bar(stat="identity", width=0.6) +
  scale_fill_manual(values=cond_colors) +
  theme_classic(base_size=10) +
  labs(y="% of Cluster 0", x="", title="HIF+M1 Cluster (C0) Composition") +
  theme(legend.position="none")

fig5 <- (p5b | p5c | p5e) + plot_annotation(tag_levels="A")
ggsave(file.path(figdir, "Figure5.pdf"), fig5, width=16, height=5, dpi=300)
ggsave(file.path(figdir, "Figure5.png"), fig5, width=16, height=5, dpi=150)

# ========== Figure 6: Bulk Validation ==========
cat("Figure 6...\n")

# Read Phase2 data
deg_199 <- read.csv("/home/moog/test/skin/Phase_output/Phase2/DEGs_GSE199939.csv")

# Volcano plot
p6e <- ggplot(deg_199, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(aes(color = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1,
                                ifelse(logFC > 0, "Up", "Down"), "NS")),
             size=0.3, alpha=0.5) +
  scale_color_manual(values=c("Up"="#E41A1C","Down"="#377EB8","NS"="grey70"), name="") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50") +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey50") +
  theme_classic(base_size=10) +
  labs(x="log2 Fold Change", y="-log10(FDR)", title="GSE199939: DFU vs Normal")

# Annotate key genes
key_genes <- c("HIF1A","NFKB1","IL1B","CXCL8","VEGFA","PTGS2","LDHA","BNIP3","ADM")
deg_key <- deg_199[deg_199$gene %in% key_genes, ]
if (nrow(deg_key) > 0) {
  p6e <- p6e + ggrepel::geom_text_repel(
    data=deg_key, aes(label=gene), size=2.5, max.overlaps=20
  )
}

ggsave(file.path(figdir, "Figure6.pdf"), p6e, width=8, height=6, dpi=300)
ggsave(file.path(figdir, "Figure6.png"), p6e, width=8, height=6, dpi=150)

# ========== Supplementary Table ==========
cat("Supplementary tables...\n")

# TableS1: Cell counts per sample
ts1 <- seu@meta.data %>%
  group_by(sample, condition) %>%
  summarise(n_cells=n(), median_nFeature=median(nFeature_RNA),
            median_nCount=median(nCount_RNA), median_mt=median(percent.mt), .groups="drop")
write.csv(ts1, file.path(figdir, "TableS1_sample_stats.csv"), row.names=FALSE)

# TableS4: Pathway scores
ts4 <- read.csv("/home/moog/test/skin/Phase_output/Phase1/pathway/pathway_scores_summary.csv")
write.csv(ts4, file.path(figdir, "TableS4_pathway_scores.csv"), row.names=FALSE)

# TableS5: Hypoxia signature
ts5 <- read.csv("/home/moog/test/skin/Phase_output/Phase4/hypoxia_signature_genes.csv")
write.csv(ts5, file.path(figdir, "TableS5_hypoxia_signature.csv"), row.names=FALSE)

cat("\n=== Step 4 complete: Figures generated ===\n")
