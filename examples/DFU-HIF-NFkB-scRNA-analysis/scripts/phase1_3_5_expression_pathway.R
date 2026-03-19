# ============================================================
# Phase 1.3-1.5: HIF1A/NFKB1 Expression + Pathway Scoring + Differential Analysis
# ============================================================

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(AUCell)
  library(GSEABase)
})

cat("=== Phase 1.3: HIF1A/NFKB1 Expression Analysis ===\n")

outdir_expr <- "/home/moog/test/skin/Phase_output/Phase1/expression"
outdir_path <- "/home/moog/test/skin/Phase_output/Phase1/pathway"
outdir_diff <- "/home/moog/test/skin/Phase_output/Phase1/differential"
dir.create(outdir_expr, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_path, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_diff, recursive = TRUE, showWarnings = FALSE)

seu <- readRDS("/home/moog/test/skin/Phase_output/Phase1/clustering/seurat_annotated.rds")
cat(sprintf("Cell count: %d\n", ncol(seu)))

# ============ 1.3 HIF1A/NFKB1 Expression ============

target_genes <- c("HIF1A", "NFKB1", "RELA", "EPAS1", "HIF1AN")
existing_genes <- target_genes[target_genes %in% rownames(seu)]
cat("Available target genes:", paste(existing_genes, collapse=", "), "\n")

# VlnPlot by celltype
p_vln <- VlnPlot(seu, features = c("HIF1A", "NFKB1"), group.by = "celltype",
                  pt.size = 0, ncol = 1) &
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(outdir_expr, "HIF1A_NFKB1_violin_celltype.png"), p_vln,
       width = 12, height = 10, dpi = 150)

# DotPlot
p_dot <- DotPlot(seu, features = existing_genes, group.by = "celltype") +
  coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("HIF/NF-kB Gene Expression by Cell Type")
ggsave(file.path(outdir_expr, "HIF_NFKB_dotplot.png"), p_dot, width = 10, height = 6, dpi = 150)

# DotPlot faceted by condition
seu@meta.data$ct_cond <- paste(seu@meta.data$celltype, seu@meta.data$condition, sep = "_")
p_dot2 <- DotPlot(seu, features = c("HIF1A", "NFKB1", "RELA"), group.by = "ct_cond") +
  coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
  ggtitle("HIF1A/NFKB1/RELA by CellType×Condition")
ggsave(file.path(outdir_expr, "HIF_NFKB_dotplot_ct_cond.png"), p_dot2,
       width = 18, height = 6, dpi = 150)

# Co-expression analysis
cat("\nCo-expression analysis...\n")
expr_data <- FetchData(seu, vars = c("HIF1A", "NFKB1", "celltype", "condition"))

# HIF1A+/NFKB1+/dual-positive proportions per cell type
coexpr_stats <- expr_data %>%
  group_by(celltype) %>%
  summarise(
    n_cells = n(),
    HIF1A_pos = sum(HIF1A > 0) / n() * 100,
    NFKB1_pos = sum(NFKB1 > 0) / n() * 100,
    dual_pos = sum(HIF1A > 0 & NFKB1 > 0) / n() * 100,
    cor_pearson = cor(HIF1A, NFKB1, method = "pearson"),
    cor_spearman = cor(HIF1A, NFKB1, method = "spearman"),
    .groups = "drop"
  )
write.csv(coexpr_stats, file.path(outdir_expr, "coexpression_stats.csv"), row.names = FALSE)
cat("Co-expression statistics:\n")
print(as.data.frame(coexpr_stats))

# Co-expression by condition
coexpr_cond <- expr_data %>%
  group_by(celltype, condition) %>%
  summarise(
    n = n(),
    dual_pos_pct = sum(HIF1A > 0 & NFKB1 > 0) / n() * 100,
    cor_spearman = cor(HIF1A, NFKB1, method = "spearman"),
    .groups = "drop"
  )
write.csv(coexpr_cond, file.path(outdir_expr, "coexpression_by_condition.csv"), row.names = FALSE)

# Scatter plot: HIF1A vs NFKB1 (colored by celltype)
p_scatter <- ggplot(expr_data %>% filter(HIF1A > 0 | NFKB1 > 0),
                    aes(x = HIF1A, y = NFKB1, color = celltype)) +
  geom_point(size = 0.1, alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  facet_wrap(~celltype, scales = "free") +
  theme_minimal() + ggtitle("HIF1A vs NFKB1 Co-expression")
ggsave(file.path(outdir_expr, "HIF1A_NFKB1_scatter.png"), p_scatter,
       width = 16, height = 12, dpi = 150)

# ============ 1.4 Pathway Activity Scoring ============
cat("\n=== Phase 1.4: Pathway Activity Scoring ===\n")

# HIF pathway gene set
hif_genes <- c("HIF1A", "VEGFA", "SLC2A1", "LDHA", "PGK1", "ENO1", "ALDOA",
               "PKM", "HK2", "PDK1", "BNIP3", "BNIP3L", "CA9", "ADM",
               "ANGPTL4", "LOX", "P4HA1", "P4HA2", "EGLN1", "EGLN3")

# NF-κB pathway gene set
nfkb_genes <- c("NFKB1", "NFKB2", "RELA", "RELB", "REL", "TNF", "IL1B",
                "IL6", "CXCL8", "CCL2", "ICAM1", "VCAM1", "PTGS2",
                "MMP9", "NFKBIA", "NFKBIB", "BCL2", "BIRC3", "TRAF1")

# M1/M2 polarization gene sets
m1_genes <- c("CD86", "NOS2", "IL1B", "TNF", "IL6", "CXCL8", "CCL2",
              "IRF5", "STAT1", "IDO1", "PTGS2", "SOCS3")
m2_genes <- c("CD163", "MRC1", "ARG1", "IL10", "TGFB1", "CCL18",
              "STAB1", "FOLR2", "CD209", "IRF4", "STAT6", "PPARG")

# Filter to genes present in dataset
hif_use <- hif_genes[hif_genes %in% rownames(seu)]
nfkb_use <- nfkb_genes[nfkb_genes %in% rownames(seu)]
m1_use <- m1_genes[m1_genes %in% rownames(seu)]
m2_use <- m2_genes[m2_genes %in% rownames(seu)]

cat(sprintf("HIF genes: %d/%d, NF-kB genes: %d/%d\n",
            length(hif_use), length(hif_genes), length(nfkb_use), length(nfkb_genes)))
cat(sprintf("M1 genes: %d/%d, M2 genes: %d/%d\n",
            length(m1_use), length(m1_genes), length(m2_use), length(m2_genes)))

# AddModuleScore
seu <- AddModuleScore(seu, features = list(hif_use), name = "HIF_score")
seu <- AddModuleScore(seu, features = list(nfkb_use), name = "NFKB_score")
seu <- AddModuleScore(seu, features = list(m1_use), name = "M1_score")
seu <- AddModuleScore(seu, features = list(m2_use), name = "M2_score")

# Pathway score visualization
p_hif <- FeaturePlot(seu, features = "HIF_score1", reduction = "umap", pt.size = 0.1) +
  scale_color_viridis_c() + ggtitle("HIF Pathway Score")
p_nfkb <- FeaturePlot(seu, features = "NFKB_score1", reduction = "umap", pt.size = 0.1) +
  scale_color_viridis_c() + ggtitle("NF-kB Pathway Score")
ggsave(file.path(outdir_path, "pathway_scores_UMAP.png"), p_hif + p_nfkb,
       width = 16, height = 6, dpi = 150)

# Pathway scores per cell type
path_data <- FetchData(seu, vars = c("HIF_score1", "NFKB_score1", "M1_score1", "M2_score1",
                                     "celltype", "condition"))

path_summary <- path_data %>%
  group_by(celltype, condition) %>%
  summarise(
    HIF_mean = mean(HIF_score1),
    NFKB_mean = mean(NFKB_score1),
    M1_mean = mean(M1_score1),
    M2_mean = mean(M2_score1),
    HIF_NFKB_cor = cor(HIF_score1, NFKB_score1, method = "spearman"),
    n = n(),
    .groups = "drop"
  )
write.csv(path_summary, file.path(outdir_path, "pathway_scores_summary.csv"), row.names = FALSE)

# Pathway score VlnPlot by condition (per celltype)
p_path_vln <- VlnPlot(seu, features = c("HIF_score1", "NFKB_score1"),
                       group.by = "celltype", split.by = "condition",
                       pt.size = 0, ncol = 1,
                       cols = c("Healthy"="#4DAF4A","DM"="#377EB8",
                                "Healing"="#FF7F00","NonHealing"="#E41A1C"))
ggsave(file.path(outdir_path, "pathway_scores_violin.png"), p_path_vln,
       width = 16, height = 12, dpi = 150)

# HIF-NF-kB coupling score scatter plot
p_coupling <- ggplot(path_data, aes(x = HIF_score1, y = NFKB_score1, color = condition)) +
  geom_point(size = 0.05, alpha = 0.2) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~celltype, scales = "free") +
  scale_color_manual(values = c("Healthy"="#4DAF4A","DM"="#377EB8",
                                "Healing"="#FF7F00","NonHealing"="#E41A1C")) +
  theme_minimal() + ggtitle("HIF-NF-kB Coupling by Cell Type and Condition")
ggsave(file.path(outdir_path, "HIF_NFKB_coupling.png"), p_coupling,
       width = 16, height = 12, dpi = 150)

# ============ 1.5 Healing vs Non-healing Differential Analysis ============
cat("\n=== Phase 1.5: Healing vs Non-healing Differential Analysis ===\n")

# Cell proportion differences
prop_test <- seu@meta.data %>%
  filter(condition %in% c("Healing", "NonHealing")) %>%
  group_by(condition, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(pct = n / sum(n) * 100) %>%
  tidyr::pivot_wider(names_from = condition, values_from = c(n, pct), values_fill = 0) %>%
  mutate(delta_pct = pct_Healing - pct_NonHealing)

write.csv(prop_test, file.path(outdir_diff, "proportion_healing_vs_nonhealing.csv"),
          row.names = FALSE)
cat("\nCell proportion differences (Healing - NonHealing):\n")
print(as.data.frame(prop_test[, c("celltype", "pct_Healing", "pct_NonHealing", "delta_pct")]))

# Pathway score differences (Healing vs NonHealing)
path_diff <- path_data %>%
  filter(condition %in% c("Healing", "NonHealing")) %>%
  group_by(celltype) %>%
  summarise(
    HIF_healing = mean(HIF_score1[condition == "Healing"]),
    HIF_nonhealing = mean(HIF_score1[condition == "NonHealing"]),
    HIF_pval = wilcox.test(HIF_score1[condition == "Healing"],
                           HIF_score1[condition == "NonHealing"])$p.value,
    NFKB_healing = mean(NFKB_score1[condition == "Healing"]),
    NFKB_nonhealing = mean(NFKB_score1[condition == "NonHealing"]),
    NFKB_pval = wilcox.test(NFKB_score1[condition == "Healing"],
                            NFKB_score1[condition == "NonHealing"])$p.value,
    M1_healing = mean(M1_score1[condition == "Healing"]),
    M1_nonhealing = mean(M1_score1[condition == "NonHealing"]),
    M1_pval = wilcox.test(M1_score1[condition == "Healing"],
                          M1_score1[condition == "NonHealing"])$p.value,
    M2_healing = mean(M2_score1[condition == "Healing"]),
    M2_nonhealing = mean(M2_score1[condition == "NonHealing"]),
    M2_pval = wilcox.test(M2_score1[condition == "Healing"],
                          M2_score1[condition == "NonHealing"])$p.value,
    .groups = "drop"
  )
write.csv(path_diff, file.path(outdir_diff, "pathway_diff_healing_vs_nonhealing.csv"),
          row.names = FALSE)

cat("\nPathway score differences (Healing vs NonHealing):\n")
print(as.data.frame(path_diff[, c("celltype", "HIF_pval", "NFKB_pval", "M1_pval", "M2_pval")]))

# Differential genes per cell type (Healing vs NonHealing)
cat("\nDifferential gene analysis per cell type...\n")
Idents(seu) <- "celltype"
key_celltypes <- c("Macrophage", "Fibroblast", "Keratinocyte", "Endothelial", "T_cell", "Myeloid")

for (ct in key_celltypes) {
  cat(sprintf("  %s: ", ct))
  sub <- subset(seu, celltype == ct & condition %in% c("Healing", "NonHealing"))
  Idents(sub) <- "condition"
  tryCatch({
    degs <- FindMarkers(sub, ident.1 = "Healing", ident.2 = "NonHealing",
                        min.pct = 0.1, logfc.threshold = 0.25, test.use = "wilcox")
    degs$gene <- rownames(degs)
    degs$celltype <- ct
    write.csv(degs, file.path(outdir_diff, paste0("DEGs_", ct, "_H_vs_NH.csv")))
    cat(sprintf("%d DEGs (up: %d, down: %d)\n",
                nrow(degs), sum(degs$avg_log2FC > 0 & degs$p_val_adj < 0.05),
                sum(degs$avg_log2FC < 0 & degs$p_val_adj < 0.05)))
  }, error = function(e) cat("Skipping:", conditionMessage(e), "\n"))
}

# Save
saveRDS(seu, "/home/moog/test/skin/Phase_output/Phase1/clustering/seurat_annotated.rds")

cat("\n=== Phase 1.3-1.5 Complete ===\n")
