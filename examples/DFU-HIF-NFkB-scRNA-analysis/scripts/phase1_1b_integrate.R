# ============================================================
# Phase 1.1b: Normalization, Integration, Dimensionality Reduction & Clustering
# ============================================================

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)
library(future)

plan("multicore", workers = 16)
options(future.globals.maxSize = 100 * 1024^3)

cat("=== Phase 1.1b: Normalization, Integration, Clustering ===\n")

outdir <- "/home/moog/test/skin/Phase_output/Phase1/clustering"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load QC-filtered object
cat("Loading QC-filtered Seurat object...\n")
seu <- readRDS("/home/moog/test/skin/Phase_output/Phase1/QC/seurat_qc_filtered.rds")
cat(sprintf("Cell count: %d, Gene count: %d\n", ncol(seu), nrow(seu)))

# --- Normalization ---
cat("LogNormalize normalization...\n")
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

# --- Highly Variable Genes ---
cat("Finding highly variable genes...\n")
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)

# --- Scale ---
cat("Scaling data...\n")
seu <- ScaleData(seu, vars.to.regress = "percent.mt")

# --- PCA ---
cat("Running PCA...\n")
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

# PCA ElbowPlot
p_elbow <- ElbowPlot(seu, ndims = 50) + ggtitle("PCA Elbow Plot")
ggsave(file.path(outdir, "PCA_elbow.png"), p_elbow, width = 8, height = 5, dpi = 150)

# --- Harmony Integration ---
cat("Harmony batch correction (by sample)...\n")
seu <- RunHarmony(seu, group.by.vars = "sample", dims.use = 1:30,
                  max.iter.harmony = 20, verbose = TRUE)

# --- UMAP ---
cat("Running UMAP...\n")
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30)

# --- Clustering ---
cat("FindNeighbors + FindClusters...\n")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30)
seu <- FindClusters(seu, resolution = 0.8)

cat(sprintf("Number of clusters: %d\n", length(levels(Idents(seu)))))

# --- Visualization ---
cat("Generating visualizations...\n")

# UMAP by cluster
p1 <- DimPlot(seu, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.1) +
  ggtitle("Clusters (res=0.8)") + NoLegend()

# UMAP by condition
p2 <- DimPlot(seu, reduction = "umap", group.by = "condition", pt.size = 0.1,
              cols = c("Healthy" = "#4DAF4A", "DM" = "#377EB8",
                       "Healing" = "#FF7F00", "NonHealing" = "#E41A1C")) +
  ggtitle("Condition")

# UMAP by sample
p3 <- DimPlot(seu, reduction = "umap", group.by = "sample", pt.size = 0.1) +
  ggtitle("Sample") + NoLegend()

ggsave(file.path(outdir, "UMAP_clusters.png"), p1, width = 10, height = 8, dpi = 150)
ggsave(file.path(outdir, "UMAP_condition.png"), p2, width = 12, height = 8, dpi = 150)
ggsave(file.path(outdir, "UMAP_sample.png"), p3, width = 10, height = 8, dpi = 150)

# Combined plot
p_combined <- (p1 | p2) / (p3 | plot_spacer())
ggsave(file.path(outdir, "UMAP_overview.png"), p_combined, width = 20, height = 16, dpi = 150)

# Cell counts per cluster
cluster_stats <- seu@meta.data %>%
  group_by(seurat_clusters, condition) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0)
write.csv(cluster_stats, file.path(outdir, "cluster_condition_counts.csv"), row.names = FALSE)

cat("\nCell counts per cluster:\n")
print(table(Idents(seu)))

# --- Save ---
cat("\nSaving integrated Seurat object...\n")
saveRDS(seu, file.path(outdir, "seurat_integrated.rds"))

cat("\n=== Phase 1.1b Complete ===\n")
