# ============================================================
# Phase 1.2: Cell Type Annotation
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(future)

plan("multicore", workers = 16)
options(future.globals.maxSize = 100 * 1024^3)

cat("=== Phase 1.2: Cell Type Annotation ===\n")

outdir <- "/home/moog/test/skin/Phase_output/Phase1/clustering"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load integrated object
seu <- readRDS(file.path(outdir, "seurat_integrated.rds"))
cat(sprintf("Cell count: %d, Cluster count: %d\n", ncol(seu), length(levels(Idents(seu)))))

# --- JoinLayers (Seurat v5) ---
cat("JoinLayers...\n")
seu <- JoinLayers(seu)

# --- FindAllMarkers ---
cat("Finding marker genes for each cluster...\n")
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,
                          test.use = "wilcox")
write.csv(markers, file.path(outdir, "all_markers.csv"), row.names = FALSE)

# Top5 markers per cluster
top5 <- markers %>% group_by(cluster) %>% top_n(5, wt = avg_log2FC)
cat("Top5 markers per cluster (first 10 clusters):\n")
print(top5 %>% filter(as.numeric(as.character(cluster)) < 10) %>%
        select(cluster, gene, avg_log2FC, pct.1, pct.2))

# --- Known cell type markers ---
cell_markers <- list(
  "Keratinocyte"     = c("KRT14", "KRT5", "KRT1", "KRT10", "KRT6A", "KRT16", "KRT17"),
  "Fibroblast"       = c("COL1A1", "COL3A1", "DCN", "LUM", "VIM", "PDGFRA"),
  "Endothelial"      = c("PECAM1", "VWF", "CDH5", "CLDN5", "EMCN"),
  "Macrophage"       = c("CD68", "CD14", "CSF1R", "MARCO", "C1QA", "C1QB"),
  "T_cell"           = c("CD3D", "CD3E", "CD2", "IL7R", "CD8A"),
  "B_cell"           = c("CD79A", "MS4A1", "CD19"),
  "Neutrophil"       = c("CSF3R", "S100A8", "S100A9", "FCGR3B"),
  "Mast_cell"        = c("TPSB2", "CPA3", "KIT", "TPSAB1"),
  "Pericyte_SMC"     = c("RGS5", "ACTA2", "MYH11", "PDGFRB"),
  "Melanocyte"       = c("PMEL", "MLANA", "TYR", "DCT"),
  "DC"               = c("CD1A", "CD1C", "CLEC10A", "FCER1A"),
  "Lymphatic_Endo"   = c("PROX1", "LYVE1", "FLT4"),
  "NK_cell"          = c("GNLY", "NKG7", "KLRD1")
)

# DotPlot
all_markers_flat <- unique(unlist(cell_markers))
existing <- all_markers_flat[all_markers_flat %in% rownames(seu)]

p_dot <- DotPlot(seu, features = existing, cluster.idents = TRUE) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7)) +
  ggtitle("Cell Type Marker Expression")
ggsave(file.path(outdir, "marker_dotplot.png"), p_dot, width = 18, height = 14, dpi = 150)

# --- Calculate module score for each cluster ---
cat("\nCalculating marker gene module scores...\n")
for (ct in names(cell_markers)) {
  genes_use <- cell_markers[[ct]][cell_markers[[ct]] %in% rownames(seu)]
  if (length(genes_use) >= 2) {
    seu <- AddModuleScore(seu, features = list(genes_use), name = paste0(ct, "_score"))
  }
}

# Summarize scores for each cluster
score_cols <- grep("_score1$", colnames(seu@meta.data), value = TRUE)
cluster_scores <- seu@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols), mean), n_cells = n(), .groups = "drop")

cat("\nCluster module score matrix:\n")
print(as.data.frame(cluster_scores[, c("seurat_clusters", "n_cells",
                                        head(score_cols, 6))]))

# --- Auto annotation ---
cat("\nAuto annotation...\n")

# Find the highest-scoring cell type for each cluster
assign_celltype <- function(row, score_cols) {
  scores <- as.numeric(row[score_cols])
  names(scores) <- gsub("_score1$", "", score_cols)
  top <- names(sort(scores, decreasing = TRUE))[1]
  return(top)
}

cluster_annotations <- data.frame(
  cluster = cluster_scores$seurat_clusters,
  n_cells = cluster_scores$n_cells
)

# Assign the highest-scoring cell type to each cluster
for (i in 1:nrow(cluster_scores)) {
  scores <- as.numeric(cluster_scores[i, score_cols])
  names(scores) <- gsub("_score1$", "", score_cols)
  cluster_annotations$auto_type[i] <- names(sort(scores, decreasing = TRUE))[1]
  cluster_annotations$top_score[i] <- max(scores)
  cluster_annotations$second_type[i] <- names(sort(scores, decreasing = TRUE))[2]
  cluster_annotations$second_score[i] <- sort(scores, decreasing = TRUE)[2]
}

# Manual correction based on marker gene inspection
# Check top marker genes for each cluster to assist judgment
for (i in 1:nrow(cluster_annotations)) {
  cl <- as.character(cluster_annotations$cluster[i])
  cl_markers <- markers %>% filter(cluster == cl) %>% top_n(5, wt = avg_log2FC)
  cluster_annotations$top_markers[i] <- paste(cl_markers$gene, collapse = ",")
}

cat("\nPreliminary annotation results:\n")
print(cluster_annotations[, c("cluster", "n_cells", "auto_type", "top_score", "top_markers")])

write.csv(cluster_annotations, file.path(outdir, "cluster_annotations_auto.csv"), row.names = FALSE)

# --- Fine annotation based on scores + markers ---
# Save score matrix for subsequent manual validation
write.csv(as.data.frame(cluster_scores), file.path(outdir, "cluster_module_scores.csv"),
          row.names = FALSE)

# Apply annotation
annotation_map <- setNames(cluster_annotations$auto_type,
                           as.character(cluster_annotations$cluster))
seu$celltype_auto <- annotation_map[as.character(seu$seurat_clusters)]

# UMAP by celltype
p_ct <- DimPlot(seu, reduction = "umap", group.by = "celltype_auto",
                label = TRUE, repel = TRUE, pt.size = 0.1) +
  ggtitle("Auto Cell Type Annotation")
ggsave(file.path(outdir, "UMAP_celltype_auto.png"), p_ct, width = 12, height = 8, dpi = 150)

# Cell type proportions
ct_prop <- seu@meta.data %>%
  group_by(condition, celltype_auto) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(pct = n / sum(n) * 100)

p_bar <- ggplot(ct_prop, aes(x = condition, y = pct, fill = celltype_auto)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(y = "Percentage (%)", x = "", fill = "Cell Type") +
  ggtitle("Cell Type Proportions by Condition")
ggsave(file.path(outdir, "celltype_proportion_bar.png"), p_bar, width = 10, height = 6, dpi = 150)

# Save
cat("\nSaving annotated Seurat object...\n")
saveRDS(seu, file.path(outdir, "seurat_annotated.rds"))

cat("\n=== Phase 1.2 Complete ===\n")

# Output cell type distribution
cat("\nCell type distribution:\n")
print(table(seu$celltype_auto))
