# ============================================================
# Phase 1.2b: Apply Cell Type Annotation (Fix)
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

cat("=== Phase 1.2b: Apply Annotation ===\n")

outdir <- "/home/moog/test/skin/Phase_output/Phase1/clustering"

# Load integrated object (original version before JoinLayers)
seu <- readRDS(file.path(outdir, "seurat_integrated.rds"))
seu <- JoinLayers(seu)

# Read auto-annotation results
annot <- read.csv(file.path(outdir, "cluster_annotations_auto.csv"))

# Refined annotation - manual correction based on marker genes
# Review top markers for each cluster to fine-tune
cat("Refining annotation based on marker genes...\n")

# Define refined annotation mapping
celltype_map <- c(
  "0"  = "Fibroblast",          # PCOLCE2,PI16,SFRP4 -> Reticular fibroblasts
  "1"  = "Fibroblast",          # IGF1,C3,C7 -> Secretory fibroblasts
  "2"  = "Pericyte",            # RGS5,FAM162B,HIGD1B
  "3"  = "Fibroblast",          # FHL2,FHL5,SSTR2
  "4"  = "Smooth_Muscle",       # RERGL,PLN,CASQ2 -> Smooth muscle
  "5"  = "T_cell",              # CD3E,KLRB1,GZMA,NKG7
  "6"  = "Keratinocyte",        # KRT1,DSC1,KRT2 -> Spinous layer keratinocytes
  "7"  = "Endothelial",         # ACKR1,SELE -> Venous endothelial
  "8"  = "Endothelial",         # RBP7,GPIHBP1 -> Capillary endothelial
  "9"  = "Myeloid",             # IL1B,LGALS2,OLR1,CD1C -> Monocyte/DC
  "10" = "Fibroblast",          # NKD2,DIO2,COL23A1 -> Papillary fibroblasts
  "11" = "Keratinocyte",        # COL17A1,ALDH3A1 -> Basal layer
  "12" = "Macrophage",          # C1QB,FOLR2,CD209,CCL18 -> Tissue-resident macrophages
  "13" = "Mast_cell",           # TPSAB1,TPSB2,CPA3
  "14" = "Fibroblast",          # MMP1,MMP3,IL11 -> HE-Fibro inflammatory fibroblasts
  "15" = "Keratinocyte",        # KRT7,KRT19 -> Glandular/sweat duct epithelium
  "16" = "Keratinocyte",        # UBE2C,PBK -> Proliferating keratinocytes
  "17" = "Smooth_Muscle",       # DES,P2RX1 -> Vascular smooth muscle
  "18" = "Keratinocyte",        # KRT6A,KRT6B,KRT16 -> Activated keratinocytes
  "19" = "Pericyte",            # RGS5,LMOD1
  "20" = "Endothelial",         # SOX18,PLVAP -> Capillary
  "21" = "Fibroblast",          # MEG3,LRP1,THBS2
  "22" = "Melanocyte",          # DCT,MLANA,TYR
  "23" = "Fibroblast",          # APOD,FOXD1
  "24" = "Lymphatic_Endo",      # CCL21,PROX1
  "25" = "Fibroblast",          # MFAP5,OGN,WNT2
  "26" = "Keratinocyte",        # DSC1,DMKN
  "27" = "T_cell",              # Proliferating T cells (FAM111B,CDCA5)
  "28" = "T_cell"               # CD3D,TRBC2 (small cluster)
)

# Apply annotation
seu$celltype <- celltype_map[as.character(seu$seurat_clusters)]

cat("\nCell type distribution:\n")
print(sort(table(seu$celltype), decreasing = TRUE))

# --- Visualization ---
cat("\nGenerating visualizations...\n")

# Color definitions
ct_colors <- c(
  "Fibroblast"     = "#E41A1C",
  "Keratinocyte"   = "#377EB8",
  "Endothelial"    = "#4DAF4A",
  "Macrophage"     = "#984EA3",
  "T_cell"         = "#FF7F00",
  "Myeloid"        = "#A65628",
  "Pericyte"       = "#F781BF",
  "Smooth_Muscle"  = "#999999",
  "Mast_cell"      = "#66C2A5",
  "Melanocyte"     = "#FC8D62",
  "Lymphatic_Endo" = "#8DA0CB"
)

# UMAP by celltype
p1 <- DimPlot(seu, reduction = "umap", group.by = "celltype",
              label = TRUE, repel = TRUE, pt.size = 0.1, cols = ct_colors) +
  ggtitle("Cell Type Annotation") + NoLegend()

p2 <- DimPlot(seu, reduction = "umap", group.by = "celltype",
              label = FALSE, pt.size = 0.1, cols = ct_colors) +
  ggtitle("Cell Types")

ggsave(file.path(outdir, "UMAP_celltype.png"), p1 + p2, width = 20, height = 8, dpi = 150)

# Cell type proportions per condition
ct_prop <- seu@meta.data %>%
  group_by(condition, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(pct = n / sum(n) * 100)

p3 <- ggplot(ct_prop, aes(x = condition, y = pct, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = ct_colors) +
  theme_minimal(base_size = 12) +
  labs(y = "Percentage (%)", x = "", fill = "Cell Type") +
  ggtitle("Cell Type Proportions by Condition")
ggsave(file.path(outdir, "celltype_proportion_stacked.png"), p3, width = 10, height = 6, dpi = 150)

# Condition x cell type count table
ct_table <- seu@meta.data %>%
  group_by(condition, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0)
write.csv(ct_table, file.path(outdir, "celltype_condition_table.csv"), row.names = FALSE)
cat("\nCell type x condition:\n")
print(as.data.frame(ct_table))

# --- Save ---
saveRDS(seu, file.path(outdir, "seurat_annotated.rds"))
cat("\n=== Phase 1.2 Complete ===\n")
