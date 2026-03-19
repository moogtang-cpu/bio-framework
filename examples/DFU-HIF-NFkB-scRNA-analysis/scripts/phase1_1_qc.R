# ============================================================
# Phase 1.1: Data Preprocessing and QC
# GSE165816 scRNA-seq — Foot skin samples
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(future)

# Parallel settings
plan("multicore", workers = 16)
options(future.globals.maxSize = 100 * 1024^3)

cat("=== Phase 1.1: Data Preprocessing and QC ===\n")

# Output directory
outdir <- "/home/moog/test/skin/Phase_output/Phase1/QC"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- Sample metadata ---
sample_info <- data.frame(
  GSM = paste0("GSM", 5050521:5050574),
  SampleID = paste0("G", c(1, "1A", 2, "2A", 3, "3A", 4, "4A", 5, 6, 7, 8, 9,
                           10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                           23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
                           36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)),
  Tissue = c("Forearm","Foot","Foot","Foot","Foot","Foot","Foot","Forearm","Foot","Foot",
             "Foot","Foot","Foot","Foot","Forearm","Forearm","Forearm","Foot","Foot","Foot",
             "Forearm","Forearm","PBMC","PBMC","PBMC","PBMC","Foot","Foot","PBMC","PBMC",
             "PBMC","Foot","Forearm","PBMC","Foot","Foot","Foot","Foot","Forearm","Forearm",
             "PBMC","Foot","Foot","Foot","Foot","Foot","Foot","Foot","Foot","Foot",
             "PBMC","Forearm","Foot","Foot"),
  Condition = c("Healthy","DM","Healing","DM","DM","DM","Healing","DM","DM","NonHealing",
                "Healing","Healing","NonHealing","Healthy","NonHealing","Healing","Healing",
                "Healthy","Healing","Healthy","Healing","Healthy","DM","Healthy","NonHealing",
                "Healing","Healing","Healthy","Healing","Healthy","Healing","Healthy","Healthy",
                "Healthy","Healthy","Healthy","NonHealing","NonHealing","NonHealing","Healthy",
                "NonHealing","DM","NonHealing","Healthy","DM","Healing","Healthy","Healthy",
                "Healing","DM","DM","DM","Healing","Healthy"),
  stringsAsFactors = FALSE
)

# Select foot skin samples only
foot_samples <- sample_info[sample_info$Tissue == "Foot", ]
cat(sprintf("Number of foot skin samples: %d\n", nrow(foot_samples)))
cat(sprintf("  Healthy: %d | DM: %d | Healing: %d | NonHealing: %d\n",
            sum(foot_samples$Condition == "Healthy"),
            sum(foot_samples$Condition == "DM"),
            sum(foot_samples$Condition == "Healing"),
            sum(foot_samples$Condition == "NonHealing")))

# --- Batch load data ---
data_dir <- "/home/moog/test/skin/data/GSE165816"
seurat_list <- list()

for (i in 1:nrow(foot_samples)) {
  gsm <- foot_samples$GSM[i]
  sid <- foot_samples$SampleID[i]
  cond <- foot_samples$Condition[i]

  # Build filename
  fname <- list.files(data_dir, pattern = paste0(gsm, "_"), full.names = TRUE)
  if (length(fname) == 0) {
    cat(sprintf("  Warning: %s (%s) file not found, skipping\n", gsm, sid))
    next
  }
  fname <- fname[grep("\\.csv\\.gz$", fname)]
  if (length(fname) == 0) next

  cat(sprintf("  Reading %s (%s, %s)... ", gsm, sid, cond))

  # Read CSV: first row is barcodes (N columns), data rows are gene + N counts (N+1 columns)
  # Read barcode row first
  con <- gzfile(fname[1], "r")
  barcode_line <- readLines(con, n = 1)
  close(con)
  barcodes <- strsplit(barcode_line, ",")[[1]]

  # Read data rows (skip barcode row)
  mat <- data.table::fread(fname[1], header = FALSE, skip = 1)
  genes <- make.unique(as.character(mat[[1]]))
  mat <- as.matrix(mat[, -1, with = FALSE])
  rownames(mat) <- genes
  colnames(mat) <- barcodes

  # Create Seurat object
  obj <- CreateSeuratObject(counts = mat, project = sid, min.cells = 3, min.features = 200)
  obj$sample <- sid
  obj$condition <- cond
  obj$gsm <- gsm

  # Mitochondrial percentage
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  cat(sprintf("%d cells, %d genes\n", ncol(obj), nrow(obj)))

  seurat_list[[sid]] <- obj
  rm(mat, obj); gc(verbose = FALSE)
}

cat(sprintf("\nSuccessfully loaded %d samples\n", length(seurat_list)))

# --- Merge ---
cat("Merging all samples...\n")
combined <- merge(seurat_list[[1]], y = seurat_list[-1],
                  add.cell.ids = names(seurat_list))
cat(sprintf("After merge: %d cells, %d genes\n", ncol(combined), nrow(combined)))

# --- QC metric summary ---
qc_stats <- combined@meta.data %>%
  group_by(sample, condition) %>%
  summarise(
    n_cells = n(),
    median_nFeature = median(nFeature_RNA),
    median_nCount = median(nCount_RNA),
    median_pct_mt = median(percent.mt),
    .groups = "drop"
  )
write.csv(qc_stats, file.path(outdir, "qc_stats_pre_filter.csv"), row.names = FALSE)
cat("\nQC statistics (before filtering):\n")
print(as.data.frame(qc_stats))

# --- QC visualization ---
cat("\nGenerating QC plots...\n")

# VlnPlot
p1 <- VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               group.by = "condition", pt.size = 0, ncol = 3) &
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(outdir, "QC_violin_by_condition.png"), p1, width = 14, height = 5, dpi = 150)

p2 <- VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               group.by = "sample", pt.size = 0, ncol = 3) &
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
ggsave(file.path(outdir, "QC_violin_by_sample.png"), p2, width = 20, height = 5, dpi = 150)

# Scatter plots
p3 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                     group.by = "condition") + NoLegend()
p4 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt",
                     group.by = "condition") + NoLegend()
ggsave(file.path(outdir, "QC_scatter.png"), p3 + p4, width = 12, height = 5, dpi = 150)

# --- Filtering ---
cat("\nApplying QC filters...\n")
cat(sprintf("Before filtering: %d cells\n", ncol(combined)))

combined <- subset(combined,
                   subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                            nCount_RNA > 500 & nCount_RNA < 40000 &
                            percent.mt < 20)

cat(sprintf("After filtering: %d cells\n", ncol(combined)))

# Post-filter statistics
qc_post <- combined@meta.data %>%
  group_by(sample, condition) %>%
  summarise(n_cells = n(), .groups = "drop")
write.csv(qc_post, file.path(outdir, "cell_counts_post_filter.csv"), row.names = FALSE)

cat("\nCell counts per condition:\n")
print(combined@meta.data %>% group_by(condition) %>% summarise(n = n(), .groups = "drop"))

# --- Save QC-filtered object ---
cat("\nSaving Seurat object...\n")
saveRDS(combined, file.path(outdir, "seurat_qc_filtered.rds"))

cat("\n=== Phase 1.1 QC Complete ===\n")

# Cleanup
rm(seurat_list); gc()
