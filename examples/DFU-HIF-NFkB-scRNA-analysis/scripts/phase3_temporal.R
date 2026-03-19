# ============================================================
# Phase 3: Temporal Dynamics Analysis
# GSE28914 + GSE50425 + GSE147890
# ============================================================

suppressMessages({
  library(GEOquery)
  library(limma)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

cat("=== Phase 3: Temporal Dynamics Analysis ===\n")

outdir <- "/home/moog/test/skin/Phase_output/Phase3"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============ 3.1 GSE28914: Acute Wound Healing D0-D1-D3-D7 ============
cat("\n--- 3.1 GSE28914 Acute Healing Timeline ---\n")

gse28 <- getGEO(filename = "/home/moog/test/skin/data/GSE28914/GSE28914_series_matrix.txt.gz",
                getGPL = FALSE)
expr28 <- exprs(gse28)
pdata28 <- pData(gse28)

cat(sprintf("GSE28914: %d probes × %d samples\n", nrow(expr28), ncol(expr28)))
cat("Sample info:\n")
print(pdata28$title)

# Extract timepoints
pdata28$timepoint <- gsub(".*intact skin.*", "D0_intact", pdata28$title)
pdata28$timepoint <- gsub(".*acute wound.*", "D0_wound", pdata28$timepoint)
pdata28$timepoint <- gsub(".*3rd post.*|.*3rd day.*", "D3", pdata28$timepoint)
pdata28$timepoint <- gsub(".*7th post.*|.*7th day.*", "D7", pdata28$timepoint)
pdata28$patient <- gsub("Patient ([0-9]+).*", "\\1", pdata28$title)

cat("Timepoint distribution:\n")
print(table(pdata28$timepoint))

# Annotate probes -> genes
# Check if annotation is needed
cat("Platform:", as.character(pdata28$platform_id[1]), "\n")

# Collapse probes to genes (take maximum value)
# Try using featureData
fdata28 <- fData(gse28)
if (ncol(fdata28) > 0) {
  cat("Feature data columns:", paste(head(colnames(fdata28)), collapse=", "), "\n")
  # Find gene name column
  gene_col <- grep("gene.*symbol|symbol|gene.name", colnames(fdata28), ignore.case = TRUE, value = TRUE)
  if (length(gene_col) > 0) {
    gene_col <- gene_col[1]
    cat("Using gene column:", gene_col, "\n")
  }
}

# ============ 3.1 GSE50425: Acute Healing D0-D14-D21 ============
cat("\n--- GSE50425 ---\n")
gse50 <- getGEO(filename = "/home/moog/test/skin/data/GSE50425/GSE50425_series_matrix.txt.gz",
                getGPL = FALSE)
expr50 <- exprs(gse50)
pdata50 <- pData(gse50)

cat(sprintf("GSE50425: %d probes × %d samples\n", nrow(expr50), ncol(expr50)))
cat("Samples:\n")
print(pdata50$title)

pdata50$timepoint <- gsub(".*intact skin.*", "D0", pdata50$title)
pdata50$timepoint <- gsub(".*14th post.*|.*14th day.*", "D14", pdata50$timepoint)
pdata50$timepoint <- gsub(".*21st post.*|.*21st day.*", "D21", pdata50$timepoint)

cat("Timepoints: ", table(pdata50$timepoint), "\n")

# ============ 3.2 GSE147890: Diabetic Wound Early Response ============
cat("\n--- 3.2 GSE147890 Diabetic Early Response ---\n")
gse147 <- getGEO(filename = "/home/moog/test/skin/data/GSE147890/GSE147890_series_matrix.txt.gz",
                 getGPL = FALSE)
expr147 <- exprs(gse147)
pdata147 <- pData(gse147)

cat(sprintf("GSE147890: %d probes × %d samples\n", nrow(expr147), ncol(expr147)))

pdata147$condition <- gsub(".*Control.*", "Control", pdata147$title)
pdata147$condition <- gsub(".*Diabetic.*", "Diabetic", pdata147$condition)
pdata147$time <- gsub(".*_0h.*", "0h", pdata147$title)
pdata147$time <- gsub(".*_24h.*", "24h", pdata147$time)
pdata147$group <- paste(pdata147$condition, pdata147$time, sep = "_")

cat("Groups:\n")
print(table(pdata147$group))

# Annotate probes -> gene names
# Need to download GPL annotation or use existing one
# Check fData first
fdata147 <- fData(gse147)
cat("GSE147890 fData columns:", paste(head(colnames(fdata147)), collapse=", "), "\n")

# For all three datasets, first use probe IDs to analyze target gene temporal dynamics
# If no fData annotation, use GPL soft file

# Read GPL annotation
cat("\nFetching GPL annotation...\n")

# GSE28914 and GSE50425 use the same platform (HuGene 2.0)
# GSE147890 may use a different platform

# Simplified approach: directly use GPL annotation from soft file
gpl_file <- "/home/moog/test/skin/data/GPL5104.soft.gz"
if (file.exists(gpl_file)) {
  cat("Reading GPL5104 annotation...\n")
  gpl <- getGEO(filename = gpl_file)
  gpl_table <- Table(gpl)
  cat("GPL columns:", paste(head(colnames(gpl_table)), collapse=", "), "\n")
}

# For each dataset, try to map target genes
target_genes <- c("HIF1A", "NFKB1", "RELA", "TNF", "IL1B", "IL6", "VEGFA", "CXCL8")

# GSE28914: check if gene annotation is already available
if (ncol(fdata28) > 0) {
  gene_cols28 <- grep("gene|symbol", colnames(fdata28), ignore.case = TRUE, value = TRUE)
  cat("\nGSE28914 available annotation columns:", paste(gene_cols28, collapse=", "), "\n")
}

# If no annotation, infer from highest-variance probes
# A more reliable approach is to download platform annotation

# Simplified approach: use known probe-gene mapping
# For HuGene 2.0 ST (GPL16686), probe IDs for common genes are known
# Here we use a statistical approach - perform temporal analysis on all probes

# GSE28914 temporal differential analysis
cat("\n--- GSE28914 Temporal Differential Analysis ---\n")

# Only compare D0_intact vs D3 vs D7
use_samples <- pdata28$timepoint %in% c("D0_intact", "D3", "D7")
expr28_sub <- expr28[, use_samples]
time_sub <- pdata28$timepoint[use_samples]

# Check whether data is on log scale
cat("Expression value range:", range(expr28_sub, na.rm=TRUE), "\n")
if (max(expr28_sub, na.rm=TRUE) > 20) {
  cat("Converting to log2...\n")
  expr28_sub <- log2(expr28_sub + 1)
}

# D3 vs D0
design28 <- model.matrix(~ 0 + factor(time_sub))
colnames(design28) <- gsub("factor\\(time_sub\\)", "", colnames(design28))
cont28 <- makeContrasts(D3_vs_D0 = D3 - D0_intact, D7_vs_D0 = D7 - D0_intact, levels = design28)
fit28 <- lmFit(expr28_sub, design28)
fit28_2 <- contrasts.fit(fit28, cont28)
fit28_2 <- eBayes(fit28_2)

deg_d3 <- topTable(fit28_2, coef = "D3_vs_D0", number = Inf)
deg_d7 <- topTable(fit28_2, coef = "D7_vs_D0", number = Inf)

cat(sprintf("D3 vs D0: %d sig probes (FDR<0.05, |logFC|>0.5)\n",
            sum(deg_d3$adj.P.Val < 0.05 & abs(deg_d3$logFC) > 0.5)))
cat(sprintf("D7 vs D0: %d sig probes (FDR<0.05, |logFC|>0.5)\n",
            sum(deg_d7$adj.P.Val < 0.05 & abs(deg_d7$logFC) > 0.5)))

# GSE147890: Control vs Diabetic early response
cat("\n--- GSE147890 Differential Analysis ---\n")

# Check scale
cat("GSE147890 expression range:", range(expr147, na.rm=TRUE), "\n")
if (max(expr147, na.rm=TRUE) > 20) {
  expr147 <- log2(expr147 + 1)
  cat("log2 transformation applied\n")
}

# Four-group comparison: Control_0h, Control_24h, Diabetic_0h, Diabetic_24h
design147 <- model.matrix(~ 0 + factor(pdata147$group))
colnames(design147) <- gsub("factor\\(pdata147\\$group\\)", "", colnames(design147))
colnames(design147) <- gsub("[^A-Za-z0-9_]", "", colnames(design147))

cat("Design matrix column names:", paste(colnames(design147), collapse=", "), "\n")

# Create contrasts: Control response (24h-0h) vs Diabetic response (24h-0h)
cont147 <- makeContrasts(
  Control_response = Control_24h - Control_0h,
  Diabetic_response = Diabetic_24h - Diabetic_0h,
  Interaction = (Diabetic_24h - Diabetic_0h) - (Control_24h - Control_0h),
  levels = design147
)

fit147 <- lmFit(expr147, design147)
fit147_2 <- contrasts.fit(fit147, cont147)
fit147_2 <- eBayes(fit147_2)

deg_ctrl <- topTable(fit147_2, coef = "Control_response", number = Inf)
deg_diab <- topTable(fit147_2, coef = "Diabetic_response", number = Inf)
deg_inter <- topTable(fit147_2, coef = "Interaction", number = Inf)

cat(sprintf("Control 24h response: %d sig probes\n",
            sum(deg_ctrl$adj.P.Val < 0.05 & abs(deg_ctrl$logFC) > 0.5)))
cat(sprintf("Diabetic 24h response: %d sig probes\n",
            sum(deg_diab$adj.P.Val < 0.05 & abs(deg_diab$logFC) > 0.5)))
cat(sprintf("Interaction effect (differential response): %d sig probes\n",
            sum(deg_inter$adj.P.Val < 0.05 & abs(deg_inter$logFC) > 0.5)))

# Save results
write.csv(deg_d3, file.path(outdir, "DEGs_GSE28914_D3_vs_D0.csv"))
write.csv(deg_d7, file.path(outdir, "DEGs_GSE28914_D7_vs_D0.csv"))
write.csv(deg_ctrl, file.path(outdir, "DEGs_GSE147890_Control_response.csv"))
write.csv(deg_diab, file.path(outdir, "DEGs_GSE147890_Diabetic_response.csv"))
write.csv(deg_inter, file.path(outdir, "DEGs_GSE147890_Interaction.csv"))

# ============ 3.3 DFU-Specific Dysregulated Genes ============
cat("\n--- 3.3 DFU-Specific Dysregulated Gene Identification ---\n")

# Acute healing response genes (significant in GSE28914 D3+D7)
acute_sig <- unique(c(
  rownames(deg_d3[deg_d3$adj.P.Val < 0.05 & abs(deg_d3$logFC) > 0.5, ]),
  rownames(deg_d7[deg_d7$adj.P.Val < 0.05 & abs(deg_d7$logFC) > 0.5, ])
))
cat(sprintf("Acute healing response probes: %d\n", length(acute_sig)))

# DFU differential probes (GSE147890 interaction effect)
dfu_sig <- rownames(deg_inter[deg_inter$adj.P.Val < 0.05 & abs(deg_inter$logFC) > 0.5, ])
cat(sprintf("Diabetic differential response probes: %d\n", length(dfu_sig)))

# DFU-specific = diabetic differential - shared with acute healing
dfu_specific <- setdiff(dfu_sig, acute_sig)
cat(sprintf("DFU-specific probes: %d\n", length(dfu_specific)))

# ============ 3.4 Temporal Pattern Summary ============
cat("\n--- 3.4 Temporal Pattern Summary ---\n")
cat("Normal healing: D0->D3/D7 significantly responsive probes:", length(acute_sig), "\n")
cat("Diabetic: probes with attenuated/delayed response (interaction effect):", length(dfu_sig), "\n")
cat("DFU-specific dysregulation (absent in normal healing):", length(dfu_specific), "\n")

cat("\n=== Phase 3 Complete ===\n")
