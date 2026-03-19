# ============================================================
# Phase 2: Bulk RNA-seq Multi-Cohort Validation
# GSE134431 + GSE199939 + GSE80178
# ============================================================

suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readxl)
  library(limma)
  library(WGCNA)
})

cat("=== Phase 2: Bulk RNA-seq Multi-Cohort Validation ===\n")

outdir <- "/home/moog/test/skin/Phase_output/Phase2"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============ 2.1 Data Loading and Preprocessing ============
cat("\n--- 2.1 Data Preprocessing ---\n")

# GSE134431 (RPKM) — 3-row header, sample names in row 2
cat("Loading GSE134431...\n")
rpkm_raw <- read_xlsx("/home/moog/test/skin/data/GSE134431/GSE134431_181016_diabetic_rna-seq_results.gene.rpkm.xlsx",
                      col_names = FALSE)

# Row 2 = sample names, data columns = 2-22 (21 samples), columns 23-25 are annotations
sample_names <- as.character(rpkm_raw[2, 2:22])
rpkm_data <- rpkm_raw[4:nrow(rpkm_raw), ]

gene_names <- as.character(rpkm_data[[1]])
rpkm_mat <- as.data.frame(lapply(rpkm_data[, 2:22], as.numeric))
rownames(rpkm_mat) <- make.unique(gene_names)
colnames(rpkm_mat) <- sample_names

cat(sprintf("GSE134431: %d genes x %d samples\n", nrow(rpkm_mat), ncol(rpkm_mat)))

# Grouping: DFS = contains "DFS", DFU = contains "UM#" or "P2#"
gse134_groups <- data.frame(
  sample = sample_names,
  group = ifelse(grepl("DFS", sample_names), "DFS", "DFU"),
  stringsAsFactors = FALSE
)
cat("GSE134431 groups:\n")
print(table(gse134_groups$group))

# GSE199939 (TPM)
cat("\nLoading GSE199939...\n")
tpm_raw <- read.delim(gzfile("/home/moog/test/skin/data/GSE199939/GSE199939_gene.tpm.matrix.annot_-_21.txt.gz"),
                      header = TRUE, check.names = FALSE)
cat(sprintf("GSE199939 raw: %d rows x %d cols\n", nrow(tpm_raw), ncol(tpm_raw)))
cat("Column names:\n")
print(head(colnames(tpm_raw), 10))

# Extract gene names and expression matrix
gene_col <- which(sapply(tpm_raw, is.character) | sapply(tpm_raw, is.factor))[1]
if (is.na(gene_col)) gene_col <- 1
tpm_genes <- as.character(tpm_raw[[gene_col]])
# Find numeric columns
num_cols <- which(sapply(tpm_raw, is.numeric))
tpm_mat <- tpm_raw[, num_cols]
rownames(tpm_mat) <- make.unique(tpm_genes)

cat(sprintf("GSE199939: %d genes x %d samples\n", nrow(tpm_mat), ncol(tpm_mat)))
cat("Sample names: ", paste(head(colnames(tpm_mat)), collapse=", "), "\n")

# Grouping
gse199_groups <- data.frame(
  sample = colnames(tpm_mat),
  group = ifelse(grepl("^DW|^D[0-9]", colnames(tpm_mat), ignore.case = TRUE) |
                   grepl("diabetic", colnames(tpm_mat), ignore.case = TRUE),
                 "DFU", "Normal"),
  stringsAsFactors = FALSE
)
# Based on data_info: DW1-DW10 are DFU, N1-N11 are Normal
gse199_groups$group <- ifelse(grepl("^DW", colnames(tpm_mat)), "DFU",
                              ifelse(grepl("^N[0-9]", colnames(tpm_mat)), "Normal", "Unknown"))
cat("GSE199939 groups:\n")
print(table(gse199_groups$group))

# GSE80178 (pre-computed data)
cat("\nLoading GSE80178...\n")
gse80 <- read.delim(gzfile("/home/moog/test/skin/data/GSE80178/GSE80178_further_processed_data.txt.gz"),
                    header = TRUE, check.names = FALSE)
cat(sprintf("GSE80178: %d rows x %d cols\n", nrow(gse80), ncol(gse80)))
cat("Column names: ", paste(head(colnames(gse80), 10), collapse=", "), "\n")

# ============ 2.2-2.3 HIF1A/NFKB1 Target Gene Expression ============
cat("\n--- 2.2-2.3 HIF1A/NFKB1 and Target Gene Expression Validation ---\n")

target_genes <- c("HIF1A", "NFKB1", "RELA", "TNF", "IL1B", "IL6", "VEGFA",
                  "CXCL8", "PTGS2", "BNIP3", "ADM", "ANGPTL4", "SLC2A1",
                  "LDHA", "NFKBIA", "ICAM1", "BCL2")

# GSE134431
cat("\nGSE134431 target gene expression:\n")
available_134 <- target_genes[target_genes %in% rownames(rpkm_mat)]
for (g in available_134) {
  dfs_vals <- as.numeric(rpkm_mat[g, gse134_groups$group == "DFS"])
  dfu_vals <- as.numeric(rpkm_mat[g, gse134_groups$group == "DFU"])
  pval <- tryCatch(wilcox.test(dfs_vals, dfu_vals)$p.value, error = function(e) NA)
  fc <- mean(dfu_vals, na.rm = TRUE) / max(mean(dfs_vals, na.rm = TRUE), 0.01)
  cat(sprintf("  %s: DFS=%.2f, DFU=%.2f, FC=%.2f, p=%.2e\n",
              g, mean(dfs_vals, na.rm=T), mean(dfu_vals, na.rm=T), fc, pval))
}

# GSE199939
cat("\nGSE199939 target gene expression:\n")
available_199 <- target_genes[target_genes %in% rownames(tpm_mat)]
for (g in available_199) {
  normal_vals <- as.numeric(tpm_mat[g, gse199_groups$group == "Normal"])
  dfu_vals <- as.numeric(tpm_mat[g, gse199_groups$group == "DFU"])
  pval <- tryCatch(wilcox.test(normal_vals, dfu_vals)$p.value, error = function(e) NA)
  fc <- mean(dfu_vals, na.rm = TRUE) / max(mean(normal_vals, na.rm = TRUE), 0.01)
  cat(sprintf("  %s: Normal=%.2f, DFU=%.2f, FC=%.2f, p=%.2e\n",
              g, mean(normal_vals, na.rm=T), mean(dfu_vals, na.rm=T), fc, pval))
}

# Correlation analysis
cat("\nHIF1A-NFKB1 correlation:\n")
if (all(c("HIF1A", "NFKB1") %in% rownames(rpkm_mat))) {
  r134 <- cor(as.numeric(rpkm_mat["HIF1A",]), as.numeric(rpkm_mat["NFKB1",]),
              method = "spearman")
  cat(sprintf("  GSE134431: rho=%.3f\n", r134))
}
if (all(c("HIF1A", "NFKB1") %in% rownames(tpm_mat))) {
  r199 <- cor(as.numeric(tpm_mat["HIF1A",]), as.numeric(tpm_mat["NFKB1",]),
              method = "spearman")
  cat(sprintf("  GSE199939: rho=%.3f\n", r199))
}

# Target gene expression visualization (GSE199939)
if (length(available_199) > 0) {
  plot_data <- data.frame()
  for (g in available_199[1:min(8, length(available_199))]) {
    df <- data.frame(
      gene = g,
      expression = as.numeric(tpm_mat[g, ]),
      group = gse199_groups$group,
      stringsAsFactors = FALSE
    )
    plot_data <- rbind(plot_data, df)
  }

  p_expr <- ggplot(plot_data, aes(x = group, y = log2(expression + 1), fill = group)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~gene, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = c("DFU" = "#E41A1C", "Normal" = "#4DAF4A")) +
    theme_minimal() + labs(y = "log2(TPM+1)", x = "") +
    ggtitle("Target Gene Expression: DFU vs Normal (GSE199939)")
  ggsave(file.path(outdir, "target_gene_expression_GSE199939.png"), p_expr,
         width = 14, height = 8, dpi = 150)
}

# ============ 2.4 WGCNA (GSE199939) ============
cat("\n--- 2.4 WGCNA Analysis ---\n")

# Filter low-expression genes
tpm_filtered <- tpm_mat[rowMeans(tpm_mat) > 1, ]
cat(sprintf("WGCNA input: %d genes\n", nrow(tpm_filtered)))

# Select top 5000 most variable genes
gene_var <- apply(tpm_filtered, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:5000]
tpm_wgcna <- as.data.frame(t(log2(tpm_filtered[top_genes, ] + 1)))

# Soft threshold selection
cat("Selecting soft threshold...\n")
allowWGCNAThreads(16)
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(tpm_wgcna, powerVector = powers, verbose = 0)
cat("R2 values:\n")
print(data.frame(Power = sft$fitIndices$Power,
                 R2 = round(sft$fitIndices$SFT.R.sq, 3))[1:10, ])

# Select appropriate power
power <- sft$powerEstimate
if (is.na(power)) power <- 6
cat(sprintf("Selected power: %d\n", power))

# Build network
cat("Building WGCNA network...\n")
net <- blockwiseModules(tpm_wgcna, power = power, maxBlockSize = 5000,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        verbose = 0)

module_colors <- labels2colors(net$colors)
cat(sprintf("Number of modules: %d\n", length(unique(module_colors)) - 1))  # subtract grey
cat("Module sizes:\n")
print(sort(table(module_colors), decreasing = TRUE))

# Find modules containing HIF1A/NFKB1
hif_module <- module_colors[names(net$colors) == "HIF1A"]
nfkb_module <- module_colors[names(net$colors) == "NFKB1"]
cat(sprintf("\nHIF1A module: %s\n", ifelse(length(hif_module)>0, hif_module, "not in top5000")))
cat(sprintf("NFKB1 module: %s\n", ifelse(length(nfkb_module)>0, nfkb_module, "not in top5000")))

# Module-trait association
trait <- data.frame(DFU = as.numeric(gse199_groups$group == "DFU"))
rownames(trait) <- rownames(tpm_wgcna)

MEs <- net$MEs
moduleTraitCor <- cor(MEs, trait, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(tpm_wgcna))

cat("\nModule-DFU associations:\n")
for (i in 1:ncol(MEs)) {
  if (moduleTraitPvalue[i, 1] < 0.05) {
    cat(sprintf("  ME%d: r=%.3f, p=%.2e\n", i-1, moduleTraitCor[i,1], moduleTraitPvalue[i,1]))
  }
}

# Save WGCNA results
wgcna_results <- data.frame(
  gene = top_genes,
  module = module_colors,
  stringsAsFactors = FALSE
)
write.csv(wgcna_results, file.path(outdir, "WGCNA_gene_modules.csv"), row.names = FALSE)

# ============ 2.5 Pathway Enrichment (simplified) ============
cat("\n--- 2.5 Pathway Analysis ---\n")

# GSE199939 differential expression analysis
cat("limma differential expression analysis (GSE199939)...\n")
expr_log <- log2(tpm_filtered + 1)
design <- model.matrix(~ 0 + factor(gse199_groups$group))
colnames(design) <- c("DFU", "Normal")
contrast <- makeContrasts(DFU - Normal, levels = design)

fit <- lmFit(expr_log, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
deg_199 <- topTable(fit2, number = Inf)
deg_199$gene <- rownames(deg_199)

cat(sprintf("GSE199939 DEGs: up=%d, down=%d (FDR<0.05, |logFC|>1)\n",
            sum(deg_199$adj.P.Val < 0.05 & deg_199$logFC > 1),
            sum(deg_199$adj.P.Val < 0.05 & deg_199$logFC < -1)))

write.csv(deg_199, file.path(outdir, "DEGs_GSE199939_DFU_vs_Normal.csv"), row.names = FALSE)

# Check HIF/NF-kB pathway genes in DEGs
hif_pathway <- c("HIF1A","VEGFA","SLC2A1","LDHA","PGK1","ENO1","PDK1","BNIP3","ADM","ANGPTL4","LOX","CA9")
nfkb_pathway <- c("NFKB1","RELA","TNF","IL1B","IL6","CXCL8","NFKBIA","ICAM1","PTGS2","BCL2","MMP9")

cat("\nHIF pathway gene changes in DFU:\n")
for (g in hif_pathway) {
  if (g %in% rownames(deg_199)) {
    cat(sprintf("  %s: logFC=%.2f, p=%.2e, sig=%s\n",
                g, deg_199[g, "logFC"], deg_199[g, "adj.P.Val"],
                ifelse(deg_199[g, "adj.P.Val"] < 0.05, "*", "ns")))
  }
}

cat("\nNF-kB pathway gene changes in DFU:\n")
for (g in nfkb_pathway) {
  if (g %in% rownames(deg_199)) {
    cat(sprintf("  %s: logFC=%.2f, p=%.2e, sig=%s\n",
                g, deg_199[g, "logFC"], deg_199[g, "adj.P.Val"],
                ifelse(deg_199[g, "adj.P.Val"] < 0.05, "*", "ns")))
  }
}

# ============ 2.6 Immune Infiltration ============
cat("\n--- 2.6 Immune Infiltration Analysis ---\n")

# Using MCPcounter or simplified gene set scoring approach
# M1/M2 marker gene scores
m1_sig <- c("CD86","NOS2","IL1B","TNF","IL6","CXCL8","CCL2","IRF5","STAT1","PTGS2")
m2_sig <- c("CD163","MRC1","ARG1","IL10","TGFB1","CCL18","STAB1","FOLR2","CD209")

m1_avail <- m1_sig[m1_sig %in% rownames(expr_log)]
m2_avail <- m2_sig[m2_sig %in% rownames(expr_log)]

m1_scores <- colMeans(expr_log[m1_avail, ])
m2_scores <- colMeans(expr_log[m2_avail, ])

immune_df <- data.frame(
  sample = colnames(expr_log),
  group = gse199_groups$group,
  M1_score = m1_scores,
  M2_score = m2_scores,
  M1_M2_ratio = m1_scores / m2_scores
)

cat("M1/M2 scores:\n")
print(immune_df %>% group_by(group) %>%
        summarise(M1=mean(M1_score), M2=mean(M2_score), ratio=mean(M1_M2_ratio), .groups="drop"))

p_immune <- ggplot(immune_df, aes(x = group, fill = group)) +
  geom_boxplot(aes(y = M1_score), alpha = 0.7) +
  theme_minimal() + labs(y = "M1 Score", title = "M1 Polarization Score")
p_immune2 <- ggplot(immune_df, aes(x = group, fill = group)) +
  geom_boxplot(aes(y = M2_score), alpha = 0.7) +
  theme_minimal() + labs(y = "M2 Score", title = "M2 Polarization Score")
ggsave(file.path(outdir, "immune_infiltration.png"), p_immune + p_immune2,
       width = 12, height = 5, dpi = 150)

write.csv(immune_df, file.path(outdir, "immune_infiltration_scores.csv"), row.names = FALSE)

cat("\n=== Phase 2 Complete ===\n")
