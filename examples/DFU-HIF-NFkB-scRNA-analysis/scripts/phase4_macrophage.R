# ============================================================
# Phase 4: In-depth analysis of macrophage HIF/NF-kB mechanisms
# GSE15949 + GSE16099 + GSE165816 macrophages
# ============================================================

suppressMessages({
  library(GEOquery)
  library(limma)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

cat("=== Phase 4: Macrophage HIF/NF-kB Mechanism Analysis ===\n")

outdir <- "/home/moog/test/skin/Phase_output/Phase4"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============ 4.1 Hypoxia response gene signature ============
cat("\n--- 4.1 Macrophage hypoxia response gene signature ---\n")

# GSE15949
cat("Reading GSE15949...\n")
gse15 <- getGEO(filename = "/home/moog/test/skin/data/GSE15949/GSE15949_series_matrix.txt.gz",
                getGPL = FALSE)
expr15 <- exprs(gse15)
pdata15 <- pData(gse15)
cat(sprintf("GSE15949: %d probes × %d samples\n", nrow(expr15), ncol(expr15)))
cat("Samples:", pdata15$title, "\n")

# Detect scale
if (max(expr15, na.rm=TRUE) > 20) {
  cat("Linear scale, applying log2 transformation\n")
  expr15 <- log2(pmax(expr15, 1))
}

# GSE16099
cat("\nReading GSE16099...\n")
gse16 <- getGEO(filename = "/home/moog/test/skin/data/GSE16099/GSE16099_series_matrix.txt.gz",
                getGPL = FALSE)
expr16 <- exprs(gse16)
pdata16 <- pData(gse16)
cat(sprintf("GSE16099: %d probes × %d samples\n", nrow(expr16), ncol(expr16)))
cat("Samples:", pdata16$title, "\n")

if (max(expr16, na.rm=TRUE) > 20) {
  cat("Linear scale, applying log2 transformation\n")
  expr16 <- log2(pmax(expr16, 1))
}

# GSE15949 differential analysis: normoxia vs hypoxia (only 2 samples, cannot use limma)
cat("\nGSE15949: Direct comparison hypoxia vs normoxia\n")
fc15 <- expr15[, 2] - expr15[, 1]  # hypoxia - normoxia (assuming column 2 is hypoxia)
# Confirm based on title
if (grepl("hypoxia", pdata15$title[1], ignore.case = TRUE)) {
  fc15 <- expr15[, 1] - expr15[, 2]
}
cat(sprintf("Upregulated probes (logFC>1): %d, Downregulated: %d\n",
            sum(fc15 > 1, na.rm=TRUE), sum(fc15 < -1, na.rm=TRUE)))

# GSE16099 differential analysis (2 normoxia vs 2 hypoxia)
cat("\nGSE16099 limma analysis...\n")
groups16 <- ifelse(grepl("hypoxia", pdata16$title, ignore.case = TRUE), "Hypoxia", "Normoxia")
cat("Groups:", table(groups16), "\n")

design16 <- model.matrix(~ 0 + factor(groups16))
colnames(design16) <- c("Hypoxia", "Normoxia")
contrast16 <- makeContrasts(Hypoxia - Normoxia, levels = design16)
fit16 <- lmFit(expr16, design16)
fit16_2 <- contrasts.fit(fit16, contrast16)
fit16_2 <- eBayes(fit16_2)
deg16 <- topTable(fit16_2, number = Inf)

cat(sprintf("GSE16099 differential probes: up=%d, down=%d (p<0.05, |logFC|>1)\n",
            sum(deg16$P.Value < 0.05 & deg16$logFC > 1),
            sum(deg16$P.Value < 0.05 & deg16$logFC < -1)))

# Read GPL annotation to map probes to genes
gpl <- getGEO(filename = "/home/moog/test/skin/data/GPL5104.soft.gz")
gpl_table <- Table(gpl)
cat("GPL5104 annotation rows:", nrow(gpl_table), "\n")

# Map GSE16099 probes to genes
probe2gene <- gpl_table[, c("ID", "SYMBOL")]
probe2gene <- probe2gene[probe2gene$SYMBOL != "" & !is.na(probe2gene$SYMBOL), ]

# Merge
deg16$probe <- rownames(deg16)
deg16_annot <- merge(deg16, probe2gene, by.x = "probe", by.y = "ID", all.x = FALSE)
cat(sprintf("After annotation: %d / %d probes mapped to genes\n", nrow(deg16_annot), nrow(deg16)))

# Hypoxia upregulated genes
hypoxia_up <- deg16_annot %>%
  filter(P.Value < 0.05 & logFC > 0.5) %>%
  arrange(desc(logFC))

cat(sprintf("\nHypoxia upregulated genes: %d\n", nrow(hypoxia_up)))
cat("Top20:\n")
print(head(hypoxia_up[, c("SYMBOL", "logFC", "P.Value")], 20))

# Build hypoxia signature gene set
hypoxia_sig_genes <- unique(hypoxia_up$SYMBOL[1:min(100, nrow(hypoxia_up))])
cat(sprintf("Hypoxia signature gene count: %d\n", length(hypoxia_sig_genes)))

# ============ 4.2 Common target gene network ============
cat("\n--- 4.2 HIF/NF-kB common target genes ---\n")

# Literature-known HIF target genes
hif_targets <- c("VEGFA", "SLC2A1", "LDHA", "PGK1", "ENO1", "ALDOA", "PKM",
                 "HK2", "PDK1", "BNIP3", "BNIP3L", "CA9", "ADM", "ANGPTL4",
                 "LOX", "P4HA1", "P4HA2", "EGLN1", "EGLN3", "EPO", "HMOX1",
                 "NOS2", "SLC2A3", "PFKFB3", "GPI", "TPI1")

# Literature-known NF-kB target genes
nfkb_targets <- c("TNF", "IL1B", "IL6", "CXCL8", "CCL2", "CCL5", "ICAM1",
                  "VCAM1", "PTGS2", "MMP9", "NFKBIA", "NFKBIB", "BCL2",
                  "BIRC3", "TRAF1", "CXCL10", "IL2", "IFNG", "CSF2",
                  "LTA", "CD80", "CD86", "SOD2", "NOS2")

# Common target genes (literature intersection)
common_lit <- intersect(hif_targets, nfkb_targets)
cat("Literature-known common target genes:", paste(common_lit, collapse=", "), "\n")

# Validate in hypoxia response data
hypoxia_validated <- hypoxia_up$SYMBOL[hypoxia_up$SYMBOL %in% c(hif_targets, nfkb_targets)]
cat("\nHIF/NF-kB target genes validated in hypoxia response:\n")
for (g in hypoxia_validated) {
  row <- hypoxia_up[hypoxia_up$SYMBOL == g, ][1, ]
  is_hif <- g %in% hif_targets
  is_nfkb <- g %in% nfkb_targets
  cat(sprintf("  %s: logFC=%.2f, HIF_target=%s, NF-kB_target=%s, Co-regulated=%s\n",
              g, row$logFC, is_hif, is_nfkb, is_hif & is_nfkb))
}

# ============ 4.3 Apply hypoxia signature to DFU single cells ============
cat("\n--- 4.3 Applying hypoxia signature to DFU single cells ---\n")

mac <- readRDS("/home/moog/test/skin/Phase_output/Phase1/macrophage/macrophage_subclustered.rds")
cat(sprintf("Macrophage count: %d\n", ncol(mac)))

# Compute hypoxia signature score
hypoxia_in_sc <- hypoxia_sig_genes[hypoxia_sig_genes %in% rownames(mac)]
cat(sprintf("Hypoxia signature genes in scRNA: %d/%d\n", length(hypoxia_in_sc), length(hypoxia_sig_genes)))

mac <- AddModuleScore(mac, features = list(hypoxia_in_sc), name = "Hypoxia_sig")

# Score comparison across conditions
hyp_data <- FetchData(mac, vars = c("Hypoxia_sig1", "condition", "seurat_clusters"))

hyp_summary <- hyp_data %>%
  group_by(condition) %>%
  summarise(
    mean_score = mean(Hypoxia_sig1),
    sd_score = sd(Hypoxia_sig1),
    n = n(),
    .groups = "drop"
  )
cat("\nHypoxia signature score (macrophages):\n")
print(as.data.frame(hyp_summary))

# Statistical test
healing_scores <- hyp_data$Hypoxia_sig1[hyp_data$condition == "Healing"]
nonhealing_scores <- hyp_data$Hypoxia_sig1[hyp_data$condition == "NonHealing"]
wt <- wilcox.test(healing_scores, nonhealing_scores)
cat(sprintf("\nHealing vs NonHealing: p=%.2e\n", wt$p.value))

# FeaturePlot
p_hyp <- FeaturePlot(mac, features = "Hypoxia_sig1", reduction = "umap", pt.size = 0.3) +
  scale_color_viridis_c() + ggtitle("Macrophage Hypoxia Signature Score")
ggsave(file.path(outdir, "macrophage_hypoxia_signature.png"), p_hyp,
       width = 8, height = 6, dpi = 150)

# VlnPlot by condition
p_vln <- VlnPlot(mac, features = "Hypoxia_sig1", group.by = "condition", pt.size = 0,
                 cols = c("Healthy"="#4DAF4A","DM"="#377EB8",
                          "Healing"="#FF7F00","NonHealing"="#E41A1C")) +
  ggtitle("Hypoxia Signature in Macrophages")
ggsave(file.path(outdir, "hypoxia_signature_violin.png"), p_vln,
       width = 8, height = 5, dpi = 150)

# ============ 4.4-4.5 Macrophage subsets and polarization ============
cat("\n--- 4.4-4.5 Macrophage subsets and polarization coupling ---\n")

# M1/M2/HIF/NFKB scores already available, analyze subset-condition associations
mac_meta <- mac@meta.data

# Distribution of each subset across conditions
cluster_cond <- mac_meta %>%
  group_by(seurat_clusters, condition) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(pct = n / sum(n) * 100)

# Find Healing-enriched subsets
healing_enriched <- cluster_cond %>%
  filter(condition == "Healing") %>%
  arrange(desc(pct))
cat("Healing-enriched macrophage subsets:\n")
print(as.data.frame(healing_enriched[1:5, ]))

# Comprehensive assessment: subsets with high HIF+NF-kB and high M1
cluster_profile <- mac_meta %>%
  group_by(seurat_clusters) %>%
  summarise(
    n = n(),
    HIF = mean(HIF_score1),
    NFKB = mean(NFKB_score1),
    M1 = mean(M1_score1),
    M2 = mean(M2_score1),
    Hypoxia = mean(Hypoxia_sig1),
    Healing_pct = sum(condition == "Healing") / n() * 100,
    NonHealing_pct = sum(condition == "NonHealing") / n() * 100,
    .groups = "drop"
  ) %>%
  mutate(
    HIF_NF_coupled = HIF + NFKB,
    M1_M2_bias = M1 - M2
  ) %>%
  arrange(desc(HIF_NF_coupled))

cat("\nMacrophage subset comprehensive assessment:\n")
print(as.data.frame(cluster_profile[, c("seurat_clusters", "n", "HIF", "NFKB",
                                         "M1", "M2", "Healing_pct", "NonHealing_pct")]))

write.csv(cluster_profile, file.path(outdir, "macrophage_cluster_profile.csv"), row.names = FALSE)

# Polarization-pathway coupling scatter plot
p_coupling <- ggplot(mac_meta, aes(x = HIF_score1 + NFKB_score1, y = M1_score1 - M2_score1,
                                   color = condition)) +
  geom_point(size = 0.3, alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Healthy"="#4DAF4A","DM"="#377EB8",
                                "Healing"="#FF7F00","NonHealing"="#E41A1C")) +
  theme_minimal() + labs(x = "HIF+NF-kB Score", y = "M1-M2 Bias") +
  ggtitle("HIF/NF-kB Coupling with M1/M2 Polarization")
ggsave(file.path(outdir, "polarization_coupling.png"), p_coupling,
       width = 8, height = 6, dpi = 150)

# Global correlation
cat("\nPolarization-pathway global correlation:\n")
cat(sprintf("  (HIF+NFKB) vs (M1-M2): rho=%.3f\n",
            cor(mac_meta$HIF_score1 + mac_meta$NFKB_score1,
                mac_meta$M1_score1 - mac_meta$M2_score1, method = "spearman")))

# Save
write.csv(hypoxia_sig_genes, file.path(outdir, "hypoxia_signature_genes.csv"), row.names = FALSE)

cat("\n=== Phase 4 Complete ===\n")
