#!/usr/bin/env Rscript
# Figures for exTreg vs Stable_Treg: polyfunctional effector phenotype
library(DESeq2)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(dplyr)

project_dir <- "/home/david/projects/luke_mrna"
out_dir <- file.path(project_dir, "deseq2_shrinkage")

# ---- Load data and run DESeq2 with Stable_Treg as reference ----
raw <- read.delim(file.path(project_dir, "counts/gene_counts.txt"),
                  comment.char = "#", check.names = FALSE)
count_matrix <- as.matrix(raw[, 8:ncol(raw)])
rownames(count_matrix) <- raw$Geneid
colnames(count_matrix) <- gsub(".*/aligned/(.*)/Aligned.sortedByCoord.out.bam",
                               "\\1", colnames(count_matrix))
gene_names <- setNames(raw$gene_name, raw$Geneid)

metadata <- read.csv(file.path(project_dir, "sample_metadata.csv"))
rownames(metadata) <- metadata$sample_id
metadata <- metadata[colnames(count_matrix), ]
metadata$cell_type <- factor(metadata$cell_type,
  levels = c("Stable_Treg", "Never_Treg", "exTreg", "Recent_Treg"))
metadata$mouse <- factor(metadata$mouse)

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata,
  design = ~ mouse + cell_type)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)

norm <- counts(dds, normalized = TRUE)
res_shrunk <- lfcShrink(dds, coef = "cell_type_exTreg_vs_Stable_Treg", type = "apeglm")

# Color palette for cell types
ct_colors <- c("Never_Treg" = "#377EB8", "exTreg" = "#E41A1C",
               "Stable_Treg" = "#4DAF4A", "Recent_Treg" = "#FF7F00")
ct_order <- c("Never_Treg", "Recent_Treg", "Stable_Treg", "exTreg")

# Helper: build per-replicate data frame for a gene list
build_expr_df <- function(genes, norm, metadata, gene_names) {
  rows <- list()
  for (g in genes) {
    gid <- names(gene_names)[gene_names == g]
    gid <- gid[gid %in% rownames(norm)]
    if (length(gid) == 0) next
    gid <- gid[1]
    for (s in rownames(metadata)) {
      rows[[length(rows) + 1]] <- data.frame(
        gene = g, sample = s,
        cell_type = as.character(metadata[s, "cell_type"]),
        mouse = as.character(metadata[s, "mouse"]),
        counts = norm[gid, s],
        stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, rows)
  df$cell_type <- factor(df$cell_type, levels = ct_order)
  return(df)
}

# ===========================================================
# Figure 1: Effector cytokine panel — per-replicate dot plot
# ===========================================================
cytokines <- c("Il21", "Il4", "Ifng", "Il2", "Il17a", "Tnf")
cytokine_labels <- c("IL-21", "IL-4", "IFN-g", "IL-2", "IL-17A", "TNF")

df_cyt <- build_expr_df(cytokines, norm, metadata, gene_names)
df_cyt$gene <- factor(df_cyt$gene, levels = cytokines,
                       labels = cytokine_labels)

# Summary stats
df_cyt_sum <- df_cyt %>%
  group_by(gene, cell_type) %>%
  summarize(mean_counts = mean(counts), se = sd(counts) / sqrt(n()), .groups = "drop")

p1 <- ggplot(df_cyt_sum, aes(x = cell_type, y = mean_counts, fill = cell_type)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.3) +
  geom_jitter(data = df_cyt, aes(x = cell_type, y = counts, fill = cell_type),
              shape = 21, size = 2.5, width = 0.15, color = "black", stroke = 0.5) +
  geom_errorbar(aes(ymin = mean_counts - se, ymax = mean_counts + se),
                width = 0.2, linewidth = 0.4) +
  facet_wrap(~ gene, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = ct_colors) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank()) +
  labs(x = NULL, y = "Normalized counts",
       title = "Effector Cytokine Expression Across CD4 T Cell Subsets")

ggsave(file.path(out_dir, "cytokine_panel_barplot.pdf"), p1,
       width = 10, height = 7)
cat("Saved cytokine_panel_barplot.pdf\n")

# ===========================================================
# Figure 2: Heatmap of key genes — exTreg vs Stable_Treg focus
# ===========================================================
heatmap_genes <- c(
  # Effector cytokines
  "Il21", "Il4", "Ifng", "Il2", "Tnf",
  # Suppressive
  "Il10", "Fgl2", "Tgfb1",
  # Tfh/GC
  "Maf", "Bcl6", "Cxcr5", "Sh2d1a", "Slamf6", "Tox2", "Pou2af1",
  # GC-associated
  "Pax5", "Sostdc1", "Serpina9",
  # Th17
  "Rorc", "Ccr6", "Bhlhe40", "Hif1a",
  # Exhaustion
  "Pdcd1", "Lag3", "Tigit",
  # Treg identity
  "Foxp3", "Il2ra", "Ikzf2", "Nrp1",
  # Regulators
  "Bach2"
)

# Map gene names to IDs
hm_ids <- sapply(heatmap_genes, function(g) {
  gid <- names(gene_names)[gene_names == g]
  gid <- gid[gid %in% rownames(norm)]
  if (length(gid) > 0) gid[1] else NA
})
hm_ids <- hm_ids[!is.na(hm_ids)]

# VST for heatmap
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[hm_ids, ]
mat <- mat - rowMeans(mat)  # center
rownames(mat) <- gene_names[rownames(mat)]

# Reorder columns: group by cell type
col_order <- rownames(metadata)[order(match(metadata$cell_type, ct_order))]
mat <- mat[, col_order]

# Annotation
annot_col <- metadata[col_order, c("cell_type", "mouse"), drop = FALSE]
annot_colors <- list(
  cell_type = ct_colors,
  mouse = c("mouse1" = "#FDBF6F", "mouse2" = "#A6CEE3", "mouse3" = "#B2DF8A"))

# Row annotation by functional category
row_annot <- data.frame(
  Category = c(rep("Effector cytokine", 5),
               rep("Suppressive", 3),
               rep("Tfh/GC program", 7),
               rep("GC-associated", 3),
               rep("Th17", 4),
               rep("Inhibitory receptor", 3),
               rep("Treg identity", 4),
               rep("Regulator", 1)),
  row.names = gene_names[hm_ids])

cat_colors <- c("Effector cytokine" = "#E41A1C", "Suppressive" = "#4DAF4A",
                "Tfh/GC program" = "#984EA3", "GC-associated" = "#FF7F00",
                "Th17" = "#A65628", "Inhibitory receptor" = "#F781BF",
                "Treg identity" = "#999999", "Regulator" = "#377EB8")
annot_colors$Category <- cat_colors

pdf(file.path(out_dir, "heatmap_extreg_stable_key_genes.pdf"),
    width = 10, height = 12)
pheatmap(mat,
  annotation_col = annot_col,
  annotation_row = row_annot,
  annotation_colors = annot_colors,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 8,
  gaps_col = cumsum(table(metadata$cell_type[match(col_order, rownames(metadata))])),
  gaps_row = c(5, 8, 15, 18, 22, 25, 29),
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(-4, 4, length.out = 101),
  main = "Key Gene Expression Across CD4 T Cell Subsets (VST, centered)")
dev.off()
cat("Saved heatmap_extreg_stable_key_genes.pdf\n")

# ===========================================================
# Figure 3: Volcano plot — exTreg vs Stable_Treg (shrunken)
# ===========================================================
res_df <- as.data.frame(res_shrunk)
res_df$gene_name <- gene_names[rownames(res_df)]
res_df$significant <- "NS"
res_df$significant[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up in exTreg"
res_df$significant[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Up in Stable_Treg"
res_df$significant <- factor(res_df$significant,
  levels = c("Up in Stable_Treg", "NS", "Up in exTreg"))

# Label genes of interest
label_genes <- c("Il21", "Il4", "Ifng", "Il2", "Il17a", "Tnf",
                 "Foxp3", "Pax5", "Sostdc1", "Serpina9",
                 "Maf", "Bcl6", "Cxcr5", "Pdcd1", "Tigit",
                 "Pou2af1", "Slamf6", "Tox2", "Sh2d1a",
                 "Rorc", "Bhlhe40", "Bach2", "Fgl2", "Il10",
                 "Il2ra", "Nrp1", "Ikzf2", "Hif1a")
res_df$label <- ifelse(res_df$gene_name %in% label_genes, res_df$gene_name, NA)

p3 <- ggplot(res_df, aes(log2FoldChange, -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_color_manual(values = c("Up in Stable_Treg" = "#4DAF4A",
                                "NS" = "grey75",
                                "Up in exTreg" = "#E41A1C"),
                     name = NULL) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(aes(label = label), size = 3.2, max.overlaps = 30,
                  show.legend = FALSE, color = "black",
                  min.segment.length = 0.1, seed = 42) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.15, 0.9),
        legend.background = element_rect(fill = alpha("white", 0.8))) +
  labs(x = "Shrunken log2 Fold Change (apeglm)",
       y = "-log10(p-value)",
       title = "exTreg vs Stable_Treg (apeglm shrunken LFC)")

ggsave(file.path(out_dir, "volcano_exTreg_vs_Stable_Treg_annotated.pdf"),
       p3, width = 10, height = 8)
cat("Saved volcano_exTreg_vs_Stable_Treg_annotated.pdf\n")

# ===========================================================
# Figure 4: Polyfunctional cytokine dot plot — exTreg vs Stable only
# ===========================================================
robust_cytokines <- c("Il21", "Il4", "Ifng", "Il2")
robust_labels <- c("IL-21\n(Tfh)", "IL-4\n(Th2/Tfh)", "IFN-g\n(Th1)", "IL-2\n(activation)")

df_robust <- build_expr_df(robust_cytokines, norm, metadata, gene_names)
df_robust <- df_robust[df_robust$cell_type %in% c("exTreg", "Stable_Treg"), ]
df_robust$gene <- factor(df_robust$gene, levels = robust_cytokines,
                          labels = robust_labels)
df_robust$cell_type <- droplevels(df_robust$cell_type)

df_robust_sum <- df_robust %>%
  group_by(gene, cell_type) %>%
  summarize(mean_counts = mean(counts), se = sd(counts) / sqrt(n()), .groups = "drop")

# Add significance annotations
sig_df <- data.frame(
  gene = factor(robust_labels, levels = robust_labels),
  label = c("p=2.0e-18", "p=4.7e-06", "p=1.2e-03", "p=1.3e-03"),
  y = c(9000, 550, 1700, 1100)
)

p4 <- ggplot(df_robust_sum, aes(x = cell_type, y = mean_counts, fill = cell_type)) +
  geom_col(width = 0.65, color = "grey30", linewidth = 0.3) +
  geom_jitter(data = df_robust, aes(x = cell_type, y = counts),
              shape = 21, size = 3, width = 0.1, color = "black",
              stroke = 0.5) +
  geom_errorbar(aes(ymin = pmax(mean_counts - se, 0), ymax = mean_counts + se),
                width = 0.15, linewidth = 0.4) +
  geom_text(data = sig_df, aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE, size = 3.2, fontface = "italic") +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("exTreg" = "#E41A1C", "Stable_Treg" = "#4DAF4A")) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = NULL, y = "Normalized counts",
       title = "Robust Effector Cytokines: exTreg vs Stable_Treg",
       subtitle = "apeglm-shrunken padj values shown; individual replicates plotted")

ggsave(file.path(out_dir, "cytokine_extreg_vs_stable.pdf"), p4,
       width = 10, height = 5)
cat("Saved cytokine_extreg_vs_stable.pdf\n")

# ===========================================================
# Figure 5: Gained vs retained vs lost gene summary
# ===========================================================
# Barplot showing categories of genes in exTreg relative to Stable_Treg
categories <- data.frame(
  gene = c("IL-21", "IL-4", "IFN-g", "IL-2", "Pou2af1", "Bhlhe40", "Tox2", "PD-1",
           "CXCR5", "TIGIT", "LAG-3", "c-MAF", "Slamf6", "Bcl6",
           "Foxp3", "Fgl2", "BACH2", "BLIMP-1"),
  lfc = c(15.75, 11.67, 4.73, 3.88, 4.57, 2.74, 2.97, 2.37,
          0.42, 0.32, 0.57, 0.98, 1.46, 0.77,
          -8.49, -4.90, -1.50, -1.57),
  category = c(rep("Gained effector", 8),
               rep("Shared (Tfr origin)", 6),
               rep("Lost suppressive", 4)),
  stringsAsFactors = FALSE
)
categories$category <- factor(categories$category,
  levels = c("Gained effector", "Shared (Tfr origin)", "Lost suppressive"))
categories$gene <- factor(categories$gene, levels = rev(categories$gene))

p5 <- ggplot(categories, aes(x = lfc, y = gene, fill = category)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Gained effector" = "#E41A1C",
                               "Shared (Tfr origin)" = "#984EA3",
                               "Lost suppressive" = "#4DAF4A"),
                    name = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.78, 0.15),
        legend.background = element_rect(fill = alpha("white", 0.9)),
        panel.grid.major.y = element_blank()) +
  labs(x = "Shrunken log2 Fold Change (exTreg / Stable_Treg)",
       y = NULL,
       title = "exTreg Transcriptional Changes Relative to Stable_Treg",
       subtitle = "Dashed lines = |LFC| = 1 threshold")

ggsave(file.path(out_dir, "gained_retained_lost_barplot.pdf"), p5,
       width = 9, height = 7)
cat("Saved gained_retained_lost_barplot.pdf\n")

cat("\n=== All figures saved to deseq2_shrinkage/ ===\n")
