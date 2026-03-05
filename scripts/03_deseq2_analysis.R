#!/usr/bin/env Rscript
# DESeq2 differential expression analysis for Luke's CD4 T cell RNA-seq
# 4 populations: Never_Treg, exTreg, Stable_Treg, Recent_Treg
# Paired design: 3 mice (blocking factor)

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(UpSetR)
library(ggrepel)

# ---- Setup ----
project_dir <- "/home/david/projects/luke_mrna"
counts_file <- file.path(project_dir, "counts/gene_counts.txt")
metadata_file <- file.path(project_dir, "sample_metadata.csv")
out_dir <- file.path(project_dir, "deseq2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load data ----
# featureCounts output: first 6 cols are annotation, then counts
raw <- read.delim(counts_file, comment.char = "#", check.names = FALSE)
count_matrix <- as.matrix(raw[, 8:ncol(raw)])
rownames(count_matrix) <- raw$Geneid

# Clean column names (featureCounts uses full BAM paths)
colnames(count_matrix) <- gsub(".*/aligned/(.*)/Aligned.sortedByCoord.out.bam", "\\1", colnames(count_matrix))

# Load metadata
metadata <- read.csv(metadata_file)
rownames(metadata) <- metadata$sample_id

# Ensure metadata order matches count matrix columns
metadata <- metadata[colnames(count_matrix), ]
stopifnot(all(rownames(metadata) == colnames(count_matrix)))

# Set factor levels (Never_Treg as reference)
metadata$cell_type <- factor(metadata$cell_type,
  levels = c("Never_Treg", "exTreg", "Stable_Treg", "Recent_Treg"))
metadata$mouse <- factor(metadata$mouse)

cat("Sample counts per group:\n")
print(table(metadata$cell_type))

# ---- DESeq2 with paired design ----
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ mouse + cell_type  # mouse as blocking factor
)

# Pre-filter low-count genes
keep <- rowSums(counts(dds) >= 10) >= 3  # at least 3 samples with >= 10 counts
dds <- dds[keep, ]
cat(sprintf("\nGenes after filtering: %d\n", nrow(dds)))

# Gene ID to name mapping
gene_names <- setNames(raw$gene_name, raw$Geneid)

# Run DESeq2
dds <- DESeq(dds)

# ---- QC Plots ----
# Variance-stabilized transform for visualization
vsd <- vst(dds, blind = FALSE)

# PCA plot
pca_data <- plotPCA(vsd, intgroup = c("cell_type", "mouse"), returnData = TRUE)
pct_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = cell_type, shape = mouse)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 14) +
  ggtitle("PCA: CD4 T cell populations")
ggsave(file.path(out_dir, "pca_plot.pdf"), p_pca, width = 8, height = 6)

# Sample distance heatmap
sample_dists <- dist(t(assay(vsd)))
dist_matrix <- as.matrix(sample_dists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pdf(file.path(out_dir, "sample_distance_heatmap.pdf"), width = 10, height = 8)
pheatmap(dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  col = colors,
  annotation_col = metadata[, c("cell_type", "mouse"), drop = FALSE],
  main = "Sample-to-sample distances")
dev.off()

# ---- Pairwise comparisons ----
cell_types <- levels(metadata$cell_type)
comparisons <- combn(cell_types, 2, simplify = FALSE)

all_results <- list()
deg_lists <- list()

for (comp in comparisons) {
  contrast_name <- paste0(comp[2], "_vs_", comp[1])
  cat(sprintf("\nContrast: %s\n", contrast_name))

  res <- results(dds, contrast = c("cell_type", comp[2], comp[1]),
                 alpha = 0.05)
  res <- res[order(res$padj), ]

  # Summary
  cat(sprintf("  Total tested: %d\n", sum(!is.na(res$padj))))
  cat(sprintf("  Up (padj<0.05, LFC>1):   %d\n",
    sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE)))
  cat(sprintf("  Down (padj<0.05, LFC<-1): %d\n",
    sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE)))

  # Save results with gene names
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df$gene_name <- gene_names[rownames(res_df)]
  write.csv(res_df, file.path(out_dir, paste0("DEG_", contrast_name, ".csv")),
            row.names = FALSE)

  all_results[[contrast_name]] <- res

  # DEG list for UpSet plot (padj < 0.05, |LFC| > 1)
  sig <- rownames(res)[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)]
  deg_lists[[contrast_name]] <- sig

  # Volcano plot with gene name labels
  res_plot <- as.data.frame(res)
  res_plot$gene_name <- gene_names[rownames(res_plot)]
  res_plot$significant <- "NS"
  res_plot$significant[res_plot$padj < 0.05 & res_plot$log2FoldChange > 1] <- "Up"
  res_plot$significant[res_plot$padj < 0.05 & res_plot$log2FoldChange < -1] <- "Down"
  res_plot$significant <- factor(res_plot$significant, levels = c("Down", "NS", "Up"))

  # Label top 15 DEGs by adjusted p-value (significant only)
  res_plot$label <- NA
  sig_rows <- which(res_plot$significant != "NS")
  if (length(sig_rows) > 0) {
    top_idx <- sig_rows[order(res_plot$padj[sig_rows])][1:min(15, length(sig_rows))]
    res_plot$label[top_idx] <- res_plot$gene_name[top_idx]
  }

  p_vol <- ggplot(res_plot, aes(log2FoldChange, -log10(pvalue), color = significant)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("Down" = "blue", "NS" = "grey70", "Up" = "red")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 20,
                    show.legend = FALSE, color = "black") +
    theme_bw(base_size = 12) +
    ggtitle(contrast_name) +
    xlab("log2 Fold Change") +
    ylab("-log10(p-value)")
  ggsave(file.path(out_dir, paste0("volcano_", contrast_name, ".pdf")),
         p_vol, width = 8, height = 6)
}

# ---- UpSet plot of DEGs across comparisons ----
# Only include comparisons that have DEGs
deg_lists_nonempty <- deg_lists[sapply(deg_lists, length) > 0]

if (length(deg_lists_nonempty) >= 2) {
  pdf(file.path(out_dir, "upset_DEGs.pdf"), width = 12, height = 8)
  print(upset(fromList(deg_lists_nonempty),
    nsets = length(deg_lists_nonempty),
    order.by = "freq",
    text.scale = 1.3,
    mainbar.y.label = "Shared DEGs",
    sets.x.label = "DEGs per comparison"))
  dev.off()
  cat("\nUpSet plot saved.\n")
} else {
  cat("\nNot enough comparisons with DEGs for UpSet plot.\n")
}

# ---- Top genes heatmap ----
# Union of top 50 DEGs from each comparison
top_genes <- unique(unlist(lapply(all_results, function(res) {
  head(rownames(res[order(res$padj), ]), 50)
})))

top_labels <- gene_names[top_genes]

mat <- assay(vsd)[top_genes[top_genes %in% rownames(vsd)], ]
mat <- mat - rowMeans(mat)  # center

pdf(file.path(out_dir, "heatmap_top_DEGs.pdf"), width = 10, height = 14)
pheatmap(mat,
  annotation_col = metadata[, c("cell_type", "mouse"), drop = FALSE],
  show_rownames = TRUE,
  labels_row = gene_names[rownames(mat)],
  fontsize_row = 6,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  main = "Top DEGs across all comparisons (VST, centered)")
dev.off()

# ---- Foxp3 and key gene check ----
# Sanity check: Foxp3 should be high in Stable_Treg and Recent_Treg
key_genes <- c("Foxp3", "Pdcd1", "Ctla4", "Tigit", "Lag3", "Il2ra", "Ikzf2")
cat("\n=== Key gene normalized counts (sanity check) ===\n")
norm_counts <- counts(dds, normalized = TRUE)
for (gene_name in key_genes) {
  gene_ids <- names(gene_names)[gene_names == gene_name]
  gene_ids <- gene_ids[gene_ids %in% rownames(norm_counts)]
  if (length(gene_ids) > 0) {
    cat(sprintf("\n%s (%s):\n", gene_name, gene_ids[1]))
    vals <- norm_counts[gene_ids[1], ]
    for (ct in levels(metadata$cell_type)) {
      samples <- metadata$sample_id[metadata$cell_type == ct]
      cat(sprintf("  %-12s mean=%.1f\n", ct, mean(vals[samples])))
    }
  }
}

# Save normalized counts matrix
norm_df <- as.data.frame(norm_counts)
norm_df$gene_id <- rownames(norm_df)
norm_df$gene_name <- gene_names[rownames(norm_df)]
write.csv(norm_df, file.path(out_dir, "normalized_counts.csv"), row.names = FALSE)

cat(sprintf("\n=== Analysis complete ===\nResults in: %s\n", out_dir))
