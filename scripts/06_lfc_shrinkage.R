#!/usr/bin/env Rscript
# LFC shrinkage analysis with baseMean filtering and MA plots
# Companion to 03_deseq2_analysis.R — uses apeglm shrinkage to tame
# inflated fold changes from low-count genes

library(DESeq2)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

# ---- Setup ----
project_dir <- "/home/david/projects/luke_mrna"
counts_file <- file.path(project_dir, "counts/gene_counts.txt")
metadata_file <- file.path(project_dir, "sample_metadata.csv")
out_dir <- file.path(project_dir, "deseq2_shrinkage")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

basemean_cutoff <- 50  # minimum average expression across all samples

# ---- Load data (same as 03) ----
raw <- read.delim(counts_file, comment.char = "#", check.names = FALSE)
count_matrix <- as.matrix(raw[, 8:ncol(raw)])
rownames(count_matrix) <- raw$Geneid
colnames(count_matrix) <- gsub(".*/aligned/(.*)/Aligned.sortedByCoord.out.bam",
                               "\\1", colnames(count_matrix))

metadata <- read.csv(metadata_file)
rownames(metadata) <- metadata$sample_id

# Ensure order matches
metadata <- metadata[colnames(count_matrix), ]
stopifnot(all(rownames(metadata) == colnames(count_matrix)))

metadata$cell_type <- factor(metadata$cell_type,
  levels = c("Never_Treg", "exTreg", "Stable_Treg", "Recent_Treg"))
metadata$mouse <- factor(metadata$mouse)

# Gene ID to name mapping
gene_names <- setNames(raw$gene_name, raw$Geneid)

# ---- Helper: run DESeq2 with a given reference level ----
run_deseq <- function(count_matrix, metadata, ref_level) {
  md <- metadata
  md$cell_type <- relevel(md$cell_type, ref = ref_level)
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = md,
    design = ~ mouse + cell_type
  )
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  dds <- DESeq(dds)
  return(dds)
}

# ---- Define all 6 pairwise comparisons ----
# apeglm requires coefficient names, which are relative to the reference level.
# Group comparisons by which reference level they need.
comparisons <- list(
  list(ref = "Never_Treg", coef = "cell_type_exTreg_vs_Never_Treg",
       name = "exTreg_vs_Never_Treg"),
  list(ref = "Never_Treg", coef = "cell_type_Stable_Treg_vs_Never_Treg",
       name = "Stable_Treg_vs_Never_Treg"),
  list(ref = "Never_Treg", coef = "cell_type_Recent_Treg_vs_Never_Treg",
       name = "Recent_Treg_vs_Never_Treg"),
  list(ref = "exTreg", coef = "cell_type_Stable_Treg_vs_exTreg",
       name = "Stable_Treg_vs_exTreg"),
  list(ref = "exTreg", coef = "cell_type_Recent_Treg_vs_exTreg",
       name = "Recent_Treg_vs_exTreg"),
  list(ref = "Stable_Treg", coef = "cell_type_Recent_Treg_vs_Stable_Treg",
       name = "Recent_Treg_vs_Stable_Treg")
)

# ---- Fit models for each unique reference level ----
ref_levels <- unique(sapply(comparisons, `[[`, "ref"))
dds_models <- list()

for (ref in ref_levels) {
  cat(sprintf("Fitting DESeq2 with reference: %s\n", ref))
  dds_models[[ref]] <- run_deseq(count_matrix, metadata, ref)
}

# ---- Run shrinkage for each comparison ----
all_shrunk <- list()
all_unshrunk <- list()
summary_rows <- list()

for (comp in comparisons) {
  cat(sprintf("\n=== %s (apeglm shrinkage) ===\n", comp$name))

  dds <- dds_models[[comp$ref]]

  # Unshrunken results (for p-values and comparison)
  res_unshrunk <- results(dds, name = comp$coef, alpha = 0.05)

  # Shrunken LFC via apeglm
  res_shrunk <- lfcShrink(dds, coef = comp$coef, type = "apeglm")

  # Apply baseMean filter
  res_shrunk_filt <- res_shrunk[res_shrunk$baseMean >= basemean_cutoff, ]
  res_unshrunk_filt <- res_unshrunk[res_unshrunk$baseMean >= basemean_cutoff, ]

  # Count DEGs (using shrunken LFC + original padj)
  n_up <- sum(res_shrunk_filt$padj < 0.05 &
              res_shrunk_filt$log2FoldChange > 1, na.rm = TRUE)
  n_down <- sum(res_shrunk_filt$padj < 0.05 &
                res_shrunk_filt$log2FoldChange < -1, na.rm = TRUE)

  # Also count unshrunken for comparison
  n_up_raw <- sum(res_unshrunk_filt$padj < 0.05 &
                  res_unshrunk_filt$log2FoldChange > 1, na.rm = TRUE)
  n_down_raw <- sum(res_unshrunk_filt$padj < 0.05 &
                    res_unshrunk_filt$log2FoldChange < -1, na.rm = TRUE)

  cat(sprintf("  baseMean >= %d genes: %d\n", basemean_cutoff, nrow(res_shrunk_filt)))
  cat(sprintf("  Unshrunken DEGs: %d up, %d down (%d total)\n",
              n_up_raw, n_down_raw, n_up_raw + n_down_raw))
  cat(sprintf("  Shrunken DEGs:   %d up, %d down (%d total)\n",
              n_up, n_down, n_up + n_down))

  summary_rows[[comp$name]] <- data.frame(
    comparison = comp$name,
    genes_tested = sum(!is.na(res_shrunk_filt$padj)),
    unshrunken_up = n_up_raw,
    unshrunken_down = n_down_raw,
    shrunken_up = n_up,
    shrunken_down = n_down,
    stringsAsFactors = FALSE
  )

  # Save filtered shrunken results with gene names
  res_df <- as.data.frame(res_shrunk_filt)
  res_df$gene_id <- rownames(res_df)
  res_df$gene_name <- gene_names[rownames(res_df)]
  res_df <- res_df[order(res_df$padj), ]
  write.csv(res_df,
    file.path(out_dir, paste0("DEG_shrunk_", comp$name, ".csv")),
    row.names = FALSE)

  all_shrunk[[comp$name]] <- res_shrunk_filt
  all_unshrunk[[comp$name]] <- res_unshrunk_filt

  # ---- MA plot: unshrunken vs shrunken side by side ----
  pdf(file.path(out_dir, paste0("MA_", comp$name, ".pdf")), width = 14, height = 6)
  par(mfrow = c(1, 2))

  # Unshrunken MA
  plotMA(res_unshrunk_filt, main = paste0(comp$name, "\n(unshrunken)"),
         ylim = c(-10, 10), alpha = 0.05)
  abline(h = c(-1, 1), lty = 2, col = "grey40")

  # Shrunken MA
  plotMA(res_shrunk_filt, main = paste0(comp$name, "\n(apeglm shrunken)"),
         ylim = c(-10, 10), alpha = 0.05)
  abline(h = c(-1, 1), lty = 2, col = "grey40")

  dev.off()

  # ---- Volcano plot with shrunken LFC ----
  res_plot <- as.data.frame(res_shrunk_filt)
  res_plot$gene_name <- gene_names[rownames(res_plot)]
  res_plot$significant <- "NS"
  res_plot$significant[res_plot$padj < 0.05 &
                       res_plot$log2FoldChange > 1] <- "Up"
  res_plot$significant[res_plot$padj < 0.05 &
                       res_plot$log2FoldChange < -1] <- "Down"
  res_plot$significant <- factor(res_plot$significant,
                                 levels = c("Down", "NS", "Up"))

  # Label top 15 significant DEGs
  res_plot$label <- NA
  sig_rows <- which(res_plot$significant != "NS")
  if (length(sig_rows) > 0) {
    top_idx <- sig_rows[order(res_plot$padj[sig_rows])][
      1:min(15, length(sig_rows))]
    res_plot$label[top_idx] <- res_plot$gene_name[top_idx]
  }

  p_vol <- ggplot(res_plot,
    aes(log2FoldChange, -log10(pvalue), color = significant)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("Down" = "blue", "NS" = "grey70",
                                  "Up" = "red")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 20,
                    show.legend = FALSE, color = "black") +
    theme_bw(base_size = 12) +
    ggtitle(paste0(comp$name, " (apeglm shrunken, baseMean >= ",
                   basemean_cutoff, ")")) +
    xlab("Shrunken log2 Fold Change") +
    ylab("-log10(p-value)")
  ggsave(file.path(out_dir, paste0("volcano_shrunk_", comp$name, ".pdf")),
         p_vol, width = 8, height = 6)
}

# ---- LFC shrinkage effect: scatter plot for each comparison ----
for (comp in comparisons) {
  unshrunk <- all_unshrunk[[comp$name]]
  shrunk <- all_shrunk[[comp$name]]

  # Match gene IDs
  shared <- intersect(rownames(unshrunk), rownames(shrunk))
  plot_df <- data.frame(
    unshrunk = unshrunk[shared, "log2FoldChange"],
    shrunk = shrunk[shared, "log2FoldChange"],
    baseMean = unshrunk[shared, "baseMean"],
    stringsAsFactors = FALSE
  )

  p_scatter <- ggplot(plot_df, aes(x = unshrunk, y = shrunk, color = log10(baseMean))) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    scale_color_viridis_c(name = "log10(baseMean)") +
    theme_bw(base_size = 12) +
    ggtitle(paste0(comp$name, ": LFC shrinkage effect")) +
    xlab("Unshrunken log2FC") +
    ylab("Shrunken log2FC (apeglm)") +
    coord_fixed()
  ggsave(file.path(out_dir, paste0("shrinkage_scatter_", comp$name, ".pdf")),
         p_scatter, width = 7, height = 6)
}

# ---- Summary table ----
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(out_dir, "shrinkage_summary.csv"),
          row.names = FALSE)

cat("\n=== Shrinkage Summary ===\n")
cat(sprintf("baseMean cutoff: %d\n\n", basemean_cutoff))
print(summary_df, row.names = FALSE)

# ---- Comparison heatmap: max |LFC| before and after shrinkage ----
cat("\n\n=== Max |log2FC| per comparison (before vs after shrinkage) ===\n")
for (comp_name in names(all_shrunk)) {
  max_raw <- max(abs(all_unshrunk[[comp_name]]$log2FoldChange), na.rm = TRUE)
  max_shrunk <- max(abs(all_shrunk[[comp_name]]$log2FoldChange), na.rm = TRUE)
  cat(sprintf("  %-35s  raw: %6.1f  shrunk: %6.1f\n",
              comp_name, max_raw, max_shrunk))
}

cat(sprintf("\n=== Analysis complete ===\nResults in: %s\n", out_dir))
