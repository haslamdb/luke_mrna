#!/usr/bin/env Rscript
# Deep dive: exTreg vs Stable_Treg — evaluating polyfunctional effector phenotype
# Tests whether baseMean >= 50 filter is too aggressive for cytokine genes

library(DESeq2)
library(apeglm)

project_dir <- "/home/david/projects/luke_mrna"
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
# Stable_Treg as reference — cytokines should be near-zero there
metadata$cell_type <- factor(metadata$cell_type,
  levels = c("Stable_Treg", "Never_Treg", "exTreg", "Recent_Treg"))
metadata$mouse <- factor(metadata$mouse)

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata,
  design = ~ mouse + cell_type)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)

cat("resultsNames:\n")
print(resultsNames(dds))

# Unshrunken and shrunken results — NO baseMean filter
res_raw <- results(dds, name = "cell_type_exTreg_vs_Stable_Treg", alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = "cell_type_exTreg_vs_Stable_Treg", type = "apeglm")

# Normalized counts
norm <- counts(dds, normalized = TRUE)

# Genes of interest
genes_of_interest <- c(
  # Effector cytokines (the polyfunctional question)
  "Ifng", "Tnf", "Il21", "Il4", "Il17a", "Il2",
  # Suppressive cytokines (should be higher in Stable)
  "Il10", "Tgfb1", "Fgl2",
  # Tfh program
  "Bcl6", "Maf", "Cxcr5", "Ascl2", "Sh2d1a", "Slamf6", "Tox2", "Pou2af1",
  # Th1
  "Tbx21", "Stat4", "Cxcr3", "Il18r1",
  # Th17
  "Rorc", "Ccr6", "Bhlhe40", "Hif1a", "Il1r1", "Il23r",
  # Exhaustion
  "Pdcd1", "Lag3", "Tigit", "Tox",
  # GC-associated
  "Pax5", "Sostdc1", "Serpina9",
  # Treg identity
  "Foxp3", "Bach2", "Prdm1", "Gzmb"
)

cat("\n\n========================================================================\n")
cat("exTreg vs Stable_Treg: Key genes — NO baseMean filter applied\n")
cat("========================================================================\n\n")
cat(sprintf("%-12s %8s  %8s %10s  %8s %10s  %8s %8s  %s\n",
  "Gene", "baseMean",
  "rawLFC", "raw_padj",
  "shrLFC", "shr_padj",
  "exTreg", "StabTreg", "Status"))
cat(paste(rep("-", 115), collapse = ""), "\n")

for (g in genes_of_interest) {
  gid <- names(gene_names)[gene_names == g]
  gid <- gid[gid %in% rownames(res_raw)]
  if (length(gid) == 0) {
    cat(sprintf("%-12s  (not in filtered gene set)\n", g))
    next
  }
  gid <- gid[1]

  ex_mean <- mean(norm[gid, metadata$cell_type == "exTreg"])
  st_mean <- mean(norm[gid, metadata$cell_type == "Stable_Treg"])

  raw_lfc <- res_raw[gid, "log2FoldChange"]
  raw_padj <- res_raw[gid, "padj"]
  shr_lfc <- res_shrunk[gid, "log2FoldChange"]
  shr_padj <- res_shrunk[gid, "padj"]
  bm <- res_raw[gid, "baseMean"]

  sig_raw <- !is.na(raw_padj) & raw_padj < 0.05 & abs(raw_lfc) > 1
  sig_shr <- !is.na(shr_padj) & shr_padj < 0.05 & abs(shr_lfc) > 1

  note <- ""
  if (sig_raw & sig_shr) note <- "ROBUST"
  else if (sig_raw & !sig_shr) note <- "lost_w_shrink"
  else if (!sig_raw & sig_shr) note <- "gained"
  else note <- "NS"

  if (bm < 50) note <- paste0(note, " [bM<50]")

  raw_padj_str <- ifelse(is.na(raw_padj), "      NA", sprintf("%10.2e", raw_padj))
  shr_padj_str <- ifelse(is.na(shr_padj), "      NA", sprintf("%10.2e", shr_padj))

  cat(sprintf("%-12s %8.0f  %8.2f %10s  %8.2f %10s  %8.0f %8.0f  %s\n",
    g, bm, raw_lfc, raw_padj_str, shr_lfc, shr_padj_str,
    ex_mean, st_mean, note))
}

# Also show per-replicate counts for the cytokines
cat("\n\n========================================================================\n")
cat("Per-replicate normalized counts for effector cytokines\n")
cat("========================================================================\n\n")

cytokines <- c("Ifng", "Tnf", "Il21", "Il4", "Il17a", "Il2", "Il10")
for (g in cytokines) {
  gid <- names(gene_names)[gene_names == g]
  gid <- gid[gid %in% rownames(norm)]
  if (length(gid) == 0) {
    cat(sprintf("%-8s: not detected\n", g))
    next
  }
  gid <- gid[1]
  cat(sprintf("\n%s (baseMean=%.0f):\n", g, mean(norm[gid, ])))
  for (ct in c("exTreg", "Stable_Treg", "Never_Treg", "Recent_Treg")) {
    samps <- rownames(metadata)[metadata$cell_type == ct]
    vals <- norm[gid, samps]
    cat(sprintf("  %-12s: %s  (mean=%.0f)\n", ct,
      paste(sprintf("%6.0f", vals), collapse = ", "), mean(vals)))
  }
}
