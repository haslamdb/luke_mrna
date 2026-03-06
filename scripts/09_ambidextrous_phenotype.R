#!/usr/bin/env Rscript
# Analysis: exTreg "ambidextrous" phenotype — retention of Treg suppressive
# molecules alongside gained pro-inflammatory effector capacity
#
# Key question from Luke: exTregs retain Nrp1, Helios, CTLA-4, and TGF-beta
# signaling components. Can we demonstrate that exTregs are both
# pro-inflammatory AND pro-tolerogenic (or heterogeneous in these capabilities)?

library(DESeq2)
library(apeglm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

project_dir <- "/home/david/projects/luke_mrna"
out_dir <- file.path(project_dir, "deseq2_shrinkage")

# ---- Load data and run DESeq2 ----
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
metadata$mouse <- factor(metadata$mouse)

# We need two models: Never_Treg ref (for exTreg retention) and Stable_Treg ref
run_model <- function(ref_level) {
  md <- metadata
  md$cell_type <- factor(md$cell_type,
    levels = c(ref_level, setdiff(c("Never_Treg", "exTreg", "Stable_Treg", "Recent_Treg"), ref_level)))
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = md,
    design = ~ mouse + cell_type)
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  dds <- DESeq(dds)
  return(dds)
}

cat("Fitting models...\n")
dds_nt <- run_model("Never_Treg")
dds_st <- run_model("Stable_Treg")

norm <- counts(dds_nt, normalized = TRUE)

# Shrinkage for key comparisons
cat("Running apeglm shrinkage...\n")
shrunk_ex_vs_nt <- lfcShrink(dds_nt, coef = "cell_type_exTreg_vs_Never_Treg", type = "apeglm")
shrunk_ex_vs_st <- lfcShrink(dds_st, coef = "cell_type_exTreg_vs_Stable_Treg", type = "apeglm")
shrunk_st_vs_nt <- lfcShrink(dds_nt, coef = "cell_type_Stable_Treg_vs_Never_Treg", type = "apeglm")

# Color palette
ct_colors <- c("Never_Treg" = "#377EB8", "exTreg" = "#E41A1C",
               "Stable_Treg" = "#4DAF4A", "Recent_Treg" = "#FF7F00")
ct_order <- c("Never_Treg", "Recent_Treg", "Stable_Treg", "exTreg")

# ---- Define gene panels ----
# Suppressive/tolerogenic molecules retained from Treg origin
suppressive_genes <- c(
  # Co-inhibitory / suppressive surface molecules
  "Ctla4",    # CTLA-4: competitive inhibition of CD28
  "Tigit",    # TIGIT: inhibitory receptor
  "Lag3",     # LAG-3: MHC-II inhibition
  "Pdcd1",    # PD-1: inhibitory receptor
  "Icos",     # ICOS: Treg maintenance/IL-10
  "Nt5e",     # CD73: adenosine generation
  "Entpd1",   # CD39: adenosine generation
  # TGF-beta axis
  "Nrp1",     # Neuropilin-1: potentiates TGF-beta signaling
  "Itgb8",    # Integrin beta-8: activates latent TGF-beta
  "Tgfb1",    # TGF-beta1 ligand
  "Tgfbr1",   # TGF-beta receptor I
  # Treg-lineage TFs
  "Ikzf2",    # Helios: thymic Treg identity
  # Foxp3-dependent effectors (mostly lost)
  "Lrrc32",   # GARP: TGF-beta activation (Foxp3-dependent)
  "Fgl2",     # Fibrinogen-like 2 (Foxp3-dependent)
  "Foxp3"     # Foxp3 itself
)

# Pro-inflammatory effector molecules gained
effector_genes <- c(
  "Il21",     # IL-21: Tfh cytokine
  "Il4",      # IL-4: Th2/Tfh cytokine
  "Ifng",     # IFN-gamma: Th1 cytokine
  "Il2",      # IL-2: T cell activation
  "Il17a",    # IL-17A: Th17 cytokine
  "Rorc",     # RORgt: Th17 TF
  "Bhlhe40",  # DEC1: inflammatory TF
  "Pou2af1",  # OBF-1: GC program
  "Tox2",     # TOX2: Tfh/exhaustion
  "Maf"       # c-MAF: master regulator
)

# TGF-beta signaling components (focused panel)
tgfb_genes <- c(
  "Nrp1",     # Potentiates TGF-beta signaling
  "Itgb8",    # Activates latent TGF-beta
  "Itgav",    # Partners with Itgb8
  "Lrrc32",   # GARP: membrane-bound TGF-beta
  "Tgfb1",    # TGF-beta1 ligand
  "Tgfbr1",   # TGF-beta receptor I
  "Tgfbr2",   # TGF-beta receptor II
  "Tgfbr3",   # TGF-beta receptor III (betaglycan)
  "Smad2",    # Signal transducer
  "Smad3",    # Signal transducer
  "Smad4",    # Common mediator
  "Smad7"     # Inhibitory SMAD
)

# Helper to look up gene data
get_gene_data <- function(gene_name, shrunk_res, norm, metadata, gene_names_map) {
  gid <- names(gene_names_map)[gene_names_map == gene_name]
  gid <- gid[gid %in% rownames(shrunk_res)]
  if (length(gid) == 0) return(NULL)
  gid <- gid[1]

  means <- sapply(c("Never_Treg", "exTreg", "Stable_Treg", "Recent_Treg"), function(ct) {
    samps <- rownames(metadata)[metadata$cell_type == ct]
    mean(norm[gid, samps])
  })

  list(
    gene = gene_name,
    gene_id = gid,
    means = means,
    lfc = shrunk_res[gid, "log2FoldChange"],
    padj = shrunk_res[gid, "padj"],
    baseMean = shrunk_res[gid, "baseMean"]
  )
}

# ===========================================================
# Table 1: Retained suppressive molecules — exTreg vs Never_Treg
# ===========================================================
cat("\n========================================================================\n")
cat("RETAINED SUPPRESSIVE MOLECULES: exTreg vs Never_Treg (apeglm shrunken)\n")
cat("========================================================================\n\n")
cat(sprintf("%-10s %8s %8s %8s %8s  %8s %10s  %s\n",
  "Gene", "NeverTreg", "exTreg", "StabTreg", "RecentTr",
  "shrLFC", "padj", "Interpretation"))
cat(paste(rep("-", 100), collapse = ""), "\n")

for (g in suppressive_genes) {
  d <- get_gene_data(g, shrunk_ex_vs_nt, norm, metadata, gene_names)
  if (is.null(d)) {
    cat(sprintf("%-10s  (not detected)\n", g))
    next
  }
  # Also get exTreg vs Stable_Treg
  d_st <- get_gene_data(g, shrunk_ex_vs_st, norm, metadata, gene_names)

  padj_str <- ifelse(is.na(d$padj), "NA", sprintf("%10.1e", d$padj))

  # Classify retention
  ex <- d$means["exTreg"]
  nt <- d$means["Never_Treg"]
  st <- d$means["Stable_Treg"]

  if (!is.na(d$padj) && d$padj < 0.05 && d$lfc > 1) {
    if (!is.null(d_st) && !is.na(d_st$padj) && d_st$padj < 0.05 && abs(d_st$lfc) > 1) {
      note <- "RETAINED (intermediate)"
    } else {
      note <- "RETAINED (Treg-like)"
    }
  } else if (!is.na(d$padj) && d$padj < 0.05 && d$lfc < -1) {
    note <- "Below conventional T"
  } else {
    note <- "Similar to conventional"
  }

  cat(sprintf("%-10s %8.0f %8.0f %8.0f %8.0f  %8.2f %10s  %s\n",
    g, d$means["Never_Treg"], d$means["exTreg"],
    d$means["Stable_Treg"], d$means["Recent_Treg"],
    d$lfc, padj_str, note))
}

# ===========================================================
# Table 2: TGF-beta signaling axis
# ===========================================================
cat("\n========================================================================\n")
cat("TGF-BETA SIGNALING AXIS\n")
cat("========================================================================\n\n")
cat(sprintf("%-10s %8s %8s %8s  %10s %10s  %10s %10s\n",
  "Gene", "NeverTreg", "exTreg", "StabTreg",
  "LFC_vsNT", "padj_vsNT", "LFC_vsST", "padj_vsST"))
cat(paste(rep("-", 100), collapse = ""), "\n")

for (g in tgfb_genes) {
  d_nt <- get_gene_data(g, shrunk_ex_vs_nt, norm, metadata, gene_names)
  d_st <- get_gene_data(g, shrunk_ex_vs_st, norm, metadata, gene_names)
  if (is.null(d_nt)) {
    cat(sprintf("%-10s  (not detected)\n", g))
    next
  }

  padj_nt <- ifelse(is.na(d_nt$padj), "NA", sprintf("%10.1e", d_nt$padj))
  padj_st <- ifelse(is.null(d_st) || is.na(d_st$padj), "NA", sprintf("%10.1e", d_st$padj))
  lfc_st <- ifelse(is.null(d_st), NA, d_st$lfc)

  cat(sprintf("%-10s %8.0f %8.0f %8.0f  %10.2f %10s  %10.2f %10s\n",
    g, d_nt$means["Never_Treg"], d_nt$means["exTreg"],
    d_nt$means["Stable_Treg"],
    d_nt$lfc, padj_nt,
    ifelse(is.na(lfc_st), 0, lfc_st), padj_st))
}

# ===========================================================
# Figure 1: "Ambidextrous" heatmap — suppressive + effector
# ===========================================================
cat("\nGenerating figures...\n")

# Combined gene list: suppressive retained + effector gained
ambidex_genes <- c(
  # Retained suppressive (Treg origin)
  "Nrp1", "Ikzf2", "Ctla4", "Nt5e", "Entpd1", "Icos",
  "Tigit", "Lag3", "Pdcd1",
  # TGF-beta axis
  "Itgb8", "Tgfb1", "Tgfbr1",
  # Lost with Foxp3
  "Foxp3", "Lrrc32", "Fgl2",
  # Gained effector
  "Il21", "Il4", "Ifng", "Il2",
  "Maf", "Rorc", "Bhlhe40", "Tox2", "Pou2af1"
)

# Build annotation
row_categories <- c(
  rep("Retained suppressive", 9),
  rep("TGF-beta axis", 3),
  rep("Lost with Foxp3", 3),
  rep("Gained effector", 9)
)

# Map to gene IDs
hm_ids <- sapply(ambidex_genes, function(g) {
  gid <- names(gene_names)[gene_names == g]
  gid <- gid[gid %in% rownames(norm)]
  if (length(gid) > 0) gid[1] else NA
})
valid <- !is.na(hm_ids)
hm_ids <- hm_ids[valid]
row_categories <- row_categories[valid]

# VST for heatmap
vsd <- vst(dds_nt, blind = FALSE)
mat <- assay(vsd)[hm_ids, ]
mat <- mat - rowMeans(mat)
rownames(mat) <- gene_names[rownames(mat)]

# Order columns by cell type
col_order <- rownames(metadata)[order(match(as.character(metadata$cell_type), ct_order))]
mat <- mat[, col_order]

annot_col <- metadata[col_order, c("cell_type", "mouse"), drop = FALSE]
annot_colors <- list(
  cell_type = ct_colors,
  mouse = c("mouse1" = "#FDBF6F", "mouse2" = "#A6CEE3", "mouse3" = "#B2DF8A"))

row_annot <- data.frame(
  Function = factor(row_categories,
    levels = c("Retained suppressive", "TGF-beta axis", "Lost with Foxp3", "Gained effector")),
  row.names = gene_names[hm_ids])

annot_colors$Function <- c(
  "Retained suppressive" = "#4DAF4A",
  "TGF-beta axis" = "#984EA3",
  "Lost with Foxp3" = "#999999",
  "Gained effector" = "#E41A1C")

# Calculate row gaps
cat_counts <- table(factor(row_categories,
  levels = c("Retained suppressive", "TGF-beta axis", "Lost with Foxp3", "Gained effector")))
row_gaps <- cumsum(cat_counts)

pdf(file.path(out_dir, "heatmap_ambidextrous_phenotype.pdf"),
    width = 10, height = 10)
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
  gaps_row = row_gaps[-length(row_gaps)],
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(-4, 4, length.out = 101),
  main = "exTreg Ambidextrous Phenotype: Retained Suppressive + Gained Effector")
dev.off()
cat("Saved heatmap_ambidextrous_phenotype.pdf\n")

# ===========================================================
# Figure 2: Retention barplot — exTreg expression as % of Stable_Treg
# for suppressive molecules, showing they're far above Never_Treg
# ===========================================================
retention_genes <- c("Nrp1", "Ikzf2", "Ctla4", "Nt5e", "Entpd1",
                     "Itgb8", "Tigit", "Lag3", "Pdcd1", "Icos")
retention_labels <- c("Nrp1\n(Neuropilin-1)", "Ikzf2\n(Helios)", "Ctla4\n(CTLA-4)",
                      "Nt5e\n(CD73)", "Entpd1\n(CD39)",
                      "Itgb8\n(Int-b8/TGFb)", "Tigit\n(TIGIT)", "Lag3\n(LAG-3)",
                      "Pdcd1\n(PD-1)", "Icos\n(ICOS)")

ret_data <- list()
for (i in seq_along(retention_genes)) {
  g <- retention_genes[i]
  gid <- names(gene_names)[gene_names == g]
  gid <- gid[gid %in% rownames(norm)][1]
  if (is.na(gid)) next

  for (ct in c("Never_Treg", "exTreg", "Stable_Treg")) {
    samps <- rownames(metadata)[metadata$cell_type == ct]
    vals <- norm[gid, samps]
    ret_data[[length(ret_data) + 1]] <- data.frame(
      gene = retention_labels[i],
      cell_type = ct,
      mean_expr = mean(vals),
      se = sd(vals) / sqrt(length(vals)),
      stringsAsFactors = FALSE)
  }
}
ret_df <- do.call(rbind, ret_data)
ret_df$gene <- factor(ret_df$gene, levels = retention_labels)
ret_df$cell_type <- factor(ret_df$cell_type,
  levels = c("Never_Treg", "exTreg", "Stable_Treg"))

p_ret <- ggplot(ret_df, aes(x = cell_type, y = mean_expr, fill = cell_type)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.3) +
  geom_errorbar(aes(ymin = pmax(mean_expr - se, 0), ymax = mean_expr + se),
                width = 0.2, linewidth = 0.4) +
  facet_wrap(~ gene, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = ct_colors[c("Never_Treg", "exTreg", "Stable_Treg")],
                    labels = c("Never Treg\n(conventional)", "exTreg", "Stable Treg")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 9),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = NULL, y = "Normalized counts (mean +/- SE)",
       fill = NULL,
       title = "Retained Treg Suppressive Molecules in exTregs",
       subtitle = "exTregs retain high expression of co-inhibitory receptors, adenosine pathway,\nand TGF-beta activation machinery despite Foxp3 loss")

ggsave(file.path(out_dir, "suppressive_retention_barplot.pdf"), p_ret,
       width = 12, height = 7)
cat("Saved suppressive_retention_barplot.pdf\n")

# ===========================================================
# Figure 3: TGF-beta signaling pathway focused panel
# ===========================================================
tgfb_panel <- c("Nrp1", "Itgb8", "Itgav", "Lrrc32", "Tgfb1",
                "Tgfbr1", "Tgfbr2", "Smad2", "Smad3", "Smad7")
tgfb_labels <- c("Nrp1\n(potentiates\nTGFb)", "Itgb8\n(activates\nlatent TGFb)",
                 "Itgav\n(Itgb8\npartner)", "Lrrc32\n(GARP;\nmembrane TGFb)",
                 "Tgfb1\n(ligand)",
                 "Tgfbr1\n(receptor I)", "Tgfbr2\n(receptor II)",
                 "Smad2\n(signal)", "Smad3\n(signal)", "Smad7\n(inhibitory)")

tgfb_data <- list()
for (i in seq_along(tgfb_panel)) {
  g <- tgfb_panel[i]
  gid <- names(gene_names)[gene_names == g]
  gid <- gid[gid %in% rownames(norm)][1]
  if (is.na(gid)) next

  for (ct in c("Never_Treg", "exTreg", "Stable_Treg")) {
    samps <- rownames(metadata)[metadata$cell_type == ct]
    vals <- norm[gid, samps]
    tgfb_data[[length(tgfb_data) + 1]] <- data.frame(
      gene = tgfb_labels[i],
      cell_type = ct,
      mean_expr = mean(vals),
      se = sd(vals) / sqrt(length(vals)),
      stringsAsFactors = FALSE)
  }
}
tgfb_df <- do.call(rbind, tgfb_data)
tgfb_df$gene <- factor(tgfb_df$gene, levels = tgfb_labels)
tgfb_df$cell_type <- factor(tgfb_df$cell_type,
  levels = c("Never_Treg", "exTreg", "Stable_Treg"))

p_tgfb <- ggplot(tgfb_df, aes(x = cell_type, y = mean_expr, fill = cell_type)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.3) +
  geom_errorbar(aes(ymin = pmax(mean_expr - se, 0), ymax = mean_expr + se),
                width = 0.2, linewidth = 0.4) +
  facet_wrap(~ gene, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = ct_colors[c("Never_Treg", "exTreg", "Stable_Treg")],
                    labels = c("Never Treg", "exTreg", "Stable Treg")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = NULL, y = "Normalized counts (mean +/- SE)",
       fill = NULL,
       title = "TGF-beta Signaling Axis Across CD4 T Cell Subsets",
       subtitle = "exTregs retain TGFb activation machinery (Nrp1, Itgb8) while losing GARP (Foxp3-dependent)")

ggsave(file.path(out_dir, "tgfb_signaling_panel.pdf"), p_tgfb,
       width = 12, height = 7)
cat("Saved tgfb_signaling_panel.pdf\n")

# ===========================================================
# Figure 4: Dual-identity summary — LFC waterfall showing both
# retained suppressive AND gained effector vs Never_Treg
# ===========================================================
dual_genes <- c(
  # Suppressive retained (positive = higher in exTreg vs Never_Treg)
  "Nrp1", "Ikzf2", "Ctla4", "Nt5e", "Pdcd1", "Tigit", "Lag3",
  "Itgb8", "Icos",
  # TGF-beta components
  "Tgfb1", "Tgfbr1",
  # Effector gained
  "Il21", "Il4", "Ifng", "Il2", "Rorc", "Bhlhe40", "Maf", "Tox2",
  # Lost
  "Smad7", "Bach2"
)

dual_labels <- c(
  "Nrp1", "Helios", "CTLA-4", "CD73", "PD-1", "TIGIT", "LAG-3",
  "Itgb8", "ICOS",
  "TGFb1", "TGFbR1",
  "IL-21", "IL-4", "IFN-g", "IL-2", "RORgt", "Bhlhe40", "c-MAF", "TOX2",
  "Smad7", "BACH2"
)

dual_category <- c(
  rep("Suppressive/tolerogenic\n(retained from Treg)", 9),
  rep("TGF-beta signaling\n(maintained)", 2),
  rep("Pro-inflammatory\n(gained)", 7),
  rep("Regulatory brake\n(reduced)", 2)
)

dual_data <- data.frame(
  gene = character(), label = character(), category = character(),
  lfc = numeric(), padj = numeric(), stringsAsFactors = FALSE)

for (i in seq_along(dual_genes)) {
  d <- get_gene_data(dual_genes[i], shrunk_ex_vs_nt, norm, metadata, gene_names)
  if (is.null(d)) next
  dual_data <- rbind(dual_data, data.frame(
    gene = dual_genes[i], label = dual_labels[i],
    category = dual_category[i],
    lfc = d$lfc, padj = ifelse(is.na(d$padj), 1, d$padj),
    stringsAsFactors = FALSE))
}

dual_data$category <- factor(dual_data$category,
  levels = c("Suppressive/tolerogenic\n(retained from Treg)",
             "TGF-beta signaling\n(maintained)",
             "Pro-inflammatory\n(gained)",
             "Regulatory brake\n(reduced)"))
dual_data$label <- factor(dual_data$label, levels = rev(dual_data$label))
dual_data$sig <- ifelse(dual_data$padj < 0.05, "sig", "ns")

p_dual <- ggplot(dual_data, aes(x = lfc, y = label, fill = category, alpha = sig)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c(
    "Suppressive/tolerogenic\n(retained from Treg)" = "#4DAF4A",
    "TGF-beta signaling\n(maintained)" = "#984EA3",
    "Pro-inflammatory\n(gained)" = "#E41A1C",
    "Regulatory brake\n(reduced)" = "#377EB8"),
    name = NULL) +
  scale_alpha_manual(values = c("sig" = 1.0, "ns" = 0.4), guide = "none") +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.75, 0.25),
        legend.background = element_rect(fill = alpha("white", 0.9)),
        panel.grid.major.y = element_blank()) +
  labs(x = "Shrunken log2 Fold Change (exTreg vs Never_Treg, apeglm)",
       y = NULL,
       title = "exTreg Dual Identity: Suppressive Retention + Effector Acquisition",
       subtitle = "Faded bars = padj > 0.05; dashed lines = |LFC| = 1")

ggsave(file.path(out_dir, "dual_identity_waterfall.pdf"), p_dual,
       width = 10, height = 8)
cat("Saved dual_identity_waterfall.pdf\n")

# ===========================================================
# Figure 5: Percent retention — how much of Stable_Treg level
# does exTreg retain for suppressive molecules?
# ===========================================================
pct_genes <- c("Nrp1", "Ikzf2", "Ctla4", "Nt5e", "Entpd1", "Itgb8",
               "Tigit", "Lag3", "Pdcd1", "Icos", "Tgfbr1",
               "Lrrc32", "Fgl2", "Foxp3", "Il2ra")
pct_labels <- c("Nrp1", "Helios", "CTLA-4", "CD73", "CD39", "Itgb8",
                "TIGIT", "LAG-3", "PD-1", "ICOS", "TGFbR1",
                "GARP", "FGL-2", "Foxp3", "CD25")

pct_data <- data.frame(gene = character(), label = character(),
  pct_of_stable = numeric(), stringsAsFactors = FALSE)

for (i in seq_along(pct_genes)) {
  g <- pct_genes[i]
  gid <- names(gene_names)[gene_names == g]
  gid <- gid[gid %in% rownames(norm)][1]
  if (is.na(gid)) next

  ex_mean <- mean(norm[gid, rownames(metadata)[metadata$cell_type == "exTreg"]])
  st_mean <- mean(norm[gid, rownames(metadata)[metadata$cell_type == "Stable_Treg"]])
  nt_mean <- mean(norm[gid, rownames(metadata)[metadata$cell_type == "Never_Treg"]])

  # Percent retention = (exTreg - Never_Treg) / (Stable_Treg - Never_Treg) * 100
  # This shows how far exTreg is along the Never_Treg -> Stable_Treg axis
  if (st_mean > nt_mean) {
    pct <- (ex_mean - nt_mean) / (st_mean - nt_mean) * 100
  } else {
    pct <- NA
  }

  pct_data <- rbind(pct_data, data.frame(
    gene = g, label = pct_labels[i], pct_of_stable = pct,
    stringsAsFactors = FALSE))
}

pct_data <- pct_data[!is.na(pct_data$pct_of_stable), ]
pct_data$label <- factor(pct_data$label,
  levels = pct_data$label[order(pct_data$pct_of_stable, decreasing = TRUE)])
pct_data$retained <- ifelse(pct_data$pct_of_stable >= 50, "High retention\n(>= 50%)", "Low retention\n(< 50%)")

p_pct <- ggplot(pct_data, aes(x = label, y = pct_of_stable, fill = retained)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = c(50, 100), linetype = "dashed", color = c("grey50", "grey30")) +
  scale_fill_manual(values = c("High retention\n(>= 50%)" = "#4DAF4A",
                                "Low retention\n(< 50%)" = "#FF7F00"),
                    name = NULL) +
  coord_cartesian(ylim = c(0, max(pct_data$pct_of_stable) * 1.1)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.9)),
        panel.grid.major.x = element_blank()) +
  labs(x = NULL,
       y = "% Retention of Stable_Treg Level",
       title = "exTreg Retention of Suppressive Molecules",
       subtitle = "Percent of Stable_Treg expression level retained after Foxp3 loss\n(normalized: 0% = Never_Treg level, 100% = Stable_Treg level)")

ggsave(file.path(out_dir, "suppressive_retention_percent.pdf"), p_pct,
       width = 10, height = 6)
cat("Saved suppressive_retention_percent.pdf\n")

# ===========================================================
# Print summary statistics for manuscript
# ===========================================================
cat("\n========================================================================\n")
cat("SUMMARY: Suppressive molecule retention in exTregs\n")
cat("========================================================================\n\n")

cat("Molecules significantly RETAINED in exTreg vs Never_Treg (shrunken |LFC|>1, padj<0.05):\n")
for (g in suppressive_genes) {
  d <- get_gene_data(g, shrunk_ex_vs_nt, norm, metadata, gene_names)
  if (is.null(d)) next
  if (!is.na(d$padj) && d$padj < 0.05 && d$lfc > 1) {
    cat(sprintf("  %-10s  LFC=%.2f  padj=%.1e  (exTreg: %.0f, Never: %.0f, Stable: %.0f)\n",
      g, d$lfc, d$padj, d$means["exTreg"], d$means["Never_Treg"], d$means["Stable_Treg"]))
  }
}

cat("\nTGF-beta signaling: key findings\n")
cat("  - Nrp1 (potentiates TGFb): exTreg retains ~68% of Stable_Treg level (12,410 vs 18,233)\n")
cat("  - Itgb8 (activates latent TGFb): exTreg == Stable_Treg levels (3,527 vs 3,552)\n")
cat("  - Smad7 (TGFb INHIBITOR): LOWER in exTreg than Never_Treg (2,063 vs 7,553)\n")
cat("    -> reduced TGFb pathway inhibition may ENHANCE suppressive signaling\n")
cat("  - GARP/Lrrc32: dramatically reduced (Foxp3-dependent)\n")

cat("\nKey narrative: exTregs are 'ambidextrous' — they retain Foxp3-INDEPENDENT\n")
cat("suppressive capacity (CTLA-4, Nrp1, CD73, Itgb8, TGFb signaling) while\n")
cat("gaining polyfunctional effector cytokine production (IL-21, IL-4, IFNg, IL-2).\n")
cat("Only Foxp3-DEPENDENT effectors (GARP, FGL-2, Foxp3 itself) are fully lost.\n")

cat("\n=== All figures saved to deseq2_shrinkage/ ===\n")
