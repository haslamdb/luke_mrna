#!/usr/bin/env Rscript
# Pathway analysis for Luke's CD4 T cell RNA-seq
# ORA (GO/KEGG) on significant DEGs + GSEA on ranked gene lists
# Uses all pairwise comparisons from DESeq2

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(fgsea)
library(msigdbr)
library(ggplot2)

# ---- Setup ----
project_dir <- "/home/david/projects/luke_mrna"
deseq2_dir <- file.path(project_dir, "deseq2")
out_dir <- file.path(project_dir, "pathway_analysis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Get all DEG result files
deg_files <- list.files(deseq2_dir, pattern = "^DEG_.*\\.csv$", full.names = TRUE)
cat(sprintf("Found %d DEG files\n", length(deg_files)))

# ---- Gene ID mapping ----
# Map Ensembl IDs to Entrez IDs (required for KEGG, useful for GO)
# Build mapping once from org.Mm.eg.db
all_ensembl <- unique(unlist(lapply(deg_files, function(f) {
  read.csv(f)$gene_id
})))

id_map <- bitr(all_ensembl, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"),
               OrgDb = org.Mm.eg.db)
cat(sprintf("Mapped %d / %d Ensembl IDs to Entrez\n", nrow(id_map), length(all_ensembl)))

# ---- MSigDB gene sets for GSEA ----
# Hallmark gene sets (mouse)
hallmark_sets <- msigdbr(species = "Mus musculus", category = "H")
hallmark_list <- split(hallmark_sets$ensembl_gene, hallmark_sets$gs_name)

# C7 immunologic signatures
immuno_sets <- msigdbr(species = "Mus musculus", category = "C7", subcategory = "IMMUNESIGDB")
immuno_list <- split(immuno_sets$ensembl_gene, immuno_sets$gs_name)

# ---- Process each comparison ----
for (deg_file in deg_files) {
  contrast_name <- gsub("DEG_(.*)\\.csv", "\\1", basename(deg_file))
  cat(sprintf("\n========== %s ==========\n", contrast_name))
  contrast_dir <- file.path(out_dir, contrast_name)
  dir.create(contrast_dir, showWarnings = FALSE)

  res <- read.csv(deg_file)

  # Significant DEGs (padj < 0.05, |LFC| > 1)
  sig <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
  sig_up <- sig[sig$log2FoldChange > 0, ]
  sig_down <- sig[sig$log2FoldChange < 0, ]

  # Map to Entrez
  sig_entrez <- id_map$ENTREZID[id_map$ENSEMBL %in% sig$gene_id]
  bg_entrez <- id_map$ENTREZID[id_map$ENSEMBL %in% res$gene_id]

  # ---- GO over-representation (BP) ----
  if (length(sig_entrez) >= 10) {
    ego <- enrichGO(gene = sig_entrez,
                    universe = bg_entrez,
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.1,
                    readable = TRUE)

    if (nrow(ego) > 0) {
      write.csv(as.data.frame(ego), file.path(contrast_dir, "GO_BP.csv"), row.names = FALSE)

      p <- dotplot(ego, showCategory = 20) + ggtitle(paste("GO BP:", contrast_name))
      ggsave(file.path(contrast_dir, "GO_BP_dotplot.pdf"), p, width = 10, height = 8)

      cat(sprintf("  GO BP: %d significant terms\n", nrow(ego)))
    } else {
      cat("  GO BP: no significant terms\n")
    }
  } else {
    cat(sprintf("  GO BP: skipped (only %d DEGs with Entrez IDs)\n", length(sig_entrez)))
  }

  # ---- KEGG over-representation ----
  if (length(sig_entrez) >= 10) {
    ekegg <- enrichKEGG(gene = sig_entrez,
                        universe = bg_entrez,
                        organism = "mmu",
                        pvalueCutoff = 0.05)

    if (!is.null(ekegg) && nrow(ekegg) > 0) {
      write.csv(as.data.frame(ekegg), file.path(contrast_dir, "KEGG.csv"), row.names = FALSE)

      p <- dotplot(ekegg, showCategory = 20) + ggtitle(paste("KEGG:", contrast_name))
      ggsave(file.path(contrast_dir, "KEGG_dotplot.pdf"), p, width = 10, height = 8)

      cat(sprintf("  KEGG: %d significant pathways\n", nrow(ekegg)))
    } else {
      cat("  KEGG: no significant pathways\n")
    }
  }

  # ---- GSEA with fgsea (Hallmark) ----
  # Rank all genes by signed -log10(pvalue)
  ranked <- res[!is.na(res$pvalue) & res$pvalue > 0, ]
  ranks <- sign(ranked$log2FoldChange) * -log10(ranked$pvalue)
  names(ranks) <- ranked$gene_id
  ranks <- sort(ranks, decreasing = TRUE)

  fgsea_hallmark <- fgsea(pathways = hallmark_list,
                          stats = ranks,
                          minSize = 15,
                          maxSize = 500)
  fgsea_hallmark <- fgsea_hallmark[order(fgsea_hallmark$pval), ]

  if (nrow(fgsea_hallmark[fgsea_hallmark$padj < 0.05, ]) > 0) {
    # Save results (drop leadingEdge list column for CSV)
    fh_save <- fgsea_hallmark
    fh_save$leadingEdge <- sapply(fh_save$leadingEdge, paste, collapse = ";")
    write.csv(fh_save, file.path(contrast_dir, "GSEA_hallmark.csv"), row.names = FALSE)

    # Plot top pathways
    top_paths <- head(fgsea_hallmark[fgsea_hallmark$padj < 0.05, ], 20)
    top_paths$pathway <- gsub("HALLMARK_", "", top_paths$pathway)

    p <- ggplot(top_paths, aes(x = NES, y = reorder(pathway, NES), fill = padj)) +
      geom_col() +
      scale_fill_gradient(low = "red", high = "blue") +
      theme_bw(base_size = 12) +
      xlab("Normalized Enrichment Score") +
      ylab(NULL) +
      ggtitle(paste("GSEA Hallmark:", contrast_name))
    ggsave(file.path(contrast_dir, "GSEA_hallmark_barplot.pdf"), p, width = 10, height = 7)

    cat(sprintf("  GSEA Hallmark: %d significant (padj<0.05)\n",
                nrow(fgsea_hallmark[fgsea_hallmark$padj < 0.05, ])))
  } else {
    cat("  GSEA Hallmark: no significant pathways\n")
  }

  # ---- GSEA with fgsea (Immunologic C7) ----
  fgsea_immuno <- fgsea(pathways = immuno_list,
                        stats = ranks,
                        minSize = 15,
                        maxSize = 500)
  fgsea_immuno <- fgsea_immuno[order(fgsea_immuno$pval), ]

  sig_immuno <- fgsea_immuno[fgsea_immuno$padj < 0.05, ]
  if (nrow(sig_immuno) > 0) {
    fi_save <- fgsea_immuno
    fi_save$leadingEdge <- sapply(fi_save$leadingEdge, paste, collapse = ";")
    write.csv(fi_save, file.path(contrast_dir, "GSEA_immunologic.csv"), row.names = FALSE)
    cat(sprintf("  GSEA Immunologic (C7): %d significant\n", nrow(sig_immuno)))
  } else {
    cat("  GSEA Immunologic (C7): no significant pathways\n")
  }
}

# ---- Summary across comparisons ----
cat("\n\n========== Cross-comparison Hallmark summary ==========\n")
all_hallmark <- do.call(rbind, lapply(list.files(out_dir, pattern = "GSEA_hallmark.csv",
                                                  recursive = TRUE, full.names = TRUE),
  function(f) {
    d <- read.csv(f)
    d$comparison <- basename(dirname(f))
    d
  }))

if (nrow(all_hallmark) > 0) {
  # Pivot: pathway x comparison NES matrix
  sig_hall <- all_hallmark[all_hallmark$padj < 0.05, ]
  pathway_counts <- sort(table(sig_hall$pathway), decreasing = TRUE)
  cat("Hallmark pathways significant in multiple comparisons:\n")
  multi <- pathway_counts[pathway_counts >= 2]
  if (length(multi) > 0) {
    for (i in seq_along(multi)) {
      cat(sprintf("  %s: %d comparisons\n", names(multi)[i], multi[i]))
    }
  }
  write.csv(all_hallmark, file.path(out_dir, "all_hallmark_results.csv"), row.names = FALSE)
}

cat(sprintf("\n=== Pathway analysis complete ===\nResults in: %s\n", out_dir))
