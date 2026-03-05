#!/usr/bin/env Rscript
# Venn diagram of characteristic genes per CD4 T cell population
# + GO enrichment for unique/hallmark pathways per cell type
#
# Usage: Rscript scripts/05_venn_cell_type_pathways.R
# Requires: ggVennDiagram, clusterProfiler, org.Mm.eg.db, ggplot2

library(ggVennDiagram)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)

# ---- Setup ----
project_dir <- "/home/david/projects/luke_mrna"
deseq2_dir <- file.path(project_dir, "deseq2")
cell_types <- c("Never_Treg", "exTreg", "Stable_Treg", "Recent_Treg")

# ---- Find characteristic genes per cell type ----
# For each cell type, find genes significantly UP (padj<0.05, LFC>1)
# vs at least 2 of the 3 other populations
sig_genes <- list()

for (ct in cell_types) {
  others <- setdiff(cell_types, ct)
  up_lists <- list()

  for (other in others) {
    # Try both file name orderings
    f1 <- file.path(deseq2_dir, paste0("DEG_", ct, "_vs_", other, ".csv"))
    f2 <- file.path(deseq2_dir, paste0("DEG_", other, "_vs_", ct, ".csv"))

    if (file.exists(f1)) {
      d <- read.csv(f1)
      # UP in ct = positive LFC (ct is numerator)
      up_lists[[other]] <- d$gene_name[!is.na(d$padj) & d$padj < 0.05 & d$log2FoldChange > 1]
    } else if (file.exists(f2)) {
      d <- read.csv(f2)
      # UP in ct = negative LFC (ct is denominator)
      up_lists[[other]] <- d$gene_name[!is.na(d$padj) & d$padj < 0.05 & d$log2FoldChange < -1]
    }
  }

  # Genes UP vs at least 2 of 3 others
  all_up <- unlist(up_lists)
  freq <- table(all_up)
  sig_genes[[ct]] <- names(freq[freq >= 2])
  cat(sprintf("%s: %d characteristic genes (UP vs >=2 others)\n", ct, length(sig_genes[[ct]])))
}

names(sig_genes) <- c("Never Treg", "exTreg", "Stable Treg", "Recent Treg")

# ---- Venn diagram ----
p_venn <- ggVennDiagram(sig_genes,
  label = "count",
  label_alpha = 0,
  label_size = 5,
  set_color = c("#4DAF4A", "#E41A1C", "#377EB8", "#FF7F00"),
  edge_size = 1.5) +
  scale_fill_gradient(low = "white", high = "#d9e8fb") +
  scale_color_manual(values = c("#4DAF4A", "#E41A1C", "#377EB8", "#FF7F00")) +
  ggtitle("Characteristic genes per CD4 T cell population\n(UP vs at least 2 of 3 other populations, padj<0.05, LFC>1)") +
  theme(plot.title = element_text(hjust = 0.5, size = 13))

ggsave(file.path(deseq2_dir, "venn_characteristic_genes.pdf"), p_venn, width = 10, height = 8)
cat("\nVenn diagram saved\n")

# ---- Gene name to Entrez mapping ----
all_names <- unique(unlist(sig_genes))
name_map <- bitr(all_names, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# ---- Unique genes and GO enrichment per cell type ----
cat("\n========== UNIQUE GENES & PATHWAYS PER CELL TYPE ==========\n")

for (i in seq_along(sig_genes)) {
  ct <- names(sig_genes)[i]
  all_genes <- sig_genes[[i]]
  others_genes <- unique(unlist(sig_genes[-i]))
  unique_genes <- setdiff(all_genes, others_genes)
  shared_genes <- intersect(all_genes, others_genes)

  cat(sprintf("\n--- %s: %d total, %d unique, %d shared ---\n",
              ct, length(all_genes), length(unique_genes), length(shared_genes)))

  # Top unique genes
  if (length(unique_genes) > 0) {
    cat("  Unique genes: ", paste(head(sort(unique_genes), 25), collapse = ", "), "\n")
  }

  # GO enrichment on unique genes
  entrez_unique <- name_map$ENTREZID[name_map$SYMBOL %in% unique_genes]
  if (length(entrez_unique) >= 5) {
    ego <- enrichGO(gene = entrez_unique, OrgDb = org.Mm.eg.db, ont = "BP",
                    pAdjustMethod = "BH", pvalueCutoff = 0.1, readable = TRUE)
    if (nrow(ego) > 0) {
      ego_df <- as.data.frame(ego)
      top <- head(ego_df[order(ego_df$p.adjust), ], 8)
      cat("  Hallmark pathways (unique genes):\n")
      for (j in 1:nrow(top)) {
        cat(sprintf("    - %-50s (p.adj=%.1e, n=%s)\n",
                    top$Description[j], top$p.adjust[j], top$Count[j]))
      }
    }
  }

  # GO enrichment on ALL characteristic genes
  entrez_all <- name_map$ENTREZID[name_map$SYMBOL %in% all_genes]
  if (length(entrez_all) >= 10) {
    ego2 <- enrichGO(gene = entrez_all, OrgDb = org.Mm.eg.db, ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
    if (nrow(ego2) > 0) {
      ego2_df <- as.data.frame(ego2)
      top2 <- head(ego2_df[order(ego2_df$p.adjust), ], 8)
      cat("  Top pathways (all characteristic genes):\n")
      for (j in 1:nrow(top2)) {
        cat(sprintf("    - %-50s (p.adj=%.1e, n=%s)\n",
                    top2$Description[j], top2$p.adjust[j], top2$Count[j]))
      }
      ct_clean <- gsub(" ", "_", ct)
      write.csv(ego2_df, file.path(deseq2_dir, paste0("GO_characteristic_", ct_clean, ".csv")),
                row.names = FALSE)
    }
  }
}

cat(sprintf("\n=== Complete. Results in: %s ===\n", deseq2_dir))
