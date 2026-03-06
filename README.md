# Treg Subset mRNA-seq Analysis

Bulk RNA-seq analysis of four FACS-sorted CD4+ T cell populations from mouse lymphoid tissue, investigating the transcriptional reprogramming that occurs when regulatory T cells (Tregs) lose Foxp3 expression.

## Experimental Design

| Population | Description | n |
|---|---|---|
| Never_Treg | CD4+ T cells that never expressed Foxp3 (conventional T cells) | 3 |
| exTreg | CD4+ T cells that previously expressed Foxp3 but lost expression | 3 |
| Stable_Treg | CD4+ T cells with sustained Foxp3 expression (committed Tregs) | 3 |
| Recent_Treg | CD4+ T cells that recently acquired Foxp3 expression | 3 |

Organism: *Mus musculus* (GRCm39), paired-end unstranded libraries, 60-104M alignment pairs per sample.

## Analysis Pipeline

| Step | Script | Tool |
|---|---|---|
| Alignment | `scripts/01_align_star.sh` | STAR |
| Quantification | `scripts/02_featurecounts.sh` | featureCounts |
| Differential expression | `scripts/03_deseq2_analysis.R` | DESeq2 (paired design: `~ mouse + cell_type`) |
| Pathway analysis | `scripts/04_pathway_analysis.R` | clusterProfiler, fgsea, msigdbr |
| Cell type signatures | `scripts/05_venn_cell_type_pathways.R` | ggVennDiagram, clusterProfiler |
| LFC shrinkage | `scripts/06_lfc_shrinkage.R` | DESeq2 + apeglm |
| exTreg vs Stable_Treg deep dive | `scripts/07_extreg_vs_stable_deep_dive.R` | DESeq2 + apeglm |

## Project Structure

```
aligned/              # STAR alignment output (BAMs)
counts/               # featureCounts gene count matrix
deseq2/               # DESeq2 results (unshrunken), volcano plots, heatmaps, PCA
deseq2_shrinkage/     # apeglm-shrunken results, MA plots, shrinkage scatter plots
pathway_analysis/     # GO, KEGG, GSEA Hallmark & immunologic results per comparison
fastqc/               # QC reports
docs/                 # Detailed analysis results and interpretation
scripts/              # All analysis scripts (numbered in order)
sample_metadata.csv   # Sample metadata
```

## Results

See [docs/analysis_results.md](docs/analysis_results.md) for full results and interpretation.
