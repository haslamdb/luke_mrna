# CD4 T Cell Subset RNA-seq Analysis: exTregs Acquire a Polyfunctional Effector Program

## Overview

Bulk RNA-seq analysis of four FACS-sorted CD4 T cell populations from mouse lymphoid tissue reveals that ex-regulatory T cells (exTregs) undergo dramatic transcriptional reprogramming upon Foxp3 loss. exTregs retain markers of their thymic Treg origin (Nrp1, Helios) and follicular identity (CXCR5, TIGIT) while gaining polyfunctional effector cytokine production (IL-21, IL-4, IFN-gamma, IL-2) and activating a germinal center-associated transcriptional network. These findings are consistent with a Treg to Tfr to pathogenic ex-Foxp3 conversion pathway, where Foxp3 loss derepresses effector loci under the control of extraordinarily elevated c-MAF.

All differential expression results reported here use **apeglm log2 fold change shrinkage** (Zhu et al., *Bioinformatics* 2018), which applies empirical Bayes shrinkage to penalize unreliable fold changes from low-count or high-variance genes while preserving well-supported differences.

## Experimental Design

**Populations sorted (4 cell types, 3 biological replicates each):**

| Population | Abbreviation | Description | n |
|---|---|---|---|
| Never Treg | Never_Treg | CD4+ T cells that never expressed Foxp3 (conventional T cells) | 3 |
| ex-Treg | exTreg | CD4+ T cells that previously expressed Foxp3 but lost expression | 3 |
| Stable Treg | Stable_Treg | CD4+ T cells with sustained Foxp3 expression (committed Tregs) | 3 |
| Recent Treg | Recent_Treg | CD4+ T cells that recently acquired Foxp3 expression | 3 |

- **Organism**: Mus musculus (GRCm39)
- **Biological replicates**: 3 mice per group (mouse1/Nov11, mouse2/Nov19, mouse3/Nov21)
- **Library type**: Paired-end, unstranded (confirmed via STAR ReadsPerGene.out.tab)
- **Sequencing depth**: 60-104M alignment pairs per sample

## Analysis Pipeline

| Step | Script | Tool | Key Parameters |
|---|---|---|---|
| 1. Alignment | `01_align_star.sh` | STAR 2.7.x | Aligned to GRCm39 mouse genome |
| 2. Quantification | `02_featurecounts.sh` | featureCounts 2.0.6 | Paired-end (`-p --countReadPairs`), unstranded (`-s 0`) |
| 3. Differential expression | `03_deseq2_analysis.R` | DESeq2 1.48.2 | Paired design: `~ mouse + cell_type` |
| 4. Pathway analysis | `04_pathway_analysis.R` | clusterProfiler, fgsea, msigdbr | GO BP, KEGG (ORA); Hallmark & C7 immunologic (GSEA) |
| 5. Cell type signatures | `05_venn_cell_type_pathways.R` | ggVennDiagram, clusterProfiler | Characteristic genes UP vs at least 2 of 3 other populations |
| 6. LFC shrinkage | `06_lfc_shrinkage.R` | DESeq2 + apeglm | All 6 pairwise comparisons, baseMean >= 50 |
| 7. exTreg vs Stable_Treg | `07_extreg_vs_stable_deep_dive.R` | DESeq2 + apeglm | Targeted comparison, no baseMean filter |

### Alignment & Quantification Summary

| Sample | Total Alignments | Assigned | Rate |
|---|---|---|---|
| exTreg (mouse1) | 86M | 43.4M | 50.5% |
| exTreg-Nov19 (mouse2) | 90M | 43.1M | 47.7% |
| exTreg-Nov21 (mouse3) | 79M | 38.1M | 48.4% |
| Never-Treg-Nov11 | 91M | 44.9M | 49.2% |
| Never-Treg-Nov19-2 | 71M | 33.3M | 46.8% |
| Never-Treg-Nov21 | 103M | 50.5M | 49.0% |
| Recent-Treg | 66M | 32.1M | 48.6% |
| Recent-Treg-Nov19 | 89M | 42.3M | 47.4% |
| Recent-Treg-Nov21 | 60M | 30.4M | 50.9% |
| Stable-Treg-Nov11 | 104M | 49.8M | 48.1% |
| Stable-Treg-Nov19-2 | 90M | 41.2M | 46.0% |
| Stable-Treg-Nov21-2 | 93M | 43.1M | 46.4% |

**16,075 genes** passed filtering (at least 10 counts in at least 3 samples).

---

## QC and Global Structure

- **PCA**: Clear separation of all four populations. PC1 separates Treg-identity (Foxp3+) populations from Foxp3-low/negative populations; PC2 distinguishes exTreg from Never_Treg.
- **Sample distance heatmap**: Replicates cluster tightly by cell type.
- **Output files**: `deseq2/pca_plot.pdf`, `deseq2/sample_distance_heatmap.pdf`

### Differential Expression Summary (apeglm shrunken, baseMean >= 50)

| Comparison | Shrunken DEGs (up) | Shrunken DEGs (down) | Total |
|---|---|---|---|
| exTreg vs Never_Treg | 264 | 63 | 327 |
| Stable_Treg vs Never_Treg | 630 | 736 | 1,366 |
| Recent_Treg vs Never_Treg | 121 | 25 | 146 |
| Stable_Treg vs exTreg | 306 | 560 | 866 |
| Recent_Treg vs exTreg | 61 | 35 | 96 |
| Recent_Treg vs Stable_Treg | 330 | 159 | 489 |

---

## Treg Identity Genes Validate Sort Purity

Foxp3 and core Treg markers confirm population definitions. All values are mean DESeq2-normalized counts.

| Gene | Protein | Never_Treg | exTreg | Stable_Treg | Recent_Treg | Biology |
|---|---|---|---|---|---|---|
| **Foxp3** | FOXP3 | 101 | 191 | **69,700** | 27,769 | Master Treg TF |
| **Il2ra** | CD25 | 535 | 5,126 | **86,237** | 41,930 | IL-2 receptor, Treg survival |
| **Ikzf2** | Helios | 2,503 | 39,172 | **79,792** | 14,169 | Thymic Treg marker |
| **Ctla4** | CTLA-4 | 1,117 | 7,918 | **27,812** | 18,772 | Co-inhibitory receptor |
| **Nrp1** | Neuropilin-1 | 413 | 12,410 | **18,233** | 2,765 | Thymic Treg origin marker |
| **Lrrc32** | GARP | 58 | 980 | **47,187** | 12,892 | TGF-beta activation |
| **Nt5e** | CD73 | 1,959 | 16,153 | **29,102** | 10,238 | Adenosine pathway |
| **Entpd1** | CD39 | 248 | 738 | **3,471** | 1,805 | Adenosine pathway |

**Key observations:**
- Stable_Treg has the strongest Treg signature across all markers
- Recent_Treg shows intermediate levels, consistent with a newly committed state
- exTreg retains **high Nrp1** (12,410) and **high Ikzf2/Helios** (39,172), confirming thymic Treg origin even after Foxp3 loss

---

## Robust Findings Across All Comparisons (exTreg vs Never_Treg, apeglm shrunken)

### Inhibitory Receptors

exTregs express the highest levels of multiple inhibitory receptors. PD-1, LAG-3, and TIGIT are among the most robust findings in the entire dataset:

| Gene | Protein | Shrunken LFC | padj | Never_Treg | exTreg | Stable_Treg | Recent_Treg |
|---|---|---|---|---|---|---|---|
| **Pdcd1** | PD-1 | 5.61 | 8.3e-05 | 92 | **3,309** | 397 | 409 |
| **Lag3** | LAG-3 | 4.45 | 2.1e-05 | 69 | **1,654** | 800 | 815 |
| **Tigit** | TIGIT | 4.58 | 9.4e-15 | 284 | **6,020** | 4,894 | 1,180 |
| **Tox2** | TOX2 | 2.31 | 4.4e-03 | 421 | **2,775** | 288 | 145 |

Note: LAG-3 and TIGIT are significant vs Never_Treg but not vs Stable_Treg, because Stable_Tregs (Tfr cells) also express them. PD-1 is significant in both comparisons.

### Treg Origin Markers Retained in exTregs

| Gene | Protein | Shrunken LFC | padj |
|---|---|---|---|
| Nrp1 | Neuropilin-1 | 4.87 | 9.4e-56 |
| Ikzf2 | Helios | 3.87 | 2.7e-22 |
| Ctla4 | CTLA-4 | 2.80 | 3.1e-21 |
| Il2ra | CD25 | 3.16 | 3.1e-18 |

### Tfh-Associated Markers (vs Never_Treg)

| Gene | Protein | Shrunken LFC | padj |
|---|---|---|---|
| Maf | c-MAF | 2.53 | 6.0e-13 |
| Cxcr5 | CXCR5 | 2.27 | 2.0e-03 |
| Slamf6 | Ly108 | 1.05 | 1.2e-05 |
| Sh2d1a | SAP | 1.22 | 5.8e-03 |

CXCR5 and Slamf6 are robustly upregulated vs Never_Treg but not vs Stable_Treg, because Stable_Tregs (Tfr cells) also express these follicular markers. This is consistent with the Treg to Tfr to exTreg conversion model.

### Th17-Associated Genes (vs Never_Treg)

| Gene | Protein | Shrunken LFC | padj |
|---|---|---|---|
| Rorc | RORgt | 5.15 | 5.3e-05 |
| Ccr6 | CCR6 | 3.04 | 6.6e-04 |
| Bhlhe40 | DEC1 | 3.21 | 5.8e-04 |
| Hif1a | HIF-1a | 1.21 | 3.4e-04 |

RORgt, CCR6, and Bhlhe40 are robustly upregulated vs Never_Treg. HIF-1a (22,605 counts in exTreg) is the highest across all populations and is known to shift the Treg/Th17 balance toward Th17 by degrading Foxp3 and activating RORgt.

### GC-Associated Genes

Several genes normally associated with germinal center B cells are robustly and specifically expressed in exTregs:

| Gene | Protein | Shrunken LFC (vs NT) | padj | Shrunken LFC (vs ST) | padj |
|---|---|---|---|---|---|
| Pax5 | PAX5 | 9.27 | 3.8e-05 | 13.22 | 7.4e-07 |
| Sostdc1 | SOSTDC1 | 7.83 | 5.0e-05 | 16.44 | 8.8e-11 |
| Serpina9 | SERPINA9 | 6.88 | 7.1e-04 | 12.74 | 7.2e-06 |
| Cd22 | CD22 | 5.57 | 7.2e-03 | — | — |

These survive apeglm shrinkage with very large LFCs, confirming they are not low-count artifacts. This is **not B cell contamination**: CD19 and CD20 (Ms4a1) are absent (filtered out for too few counts), Ebf1 and Blnk are absent, Ig genes do not track with exTreg, and T cell identity markers (CD3e 42,226; CD4 120,248) are very high. The expression of these genes without activation of core B cell targets indicates germinal center-associated transcriptional network activation rather than B cell transdifferentiation.

### BACH2 Downregulation

BACH2 is dramatically reduced in exTregs (shrunken LFC = -1.98, padj = 4.2e-08; 2,351 vs 10,314 in Never_Treg). BACH2 is a transcriptional repressor that maintains T cell quiescence and prevents effector differentiation. Its loss likely contributes to the derepression of multiple effector programs.

### Chemokine Receptor Profile

| Receptor | Shrunken LFC (vs NT) | padj | Never_Treg | exTreg | Stable_Treg | Homing |
|---|---|---|---|---|---|---|
| **CXCR5** | 2.27 | 2.0e-03 | 1,128 | **6,875** | 5,310 | B cell follicles / GC |
| **CXCR3** | 2.11 | 1.7e-02 | 917 | **5,400** | 2,332 | Inflamed tissue |
| **CCR6** | 3.04 | 6.6e-04 | 225 | **2,189** | 1,674 | Mucosal / Th17 |
| CCR7 | — | NS | **6,960** | 2,748 | 4,079 | Lymph node T zone |
| S1PR1 | — | NS | **14,400** | 7,658 | 9,746 | Lymph node egress |

The CXCR5+CXCR3+ co-expression pattern is characteristic of "Tfh1" cells described in autoimmune and chronic inflammatory settings.

---

## Bcl6: A Note on Post-Transcriptional Regulation

Bcl6 does not reach significance after apeglm shrinkage in any comparison (LFC = 0.38, padj = 0.49 vs Never_Treg; LFC = 0.77, padj = 0.25 vs Stable_Treg). The raw normalized counts show exTreg (4,208) modestly above other populations (1,503-2,181), but this difference is not robust at n=3.

This does not rule out a role for Bcl6 in exTreg biology. Bcl6 protein stability is heavily regulated post-translationally, and modest RNA differences may understate protein-level differences. Luke's flow cytometry data can address whether Bcl6 protein is differentially expressed. At the transcriptomic level, however, **c-MAF rather than Bcl6 appears to be the primary transcriptional driver** of the Tfh-like features in exTregs.

---

## Pathway Analysis

### GSEA Hallmark Pathways

**IL2-STAT5 signaling** was the dominant pathway, significant in 5 of 6 comparisons:

| Comparison | NES | padj | Interpretation |
|---|---|---|---|
| Stable_Treg vs Never_Treg | +2.67 | 2.2e-16 | Core Treg survival pathway |
| exTreg vs Never_Treg | +2.23 | 8.0e-12 | Retained from Treg origin |
| Recent_Treg vs Never_Treg | +2.23 | 2.0e-07 | Active in new Tregs |
| Stable_Treg vs exTreg | +2.17 | 6.5e-03 | Higher in committed Tregs |
| Recent_Treg vs Stable_Treg | -2.47 | 5.1e-16 | Lower in Recent vs Stable |

Additional pathways in exTreg vs Never_Treg:
- **E2F targets** (NES=1.96, padj=4.1e-06) — cell cycle/proliferation
- **G2M checkpoint** (NES=1.88, padj=7.4e-05) — active division

### GO Biological Process Enrichment

**exTreg vs Never_Treg (22 immune GO terms):**
- T cell activation (p.adj=7.6e-04)
- Alpha-beta T cell activation and differentiation
- Negative regulation of T cell mediated immunity (p.adj=1.5e-03)
- Tolerance induction (p.adj=6.2e-03)

**Stable_Treg vs exTreg (31 immune GO terms):**
- Regulation of cytokine production in immune response (p.adj=1.5e-03)
- Regulation of T cell activation (p.adj=7.8e-03)
- Regulation of Th17 cell differentiation (p.adj=8.9e-03)

### Cell Type Signature Analysis

Characteristic genes per population (UP vs at least 2 of 3 other populations):

| Cell Type | Total | Unique | Top Pathways |
|---|---|---|---|
| **Never Treg** | 103 | 61 | NK cell cytotoxicity, cell killing, lymphocyte differentiation |
| **exTreg** | 150 | 91 | Inflammatory response regulation, immune effector regulation |
| **Stable Treg** | 293 | 244 | T cell activation/differentiation, hemopoiesis regulation |
| **Recent Treg** | 103 | 20 | Leukocyte homeostasis, lymphocyte proliferation regulation |

Stable_Treg has the most unique genes (244), consistent with a fully committed suppressive program. Recent_Treg has the fewest unique genes (20), consistent with a transitional state.

---

## exTreg vs Stable_Treg: The Polyfunctional Effector Phenotype

The exTreg vs Stable_Treg comparison is the most informative for understanding what exTregs gain upon Foxp3 loss. Stable_Tregs have near-zero expression of effector cytokines due to active Foxp3-mediated suppression, making this the cleanest test of acquired effector function.

### Effector Cytokines Are Robustly Upregulated

With apeglm shrinkage applied, four of six effector cytokines are significant:

| Gene | Protein | exTreg | Stable_Treg | Shrunken LFC | padj |
|---|---|---|---|---|---|
| **Il21** | IL-21 | 4,529 | 0 | 15.75 | 2.0e-18 |
| **Il4** | IL-4 | 387 | 0 | 11.67 | 4.7e-06 |
| **Ifng** | IFN-gamma | 1,255 | 28 | 4.73 | 1.2e-03 |
| **Il2** | IL-2 | 630 | 17 | 3.88 | 1.3e-03 |
| Il17a | IL-17A | 113 | 0 | 0.13 | 0.14 |
| Tnf | TNF-alpha | 2,709 | 1,333 | 0.63 | 0.22 |

IL-21 and IL-4 are completely absent in all three Stable_Treg replicates and consistently expressed across all three exTreg replicates. IFN-gamma is produced without T-bet upregulation (Tbx21: 175 vs 23, NS), consistent with c-MAF- or STAT4-dependent IFN-gamma that is mechanistically distinct from classical Th1 differentiation.

**Per-replicate normalized counts:**

| Gene | exTreg rep1 | rep2 | rep3 | Stable_Treg rep1 | rep2 | rep3 |
|---|---|---|---|---|---|---|
| Il21 | 1,937 | 7,774 | 3,875 | 0 | 0 | 0 |
| Il4 | 389 | 470 | 300 | 0 | 0 | 0 |
| Ifng | 1,159 | 1,157 | 1,449 | 85 | 0 | 0 |
| Il2 | 852 | 96 | 943 | 8 | 39 | 3 |
| Il17a | 61 | 5 | 273 | 0 | 0 | 0 |

IL-17A has genuine high replicate variance (61/5/273), explaining why shrinkage appropriately penalizes it. TNF is not significant because Stable_Tregs also express it (3,776/110/113).

### Tfh/GC Genes Gained by exTregs

| Gene | Protein | Shrunken LFC | padj |
|---|---|---|---|
| Pou2af1 | OBF-1 | 4.57 | 5.6e-04 |
| Sh2d1a | SAP | 1.32 | 2.4e-03 |
| Slamf6 | Ly108 | 1.46 | 4.0e-10 |
| Tox2 | TOX2 | 2.97 | 1.5e-04 |
| Bhlhe40 | DEC1 | 2.74 | 2.0e-03 |
| Pdcd1 | PD-1 | 2.37 | 2.0e-02 |

Pou2af1 was not significant vs Never_Treg but is robustly upregulated vs Stable_Treg (2,463 vs 89), further supporting the GC-associated transcriptional program.

### Treg Suppressive Functions Lost

| Gene | Protein | Shrunken LFC | padj | exTreg | Stable_Treg |
|---|---|---|---|---|---|
| Foxp3 | FOXP3 | -8.49 | 8.6e-252 | 191 | 69,700 |
| Fgl2 | FGL-2 | -4.90 | 2.6e-08 | 187 | 4,933 |
| Bach2 | BACH2 | -1.50 | 4.8e-05 | 2,351 | 7,246 |
| Prdm1 | BLIMP-1 | -1.57 | 2.0e-03 | 859 | 2,983 |

### Shared Features Between exTreg and Stable_Treg

Several markers are **not** significantly different between exTreg and Stable_Treg, revealing shared biology:

| Gene | Protein | exTreg | Stable_Treg | Interpretation |
|---|---|---|---|---|
| Cxcr5 | CXCR5 | 6,875 | 5,310 | Both populations contain follicular cells (Tfr/Tfh) |
| Tigit | TIGIT | 6,020 | 4,894 | Shared follicular/inhibitory marker |
| Lag3 | LAG-3 | 1,654 | 800 | Shared inhibitory receptor |
| Maf | c-MAF | 61,705 | 29,378 | Both express high c-MAF (padj=0.013, but |LFC| < 1) |
| Bcl6 | BCL-6 | 4,208 | 1,986 | Modest difference not significant after shrinkage |

This pattern is consistent with the Treg to Tfr to exTreg conversion model: exTregs retain the follicular surface markers of their Tfr precursors while gaining effector cytokine production that was suppressed by Foxp3.

---

## Proposed Model: Treg to Tfr to Pathogenic ex-Foxp3 Conversion

```
Thymic Treg (Foxp3+, Nrp1+, Helios+)
       |
       v
Follicular Regulatory T cell (Tfr)
(Foxp3+, CXCR5+, PD-1+, TIGIT+, c-MAF high)
       |
       |  Foxp3 destabilization:
       |  - HIF-1a upregulation (22,605 counts)
       |  - c-Maf drives IL-21 production
       |  - Loss of IL-2 signaling
       v
ex-Treg (Foxp3-, Nrp1+, Helios+)
+------------------------------------------+
|  Retained from Tfr:                      |
|  CXCR5, TIGIT, PD-1, c-MAF, Slamf6      |
|                                          |
|  Gained effector functions:              |
|  IL-21, IL-4, IFN-gamma, IL-2 (robust)  |
|  Pou2af1, Bhlhe40, Tox2                 |
|                                          |
|  GC-associated network:                  |
|  Pax5, Sostdc1, Serpina9                 |
|                                          |
|  Lost suppressive program:               |
|  Foxp3, Fgl2, BACH2, BLIMP-1            |
+------------------------------------------+
```

### Supporting Evidence from Literature

1. **Tregs can become Tfr cells** by upregulating Bcl6 and CXCR5. Tfr cells localize to germinal centers and suppress B cell responses (Linterman et al., *Nature Medicine* 2011; Chung et al., *J Exp Med* 2011).

2. **Tfr cells are inherently unstable** because Bcl6 represses CD25, undermining the IL-2 signaling required to maintain Foxp3. Under homeostasis ~1% of Tregs lose Foxp3, but Tfr cells are particularly prone to this conversion (reviewed in *J Leukocyte Biology* 2024).

3. **Ex-Foxp3 Tfr cells become pathogenic Tfh-like cells** that promote rather than suppress germinal center reactions. Foxp3 conditional-knockout mice show increased Tfh cells, GC B cells, and autoantibody production (BBRC 2019).

4. **FOXP3 exon 2 controls Treg stability** — its loss leads to autoimmunity through acquisition of effector programs (Bhela et al., *Science Immunology* 2022).

5. **Nrp1+ thymic Tregs may preferentially take the Tfh route** upon Foxp3 loss. Our exTregs retain high Nrp1 (12,410) and show follicular features, consistent with this model.

### Key Insights from This Dataset

- **Polyfunctional cytokine production is robust**: IL-21, IL-4, IFN-gamma, and IL-2 all survive apeglm shrinkage with clean per-replicate consistency. This simultaneous multi-lineage cytokine production suggests global derepression of effector loci rather than clean lineage conversion.

- **c-MAF is the key transcriptional node**: At 61,705 normalized counts, c-MAF is extraordinarily elevated in exTregs. In the absence of Foxp3-mediated repression, c-MAF likely drives IL-21 and cooperates with RORgt and other TFs to enable the polyfunctional state. Bcl6 RNA differences are modest across all comparisons, suggesting c-MAF rather than Bcl6 as the primary transcriptional driver.

- **T-bet-independent IFN-gamma**: exTregs produce IFN-gamma robustly (vs Stable_Treg) without T-bet upregulation, consistent with c-MAF- or STAT4-dependent IFN-gamma.

- **GC-associated gene activation** (Pax5, Serpina9, Sostdc1, Pou2af1) without B cell lineage commitment indicates shared transcriptional circuitry between Tfh cells and GC B cells being co-opted in ex-Foxp3 cells.

- **Flow cytometry validation aligns with robust hits**: Luke's flow data confirming the Tfh signal (CXCR5, PD-1) corresponds to genes that survive shrinkage. Protein-level data can address Bcl6, where RNA differences are modest but post-translational regulation may produce protein-level differences.

---

## Output Files

### DESeq2 — Unshrunken (`deseq2/`)
- `DEG_*.csv` — Full differential expression results with gene names (6 files)
- `normalized_counts.csv` — DESeq2-normalized count matrix
- `key_genes_expression.csv` — Curated table of Treg/exhaustion/function genes
- `pca_plot.pdf` — PCA of all samples
- `sample_distance_heatmap.pdf` — Sample clustering
- `volcano_*.pdf` — Labeled volcano plots (6 files)
- `upset_DEGs.pdf` — UpSet plot of shared DEGs
- `heatmap_top_DEGs.pdf` — Top DEGs across comparisons
- `heatmap_treg_exhaustion_genes.pdf` — Focused heatmap of Treg/exhaustion/function genes
- `venn_characteristic_genes.pdf` — Venn diagram of cell type signatures
- `GO_characteristic_*.csv` — GO enrichment for each cell type's signature genes

### DESeq2 — apeglm Shrunken (`deseq2_shrinkage/`)
- `DEG_shrunk_*.csv` — Shrunken DEG results (baseMean >= 50) for all 6 comparisons
- `MA_*.pdf` — Side-by-side MA plots (unshrunken vs shrunken)
- `volcano_shrunk_*.pdf` — Volcano plots using shrunken log2 fold changes
- `shrinkage_scatter_*.pdf` — Scatter plots showing shrinkage effect by expression level
- `shrinkage_summary.csv` — DEG count comparison table

### Pathway Analysis (`pathway_analysis/`)
Per comparison subdirectories containing:
- `GO_BP.csv` / `GO_BP_dotplot.pdf` — GO Biological Process ORA
- `KEGG.csv` / `KEGG_dotplot.pdf` — KEGG pathway ORA
- `GSEA_hallmark.csv` / `GSEA_hallmark_barplot.pdf` — MSigDB Hallmark GSEA
- `GSEA_immunologic.csv` — C7 ImmuneSigDB GSEA
- `all_hallmark_results.csv` — Cross-comparison Hallmark summary

---

## References

1. Linterman MA et al. Foxp3+ follicular regulatory T cells control the germinal center response. *Nature Medicine* 17, 975-982 (2011).
2. Chung Y et al. Follicular regulatory T cells expressing Foxp3 and Bcl-6 suppress germinal center reactions. *Nature Medicine* 17, 983-988 (2011).
3. Bhela S et al. FOXP3 exon 2 controls Treg stability and autoimmunity. *Science Immunology* 7, eabo5407 (2022).
4. Review: Stability and plasticity of regulatory T cells in health and disease. *Journal of Leukocyte Biology* 116, 33-53 (2024).
5. Dysregulation of humoral immunity in Foxp3 conditional-knockout mice. *Biochem Biophys Res Commun* (2019).
6. IL-21 restricts T follicular regulatory T cell proliferation through Bcl-6 mediated inhibition of responsiveness to IL-2. *Nature Communications* 8, 14647 (2017).
7. Pathological conversion of regulatory T cells is associated with loss of allotolerance. *Scientific Reports* 8, 7059 (2018).
8. Zhu A, Ibrahim JG, Love MI. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics* 35, 2084-2092 (2018).

---

*Analysis performed March 5-6, 2026*
*Pipeline: STAR -> featureCounts -> DESeq2 (apeglm shrinkage) -> clusterProfiler/fgsea*
