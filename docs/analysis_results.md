# CD4 T Cell Subset RNA-seq Analysis: exTregs Acquire a Pathogenic Tfh-like Program

## Overview

Bulk RNA-seq analysis of four FACS-sorted CD4 T cell populations from mouse lymphoid tissue reveals that ex-regulatory T cells (exTregs) undergo a dramatic transcriptional reprogramming upon Foxp3 loss. Rather than reverting to a naive state or committing to a single effector lineage, exTregs acquire a polyfunctional, germinal center-associated program dominated by T follicular helper (Tfh) cell features, with concurrent activation of Th1 and Th17 gene modules. These findings are consistent with a Treg → Tfr → pathogenic ex-Foxp3 conversion pathway described in recent literature, and extend it by demonstrating simultaneous multi-lineage effector cytokine capacity.

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
- **Sequencing depth**: 60–104M alignment pairs per sample

## Analysis Pipeline

| Step | Script | Tool | Key Parameters |
|---|---|---|---|
| 1. QC & trimming | `01_star_align.sh` | STAR 2.7.x | Aligned to GRCm39 mouse genome |
| 2. Quantification | `02_featurecounts.sh` | featureCounts 2.0.6 | Paired-end (`-p --countReadPairs`), unstranded (`-s 0`), `--extraAttributes gene_name` |
| 3. Differential expression | `03_deseq2_analysis.R` | DESeq2 1.48.2 | Paired design: `~ mouse + cell_type`, Never_Treg as reference |
| 4. Pathway analysis | `04_pathway_analysis.R` | clusterProfiler, fgsea, msigdbr | GO BP, KEGG (ORA); Hallmark & C7 immunologic (GSEA) |
| 5. Cell type signatures | `05_venn_cell_type_pathways.R` | ggVennDiagram, clusterProfiler | Characteristic genes UP vs ≥2 of 3 other populations |

### Strandedness Verification

STAR ReadsPerGene.out.tab showed near-equal counts in sense and antisense columns across all samples (~50/50 split), confirming an **unstranded** library. featureCounts was run with `-s 0` accordingly.

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

**16,075 genes** passed filtering (≥10 counts in ≥3 samples).

---

## Differential Expression Results

All 6 pairwise comparisons were performed using DESeq2 with a paired design (`~ mouse + cell_type`) to account for mouse-to-mouse variation. Significance thresholds: padj < 0.05, |log2FC| > 1.

| Comparison | Up | Down | Total DEGs |
|---|---|---|---|
| exTreg vs Never_Treg | 485 | 178 | 663 |
| Stable_Treg vs Never_Treg | 804 | 1,154 | 1,958 |
| Recent_Treg vs Never_Treg | 310 | 66 | 376 |
| Stable_Treg vs exTreg | 414 | 1,141 | 1,555 |
| Recent_Treg vs exTreg | 151 | 101 | 252 |
| Recent_Treg vs Stable_Treg | 968 | 214 | 1,182 |

### QC Validation

- **PCA**: Clear separation of all four populations. PC1 separates Treg-identity (Foxp3+) populations from Foxp3-low/negative populations; PC2 distinguishes exTreg from Never_Treg.
- **Sample distance heatmap**: Replicates cluster tightly by cell type.
- **Output files**: `deseq2/pca_plot.pdf`, `deseq2/sample_distance_heatmap.pdf`

---

## Sanity Check: Treg Identity Genes

Foxp3 and core Treg markers validate the sort purity and population definitions.

| Gene | Protein | Never_Treg | exTreg | Stable_Treg | Recent_Treg | Biology |
|---|---|---|---|---|---|---|
| **Foxp3** | FOXP3 | 101 | 191 | **69,700** | 27,769 | Master Treg TF |
| **Il2ra** | CD25 | 535 | 5,126 | **86,237** | 41,930 | IL-2 receptor, Treg survival |
| **Ikzf2** | Helios | 2,503 | 39,172 | **79,792** | 14,169 | Thymic Treg marker |
| **Ctla4** | CTLA-4 | 1,117 | 7,918 | **27,812** | 18,772 | Co-inhibitory receptor |
| **Nrp1** | Neuropilin-1 | 413 | 12,410 | **18,233** | 2,765 | Thymic Treg origin marker |
| **Lrrc32** | GARP | 58 | 980 | **47,187** | 12,892 | TGF-beta activation |
| **Itgae** | CD103 | 31 | 328 | **9,466** | 1,548 | Tissue-resident Treg |
| **Nt5e** | CD73 | 1,959 | 16,153 | **29,102** | 10,238 | Adenosine pathway (suppression) |
| **Entpd1** | CD39 | 248 | 738 | **3,471** | 1,805 | Adenosine pathway (suppression) |

**Key observations:**
- Stable_Treg has the strongest Treg signature across all markers
- Recent_Treg shows intermediate levels, consistent with a newly committed state
- exTreg retains **high Nrp1** (12,410), confirming thymic Treg origin even after Foxp3 loss
- exTreg retains **high Ikzf2/Helios** (39,172), suggesting recent Treg identity

---

## Key Finding: exTregs Acquire a Dominant Tfh-like Program

### Exhaustion Markers

exTregs express the highest levels of multiple inhibitory receptors:

| Gene | Protein | Never_Treg | exTreg | Stable_Treg | Recent_Treg | Fold vs next |
|---|---|---|---|---|---|---|
| **Pdcd1** | PD-1 | 92 | **3,309** | 397 | 409 | 8.1x |
| **Lag3** | LAG-3 | 69 | **1,654** | 800 | 815 | 2.0x |
| **Tigit** | TIGIT | 284 | **6,020** | 4,894 | 1,180 | 1.2x |
| **Havcr2** | TIM-3 | 0 | **63** | 23 | 41 | 1.5x |
| **Tox** | TOX | 3,603 | **7,500** | 6,031 | 4,586 | 1.2x |
| **Tox2** | TOX2 | 421 | **2,775** | 288 | 145 | 6.6x |

While PD-1 and exhaustion markers are elevated, the functional profile of exTregs reveals that these cells are not simply "exhausted" — they have acquired active effector programs.

### T Follicular Helper (Tfh) Program — Strongest Signal

The most striking finding is that exTregs express a near-complete Tfh transcriptional program:

| Gene | Protein | Never_Treg | exTreg | Stable_Treg | Recent_Treg | Role |
|---|---|---|---|---|---|---|
| **Bcl6** | BCL-6 | 2,181 | **4,208** | 1,986 | 1,503 | Master Tfh TF |
| **Cxcr5** | CXCR5 | 1,128 | **6,875** | 5,310 | 1,755 | GC homing receptor |
| **Il21** | IL-21 | 648 | **4,529** | 0 | 27 | Primary Tfh cytokine |
| **Ascl2** | ASCL2 | 37 | **380** | 0 | 1 | Tfh-specific TF |
| **Sh2d1a** | SAP | 2,445 | **6,411** | 2,458 | 3,305 | SLAM family signaling |
| **Slamf6** | Ly108 | 16,480 | **35,794** | 12,897 | 17,484 | GC T-B interaction |
| **Tox2** | TOX2 | 421 | **2,775** | 288 | 145 | Tfh differentiation |
| **Pou2af1** | OBF-1/BOB.1 | 1,383 | **2,463** | 89 | 556 | GC-associated co-activator |
| **Maf** | c-MAF | 10,761 | **61,705** | 29,378 | 32,006 | IL-21 driver, Tfh/Th17 |

**IL-21** is the defining Tfh cytokine and is almost exclusively produced by exTregs (4,529 vs 0–648 in all others). Combined with Bcl6, CXCR5, ASCL2, and SAP, this represents a functional Tfh program.

**c-MAF** at 61,705 normalized counts is extraordinarily high — nearly 6x the level in Never_Treg and 2x Stable_Treg. c-MAF is a key driver of IL-21 transcription and cooperates with Bcl6 in Tfh differentiation.

### Th1-like Properties — IFN-gamma Without T-bet

| Gene | Protein | Never_Treg | exTreg | Stable_Treg | Recent_Treg |
|---|---|---|---|---|---|
| **Ifng** | IFN-gamma | 204 | **1,255** | 28 | 84 |
| **Tnf** | TNF-alpha | 866 | **2,709** | 1,333 | 1,450 |
| **Stat4** | STAT4 | 4,969 | **8,793** | 5,248 | 5,110 |
| **Cxcr3** | CXCR3 | 917 | **5,400** | 2,332 | 1,420 |
| **Il18r1** | IL-18Ra | 6,402 | **8,665** | 4,444 | 3,031 |
| **Tbx21** | T-bet | 204 | 175 | 23 | 26 |

exTregs produce IFN-gamma (6x vs Never_Treg) and TNF despite **no upregulation of T-bet**. This T-bet-independent IFN-gamma production is consistent with a Tfh-associated mechanism where Bcl6 and c-Maf can drive IFN-gamma in the absence of classical Th1 commitment. CXCR3 co-expression with CXCR5 may indicate a CXCR3+CXCR5+ "Tfh1" phenotype described in autoimmune and chronic inflammatory settings.

### Th17-like Properties — Partial Program

| Gene | Protein | Never_Treg | exTreg | Stable_Treg | Recent_Treg |
|---|---|---|---|---|---|
| **Rorc** | RORgt | 13 | **653** | 126 | 562 |
| **Il17a** | IL-17A | 0 | **113** | 0 | 12 |
| **Il23r** | IL-23R | 61 | **121** | 0 | 30 |
| **Ccr6** | CCR6 | 225 | **2,189** | 1,674 | 1,377 |
| **Hif1a** | HIF-1alpha | 9,323 | **22,605** | 18,608 | 13,608 |
| **Il1r1** | IL-1R1 | 182 | **926** | 281 | 268 |
| **Bhlhe40** | DEC1 | 291 | **2,724** | 429 | 332 |

RORgt is activated (50x vs Never_Treg) and exTregs are the only population producing detectable IL-17A (113 counts). HIF-1alpha (22,605) is the highest across all populations — HIF-1a is known to shift the Treg/Th17 balance toward Th17 by degrading Foxp3 and activating RORgt.

### Effector Cytokine Polyfunctionality

exTregs are the only population simultaneously producing all of the following:

| Cytokine | exTreg | All others | Lineage association |
|---|---|---|---|
| IFN-gamma | 1,255 | 28–204 | Th1 |
| TNF | 2,709 | 866–1,450 | Th1/inflammatory |
| IL-21 | 4,529 | 0–648 | Tfh |
| IL-4 | 387 | 0–95 | Th2/Tfh |
| IL-17A | 113 | 0–12 | Th17 |
| IL-2 | 630 | 17–600 | T cell activation |

This polyfunctionality is unusual and suggests that exTregs are in a **dysregulated state with simultaneous derepression of multiple effector loci**, rather than undergoing clean lineage conversion.

---

## Germinal Center-Associated Genes: Not B Cell Contamination

Several genes normally associated with B cells were found to be highly expressed specifically in exTregs, raising the question of sort contamination:

| Gene | exTreg | Others | Normal expression |
|---|---|---|---|
| **Pax5** | 1,173 | 0–2 | B cell master TF |
| **Serpina9** | 1,079 | 0–5 | GC B cells |
| **Sostdc1** | 7,245 | 0–21 | GC-associated |
| **Cd22** | 488 | 8–110 | B cell surface |

### Evidence Against B Cell Contamination

1. **CD19 and CD20 (Ms4a1) are absent** — filtered out for having too few counts. These are the most highly expressed B cell surface markers and would be easily detected even at low contamination levels.
2. **Ebf1 and Blnk are absent** — core B cell transcription factors and signaling molecules not detected.
3. **Immunoglobulin genes do not track with exTreg** — Ighm is highest in Stable_Treg (53,969), not exTreg (16,282). Ig genes represent only 0.04–0.14% of the library across all populations (ambient RNA).
4. **T cell identity is rock solid** — CD3e (42,226), CD4 (120,248), TCR chains all at very high levels in exTreg. These are unambiguously T cells.
5. **Core Pax5 transcriptional targets are not activated** — if Pax5 were functionally driving a B cell program, CD19, CD20, Blnk, and Ebf1 would be its first targets. Their absence indicates Pax5 expression without functional B cell commitment.

### The Tfh–GC Connection Explains These Genes

The expression of Pax5, Serpina9, and Pou2af1 in exTregs is better explained by their shared expression in germinal center reactions. Pou2af1 (OBF-1/BOB.1) is a known Tfh-associated gene, and its co-expression with Pax5 (Spearman r = 0.81 with Sostdc1, r = 0.75 with Serpina9) suggests that exTregs are activating a **germinal center-associated transcriptional network** rather than undergoing B cell transdifferentiation.

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
- **KRAS signaling up** (NES=1.60) — MAPK signaling

Recent_Treg vs Stable_Treg revealed downregulation of:
- **IL6-JAK-STAT3 signaling** (NES=-1.95) — less inflammatory signaling in Recent
- **TGF-beta signaling** (NES=-1.73) — less suppressive program
- **IFN-gamma response** (NES=-1.57) — less effector-like

### GO Biological Process Enrichment (Treg/Immune-related terms)

**exTreg vs Never_Treg (22 immune GO terms):**
- T cell activation (p.adj=7.6e-04)
- Alpha-beta T cell activation and differentiation
- Negative regulation of T cell mediated immunity (p.adj=1.5e-03)
- Tolerance induction (p.adj=6.2e-03)
- T cell tolerance induction (p.adj=3.3e-02)

**Stable_Treg vs exTreg (31 immune GO terms):**
- Regulation of cytokine production in immune response (p.adj=1.5e-03)
- Regulation of T cell activation (p.adj=7.8e-03)
- Regulation of Th17 cell differentiation (p.adj=8.9e-03)
- T-helper 17 type immune response (p.adj=2.9e-02)

### Cell Type Signature Analysis (Venn Diagram)

Characteristic genes per population (UP vs at least 2 of 3 other populations):

| Cell Type | Total | Unique | Top Pathways |
|---|---|---|---|
| **Never Treg** | 103 | 61 | NK cell cytotoxicity, cell killing, lymphocyte differentiation |
| **exTreg** | 150 | 91 | Inflammatory response regulation, immune effector regulation, macrophage activation |
| **Stable Treg** | 293 | 244 | T cell activation/differentiation, hemopoiesis regulation, IFN-gamma response |
| **Recent Treg** | 103 | 20 | Leukocyte homeostasis, lymphocyte proliferation regulation |

Stable_Treg has the most unique genes (244), consistent with it being the most transcriptionally distinct population with a fully committed suppressive program. Recent_Treg has the fewest unique genes (20), consistent with a transitional state sharing features with multiple populations.

---

## Proposed Model: Treg → Tfr → Pathogenic ex-Foxp3 Conversion

Based on these transcriptomic data and published literature, we propose the following model for exTreg differentiation:

```
Thymic Treg (Foxp3+, Nrp1+, Helios+)
       │
       ▼
Follicular Regulatory T cell (Tfr)
(Foxp3+, Bcl6+, CXCR5+, PD-1+)
       │
       │  Foxp3 destabilization:
       │  - HIF-1a upregulation
       │  - IL-21 autocrine/paracrine
       │  - Bcl6 represses CD25 → loss of IL-2 signaling
       │  - c-Maf drives IL-21 production
       ▼
ex-Treg (Foxp3-, Nrp1+, Helios+)
┌──────────────────────────────────┐
│  Dominant Tfh program:           │
│  Bcl6, CXCR5, IL-21, SAP, ASCL2 │
│                                  │
│  Exhaustion overlay:             │
│  PD-1, LAG-3, TIGIT, TOX/TOX2   │
│                                  │
│  Polyfunctional cytokines:       │
│  IFN-g, TNF, IL-17A, IL-4, IL-21│
│                                  │
│  GC-associated genes:            │
│  Pou2af1, Pax5, Serpina9         │
└──────────────────────────────────┘
```

### Supporting Evidence from Literature

1. **Tregs can become Tfr cells** by upregulating Bcl6 and CXCR5. Tfr cells localize to germinal centers and suppress B cell responses (Linterman et al., *Nature Medicine* 2011; Chung et al., *J Exp Med* 2011).

2. **Tfr cells are inherently unstable** because Bcl6 represses CD25, undermining the IL-2 signaling required to maintain Foxp3. Under homeostasis ~1% of Tregs lose Foxp3, but Tfr cells are particularly prone to this conversion (reviewed in *J Leukocyte Biology* 2024).

3. **Ex-Foxp3 Tfr cells become pathogenic Tfh-like cells** that promote rather than suppress germinal center reactions. Foxp3 conditional-knockout mice show increased Tfh cells, GC B cells, and autoantibody production (BBRC 2019).

4. **FOXP3 exon 2 controls Treg stability** — its loss leads to autoimmunity through acquisition of effector programs (Bhela et al., *Science Immunology* 2022).

5. **Nrp1+ thymic Tregs may preferentially take the Tfh route** upon Foxp3 loss, while Nrp1- peripherally-derived Tregs are more susceptible to Th17/Th1 conversion. Our exTregs retain high Nrp1 (12,410) and show a dominant Tfh program, consistent with this model.

### Novel Aspects of These Data

- **Simultaneous multi-lineage cytokine production**: Most published studies report exTregs becoming Th1-like OR Th17-like OR Tfh-like. Our exTregs produce IFN-gamma, IL-17A, and IL-21 simultaneously, suggesting a globally dysregulated state rather than clean lineage conversion.

- **Extraordinary c-MAF expression** (61,705 counts): c-MAF is a critical driver of IL-21 and cooperates with both Bcl6 (Tfh) and RORgt (Th17). The ~6x elevation over Never_Treg may be a key node enabling the polyfunctional state.

- **T-bet-independent IFN-gamma**: exTregs produce IFN-gamma without T-bet upregulation, suggesting Bcl6/c-Maf-driven or STAT4-dependent IFN-gamma that is mechanistically distinct from classical Th1 differentiation.

- **Germinal center-associated gene activation** (Pax5, Serpina9, Pou2af1) without B cell lineage commitment, indicating shared transcriptional circuitry between Tfh cells and GC B cells being co-opted in ex-Foxp3 cells.

---

## Treg Function and Stability Genes

| Gene | Protein | Never_Treg | exTreg | Stable_Treg | Recent_Treg | Interpretation |
|---|---|---|---|---|---|---|
| Il10 | IL-10 | 0 | 142 | **2,061** | 289 | Suppressive cytokine; highest in Stable |
| Tgfb1 | TGF-beta1 | 2,937 | 3,462 | **6,661** | 2,645 | Suppressive cytokine; highest in Stable |
| Fgl2 | Fibrinogen-like 2 | 308 | 187 | **4,933** | 1,235 | Treg effector molecule |
| Prdm1 | BLIMP-1 | 677 | 859 | **2,983** | 1,492 | Effector Treg differentiation |
| Gzmb | Granzyme B | 26 | 165 | **417** | 425 | Treg-mediated killing |
| Bcl2 | BCL-2 | 15,959 | 9,103 | 12,758 | **25,868** | Survival; highest in Recent |
| Maf | c-MAF | 10,761 | **61,705** | 29,378 | 32,006 | Multi-lineage TF |
| Rorc | RORgt | 13 | **653** | 126 | 562 | Th17 master TF |
| Bach2 | BACH2 | **10,314** | 2,351 | 7,246 | 7,651 | Represses effector programs |

Notable: **BACH2** is dramatically reduced in exTregs (2,351 vs 10,314 in Never_Treg). BACH2 is a transcriptional repressor that maintains T cell quiescence and prevents effector differentiation. Its loss in exTregs likely contributes to the derepression of multiple effector programs.

---

## Master Transcription Factor Landscape

| TF | Never_Treg | exTreg | Stable_Treg | Recent_Treg | Lineage |
|---|---|---|---|---|---|
| Foxp3 | 101 | 191 | **69,700** | 27,769 | Treg |
| Bcl6 | 2,181 | **4,208** | 1,986 | 1,503 | Tfh |
| Tbx21 | 204 | 175 | 23 | 26 | Th1 |
| Gata3 | 2,504 | 3,680 | **5,268** | 2,465 | Th2/Treg |
| Rorc | 13 | **653** | 126 | 562 | Th17 |
| Eomes | **447** | 202 | 136 | 208 | Cytotoxic/Th1 |
| Maf | 10,761 | **61,705** | 29,378 | 32,006 | Tfh/Th17/Th2 |
| Prdm1 | 677 | 859 | **2,983** | 1,492 | Effector Treg |
| Tox | 3,603 | **7,500** | 6,031 | 4,586 | Exhaustion |
| Tox2 | 421 | **2,775** | 288 | 145 | Tfh/exhaustion |
| Batf | 976 | 2,089 | **2,540** | 1,009 | AP-1 family |
| Irf4 | 1,496 | 5,154 | **6,746** | 4,184 | Effector differentiation |
| Bach2 | **10,314** | 2,351 | 7,246 | 7,651 | Quiescence |
| Bhlhe40 | 291 | **2,724** | 429 | 332 | Th1/pathogenic |
| Hif1a | 9,323 | **22,605** | 18,608 | 13,608 | Th17/glycolysis |
| Ikzf2 | 2,503 | **39,172** | 79,792 | 14,169 | Treg (Helios) |

The exTreg TF landscape is dominated by **c-MAF >> Bcl6 + RORgt + Bhlhe40 + HIF-1a**, with low T-bet and absent Foxp3. This is not a classical Th1, Th2, or Th17 profile — it is a unique combination consistent with a GC-associated, polyfunctional, pathogenic ex-Treg state.

---

## Chemokine Receptor Profile — Trafficking Fate

| Receptor | Never_Treg | exTreg | Stable_Treg | Recent_Treg | Homing |
|---|---|---|---|---|---|
| **CXCR5** | 1,128 | **6,875** | 5,310 | 1,755 | B cell follicles / GC |
| **CXCR3** | 917 | **5,400** | 2,332 | 1,420 | Inflamed tissue (Th1) |
| **CCR6** | 225 | **2,189** | 1,674 | 1,377 | Mucosal / Th17 |
| CCR7 | **6,960** | 2,748 | 4,079 | 5,118 | Lymph node T zone |
| S1PR1 | **14,400** | 7,658 | 9,746 | 8,310 | Lymph node egress |
| CCR4 | 1,060 | 3,493 | **9,096** | 4,326 | Skin / Treg homing |

exTregs show **reduced CCR7 and S1PR1** (less T zone/egress), with **elevated CXCR5 + CXCR3 + CCR6**. The CXCR5+CXCR3+ co-expression pattern is characteristic of "Tfh1" cells described in autoimmune and chronic inflammatory settings, positioning exTregs to traffic to both B cell follicles and sites of inflammation.

---

## Output Files

### DESeq2 (`deseq2/`)
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

### Pathway Analysis (`pathway_analysis/`)
Per comparison subdirectories containing:
- `GO_BP.csv` / `GO_BP_dotplot.pdf` — GO Biological Process ORA
- `KEGG.csv` / `KEGG_dotplot.pdf` — KEGG pathway ORA
- `GSEA_hallmark.csv` / `GSEA_hallmark_barplot.pdf` — MSigDB Hallmark GSEA
- `GSEA_immunologic.csv` — C7 ImmuneSigDB GSEA
- `all_hallmark_results.csv` — Cross-comparison Hallmark summary

---

## References

1. Linterman MA et al. Foxp3+ follicular regulatory T cells control the germinal center response. *Nature Medicine* 17, 975–982 (2011).
2. Chung Y et al. Follicular regulatory T cells expressing Foxp3 and Bcl-6 suppress germinal center reactions. *Nature Medicine* 17, 983–988 (2011).
3. Bhela S et al. FOXP3 exon 2 controls Treg stability and autoimmunity. *Science Immunology* 7, eabo5407 (2022).
4. Review: Stability and plasticity of regulatory T cells in health and disease. *Journal of Leukocyte Biology* 116, 33–53 (2024).
5. Dysregulation of humoral immunity in Foxp3 conditional-knockout mice. *Biochem Biophys Res Commun* (2019).
6. IL-21 restricts T follicular regulatory T cell proliferation through Bcl-6 mediated inhibition of responsiveness to IL-2. *Nature Communications* 8, 14647 (2017).
7. Pathological conversion of regulatory T cells is associated with loss of allotolerance. *Scientific Reports* 8, 7059 (2018).

---

*Analysis performed March 5, 2026*
*Pipeline: STAR → featureCounts → DESeq2 → clusterProfiler/fgsea*
