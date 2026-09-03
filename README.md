# Hypertrophic Cardiomyopathy (HCM) Transcriptomic Analysis

**Integrative Transcriptomics Identifies Non-sarcomeric Dysregulated Ion Homeostasis and Innate Immunity Genes as Therapeutic Targets in Hypertrophic Cardiomyopathy (HCM)**

This repository contains the code and curated result files for the RNA-seq, network-biology, and target-based drug knowledge-mining analysis reported in the manuscript.

---

## Overview

Hypertrophic Cardiomyopathy (HCM) is the most prevalent monogenic cardiac disorder, with no definitive cure and treatment currently limited to symptomatic management. This project applies an RNA-seq–based workflow that integrates transcriptomics, network biology, and target-based drug knowledge mining to uncover non-sarcomeric molecular mechanisms and prioritise candidate therapeutic targets.

## Objectives

- Identify differentially expressed genes (DEGs) in HCM myocardial tissue
- Elucidate dysregulated biological pathways (GO / KEGG)
- Construct a protein–protein interaction (PPI) network and identify hub genes
- Mine curated drug–target knowledge (DrugBank) for candidate compounds against the prioritized hub genes

## Dataset

- **Source:** NCBI Sequence Read Archive (SRA)
- **BioProject:** PRJNA961804
- **Sample type:** Human myocardial tissue (HCM vs. control)
- **Reference genome:** GRCh38

Raw sequencing data are publicly available from SRA; they are **not** redistributed here. See `preprocessing.sh` for how to download and process them.

---

## Workflow and tools

| Stage | Tool(s) |
|-------|---------|
| Data retrieval | SRA Toolkit (`prefetch`, `fasterq-dump`) |
| Quality control | FastQC, MultiQC |
| Trimming | Trimmomatic (`SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:36`; no adapter clip) |
| Alignment & indexing | Bowtie2 (reference: GRCh38) |
| Quantification | featureCounts |
| Differential expression | DESeq2 (R) — **&#124;log2FC&#124; ≥ 3, padj < 0.05** |
| Visualisation | ggplot2, pheatmap, factoextra |
| Functional enrichment | ShinyGO, Enrichr (GO / KEGG) |
| PPI network | Metascape (integrating STRING, BioGRID, OmniPath, InWeb_IM) |
| Modules & hubs | Cytoscape — MCODE (clusters) and CytoHubba (MCC ranking) |
| Drug knowledge mining | DrugBank (target-based) |

---

## Repository contents

```
.
├── README.md                     # this file
├── LICENSE                       # MIT license
├── preprocessing.sh              # CLI steps: SRA → FASTQ → QC → trim → align → count
├── install_packages.R            # installs the R/Bioconductor packages used here
├── deseq2_analysis.R             # DESeq2 differential expression
├── volcano_plot.R                # volcano plot of DEGs
├── heatmap.R                     # clustered heatmap of high-confidence DEGs
├── enrichment_plot.R             # GO enrichment bar plot
├── kegg_gene_chord.R             # KEGG pathway–gene chord diagram
├── cleaned_counts.csv            # featureCounts matrix (genes × 8 samples)
├── metadata_counts.csv           # sample sheet (samples, condition)
├── pathways_enrichr.csv          # Enrichr KEGG enrichment export
├── go_enrichment.csv             # GO enrichment export (per-cluster representative terms)
└── DEGs_86_final.csv             # 86 high-confidence DEGs (|log2FC| ≥ 3; = DESeq2_significant_geneslog3)
```

## Input data files (included)

The processed data needed to reproduce every result are committed in this repository (the *raw* FASTQ files are not — they are available from SRA, BioProject PRJNA961804):

| File | Description | Used by |
|------|-------------|---------|
| `cleaned_counts.csv` | featureCounts matrix, gene symbols × 8 samples | `deseq2_analysis.R`, `heatmap.R` |
| `metadata_counts.csv` | sample sheet: `samples`, `condition` (`control` / `diseased`) | `deseq2_analysis.R`, `heatmap.R` |
| `pathways_enrichr.csv` | Enrichr KEGG enrichment (`Term`, `Genes`, `Adjusted.P.value`, …) | `kegg_gene_chord.R` |
| `go_enrichment.csv` | GO enrichment (`Category`, `Term`, `PValue`, `Genes`) | `enrichment_plot.R` |

The dataset is 8 samples: 6 HCM (diseased) and 2 control.

---

## Reproducing the analysis

```bash
# 1. Preprocess raw reads (edit accession list inside the script first)
bash preprocessing.sh          # produces cleaned_counts.csv via featureCounts
```

```r
# 2. Install R dependencies (once)
source("install_packages.R")

# 3. Differential expression (uses the included count matrix + metadata)
source("deseq2_analysis.R")    # writes DESeq2_results.csv and the three threshold sets

# 4. Figures
source("volcano_plot.R")
source("heatmap.R")
source("enrichment_plot.R")
source("kegg_gene_chord.R")
```

### DEG thresholds

`deseq2_analysis.R` applies progressively stricter fold-change cut-offs at `padj < 0.05`, reproducing the counts reported in the manuscript:

| Threshold | DEGs |
|-----------|------|
| &#124;log2FC&#124; ≥ 1 (discovery set) | 989 |
| &#124;log2FC&#124; ≥ 2 | 409 |
| **&#124;log2FC&#124; ≥ 3 (high-confidence)** | **86** (24 up, 62 down) |

The |log2FC| ≥ 3 output (`DESeq2_significant_geneslog3.csv`) is the 86 high-confidence DEG set used for all downstream network, enrichment, hub-gene, and drug-mining steps, and is also provided as **`DEGs_86_final.csv`** (identical to Supplementary Table S2).

---

## Key results

- **86 high-confidence DEGs** (24 up-, 62 down-regulated)
- **Dysregulated pathways:** mineral absorption / zinc-ion homeostasis (metallothioneins MT1A, MT1M) and neutrophil extracellular trap (NET) formation (innate-immune dysregulation)
- **Hub genes (CytoHubba MCC, top 10):** FPR1, FPR2, CXCR1, CCR1, SSTR5, EDN2, EDN3, S100A8, S100A9, TLR2
- **Candidate agents (DrugBank):** zinc salts (zinc chloride, zinc acetate, zinc sulfate) as a mineral-supplementation lead; formyl-peptide-receptor antagonists rebamipide and nedocromil, and the S100A9 modulator tasquinimod, among others

## Biological significance

The analysis complements the sarcomere-centric view of HCM by highlighting non-sarcomeric processes, impaired metallothionein-mediated zinc homeostasis and activated innate immunity/NET formation, supporting a dual-pathology model and nominating targets for experimental validation.

---

## Citation

If you use this code or the curated result files, please cite the associated manuscript (Current Genomics, under revision) and this repository. The BioProject (PRJNA961804) should be cited for the primary data.

## Author

Muhammad Haris

## Disclaimer

This study is computational and its findings require experimental validation.
