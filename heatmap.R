# =============================================================================
# Clustered heatmap of the high-confidence DEGs
# Uses variance-stabilising transformed (VST) counts, Z-scored per gene,
# with hierarchical clustering of rows (genes) and columns (samples).
# Inputs:
#   cleaned_counts.csv, metadata_counts.csv
#   DEGs_86_final.csv   (curated high-confidence DEGs; column 'gene')
# =============================================================================

library(DESeq2)
library(pheatmap)

counts <- read.csv("cleaned_counts.csv", row.names = 1, check.names = FALSE)
meta   <- read.csv("metadata_counts.csv")
rownames(meta)   <- toupper(meta$samples)
colnames(counts) <- toupper(colnames(counts))
counts <- counts[, rownames(meta)]
meta$condition <- factor(tolower(meta$condition), levels = c("control", "diseased"))

dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData   = meta,
                              design    = ~ condition)
dds <- DESeq(dds)
vsd <- vst(dds)

# Restrict to the curated high-confidence DEGs present in the matrix
degs <- read.csv("DEGs_86_final.csv")$gene
degs <- intersect(degs, rownames(assay(vsd)))

mat <- assay(vsd)[degs, ]
mat <- t(scale(t(mat)))            # per-gene Z-score

ann <- data.frame(Condition = meta$condition)
rownames(ann) <- rownames(meta)

pheatmap(mat,
         cluster_rows        = TRUE,
         cluster_cols        = TRUE,
         clustering_method   = "complete",
         annotation_col      = ann,
         show_rownames       = FALSE,
         main                = "High-confidence DEGs (HCM vs. Control)")
