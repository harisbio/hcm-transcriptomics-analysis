# =============================================================================
# Differential expression analysis (DESeq2)
# HCM myocardium vs. control  |  reference genome GRCh38
# =============================================================================
# Inputs:
#   cleaned_counts.csv   - featureCounts matrix (genes x samples; gene symbols as row names)
#   metadata_counts.csv  - sample sheet: columns 'samples', 'condition' (control/diseased)
# Outputs:
#   DESeq2_results.csv               - full results table (all genes)
#   DESeq2_significant_geneslog1.csv - |log2FC| >= 1 & padj < 0.05  (discovery set)
#   DESeq2_significant_geneslog2.csv - |log2FC| >= 2 & padj < 0.05
#   DESeq2_significant_geneslog3.csv - |log2FC| >= 3 & padj < 0.05  (86 high-confidence DEGs)
# =============================================================================

library(DESeq2)

counts <- read.csv("cleaned_counts.csv", row.names = 1, check.names = FALSE)
meta   <- read.csv("metadata_counts.csv")

# Match samples by name (case-insensitive) and order counts to the metadata
rownames(meta)   <- toupper(meta$samples)
colnames(counts) <- toupper(colnames(counts))
counts <- counts[, rownames(meta)]                       # reorder to metadata
stopifnot(all(colnames(counts) == rownames(meta)))

# Control is the reference level; positive log2FC = up-regulated in disease
meta$condition <- factor(tolower(meta$condition), levels = c("control", "diseased"))

dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData   = meta,
                              design    = ~ condition)
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "diseased", "control"))
write.csv(as.data.frame(res), "DESeq2_results.csv")

# Progressive fold-change thresholds (padj < 0.05)
for (lfc in 1:3) {
  sig <- subset(res, padj < 0.05 & abs(log2FoldChange) >= lfc)
  out <- data.frame(gene = rownames(sig), log2FoldChange = sig$log2FoldChange)
  write.csv(out, sprintf("DESeq2_significant_geneslog%d.csv", lfc), row.names = FALSE)
  cat(sprintf("|log2FC| >= %d & padj < 0.05 : %d genes\n", lfc, nrow(sig)))
}
# Expected: 989 (>=1), 294 (>=2), 86 (>=3). The 86 set == DEGs_86_final.csv.
