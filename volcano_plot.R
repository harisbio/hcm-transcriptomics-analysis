# =============================================================================
# Volcano plot of differential expression
# Dashed vertical guides at |log2FC| = 1, 2, 3 illustrate the progressive
# fold-change thresholds; the high-confidence set uses |log2FC| >= 3.
# Input:  DESeq2_results.csv  (from deseq2_analysis.R)
# =============================================================================

library(ggplot2)
library(ggrepel)
library(dplyr)

res_df <- read.csv("DESeq2_results.csv", row.names = 1)
res_df <- res_df[!is.na(res_df$padj), ]

res_df$Significance <- ifelse(
  res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 3,
  ifelse(res_df$log2FoldChange > 0, "Upregulated", "Downregulated"),
  "Not Significant"
)

# Label the strongest high-confidence genes (top 20 by adjusted p-value)
top_genes <- res_df %>%
  tibble::rownames_to_column("gene") %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 3) %>%
  arrange(padj) %>%
  slice_head(n = 20)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.8, size = 1) +
  scale_color_manual(values = c("Upregulated"     = "red",
                                "Downregulated"   = "blue",
                                "Not Significant" = "grey60")) +
  geom_vline(xintercept = c(-3, -2, -1, 1, 2, 3), linetype = "dashed",
             colour = "grey70", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "grey70", linewidth = 0.3) +
  geom_text_repel(data = top_genes,
                  aes(x = log2FoldChange, y = -log10(padj), label = gene),
                  inherit.aes = FALSE, size = 2, max.overlaps = 30) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Expression",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value")
