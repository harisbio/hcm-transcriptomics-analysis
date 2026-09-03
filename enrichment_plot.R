# =============================================================================
# GO enrichment bar plot
# Plots -log10(P-value) for enriched GO terms, coloured by GO domain.
# Works with a tabular GO enrichment export (ShinyGO or Enrichr) containing
# at least the columns: Category, Term, PValue.
# Input: go_enrichment.csv
# =============================================================================

library(ggplot2)
library(readr)
library(dplyr)
library(forcats)

go <- read_csv("go_enrichment.csv")
colnames(go) <- make.names(colnames(go))

go <- go %>%
  mutate(NegLog10P = -log10(PValue)) %>%
  arrange(NegLog10P) %>%
  mutate(Term = factor(Term, levels = Term))

# Domain colours (tolerant of common category labels)
domain_colors <- c(
  "Biological Process" = "#3ef0bd", "GOTERM_BP_DIRECT" = "#3ef0bd", "BP" = "#3ef0bd",
  "Cellular Component" = "#ebb48a", "GOTERM_CC_DIRECT" = "#ebb48a", "CC" = "#ebb48a",
  "Molecular Function" = "#d1f598", "GOTERM_MF_DIRECT" = "#d1f598", "MF" = "#d1f598"
)

ggplot(go, aes(x = NegLog10P, y = fct_rev(Term), fill = Category)) +
  geom_col() +
  scale_fill_manual(values = domain_colors, na.value = "grey70") +
  theme_minimal() +
  labs(title = "GO Functional Enrichment",
       x = "-log10(P-value)", y = NULL, fill = "GO domain")
