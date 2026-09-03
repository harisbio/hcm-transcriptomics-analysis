# =============================================================================
# KEGG pathway - gene chord diagram
# Links the top enriched KEGG pathways (Enrichr) to their member DEGs,
# coloured by log2 fold change. Panel label "C" corresponds to the figure.
# Inputs:
#   DEGs_86_final.csv       (columns 'gene', 'log2FoldChange')
#   pathways_enrichr.csv    (Enrichr KEGG export; columns incl. Term, Genes, Adjusted.P.value)
# =============================================================================
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(grid)

deg <- read.csv("DEGs_86_final.csv")
enrichr <- read.csv("pathways_enrichr.csv")

top_pathways <- enrichr %>%
  arrange(Adjusted.P.value) %>%
  head(6) %>%
  mutate(Term = gsub(" - COVID-19.*", "", Term))

plot_df <- top_pathways %>%
  select(Term, Genes) %>%
  separate_rows(Genes, sep = ";") %>%
  inner_join(deg, by = c("Genes" = "gene")) %>%
  select(gene = Genes, pathway = Term, log2FC = log2FoldChange) %>%
  arrange(desc(log2FC))

genes <- unique(plot_df$gene)
pathways <- unique(plot_df$pathway)

col_fun_fc <- colorRamp2(c(-2.4, 0, 2.4), c("blue", "#EEEEEE", "red"))
pathway_palette <- c("#66c2a5", "#ffd92f", "#bebada",
                     "#fb8072", "#80b1d3", "#fdb462", "#b3de69")
pathway_cols <- setNames(pathway_palette[1:length(pathways)], pathways)

add_alpha <- function(col, alpha) {
  tmp <- col2rgb(col)
  rgb(tmp[1], tmp[2], tmp[3], alpha = alpha * 255, maxColorValue = 255)
}

circos.clear()

pathway_w <- plot_df %>%
  group_by(pathway) %>%
  summarize(w = n())

gene_w <- rep(1.5, length(genes))

sector_widths <- c(setNames(gene_w, genes),
                   setNames(pathway_w$w, pathway_w$pathway))

circos.par(
  gap.after = c(rep(1, length(genes) - 1), 15,
                rep(4, length(pathways) - 1), 15),
  start.degree = 120,
  cell.padding = c(0, 0, 0, 0),
  points.overflow.warning = FALSE
)

circos.initialize(names(sector_widths),
                  xlim = cbind(0, sector_widths))

circos.track(ylim = c(0, 1), track.height = 0.25,
             bg.border = NA,
             panel.fun = function(x, y) {

  sector <- get.cell.meta.data("sector.index")
  center <- get.cell.meta.data("xcenter")

  if (sector %in% genes) {
    circos.text(center, 0, sector,
                facing = "clockwise",
                niceFacing = TRUE,
                adj = c(0, 0.5),
                cex = 0.55,
                font = 2)
  } else {
    circos.text(center, 0.3, sector,
                facing = "bending.outside",
                niceFacing = TRUE,
                adj = c(0.5, 0),
                cex = 0.8,
                font = 2)
  }
})

circos.track(ylim = c(0, 1), track.height = 0.03,
             bg.border = NA,
             panel.fun = function(x, y) {

  sector <- get.cell.meta.data("sector.index")

  if (sector %in% genes) {
    val <- plot_df$log2FC[match(sector, plot_df$gene)]
    circos.rect(0, 0, sector_widths[sector], 1,
                col = col_fun_fc(val),
                border = "white", lwd = 0.3)
  } else {
    circos.rect(0, 0, sector_widths[sector], 1,
                col = pathway_cols[sector],
                border = NA)
  }
})

pathway_offsets <- setNames(rep(0, length(pathways)), pathways)

for (i in 1:nrow(plot_df)) {

  g <- plot_df$gene[i]
  p <- plot_df$pathway[i]

  start_p <- pathway_offsets[p]
  end_p <- start_p + 1
  pathway_offsets[p] <- end_p

  circos.link(g, c(0, sector_widths[g]),
              p, c(start_p, end_p),
              col = add_alpha(pathway_cols[p], 0.35),
              border = NA,
              rou = 0.72)
}

lgd_fc <- Legend(
  title = "log2FC",
  col_fun = col_fun_fc,
  at = c(-2.4, -1.2, 0, 1.2, 2.4),
  labels = c("-2.4", "-1.2", "0", "1.2", "2.4"),
  grid_width = unit(4, "mm"),
  grid_height = unit(4, "mm"),
  title_gp = gpar(fontsize = 10, fontface = "bold")
)

draw(lgd_fc,
     x = unit(0.94, "npc"),
     y = unit(0.65, "npc"),
     just = c("right", "top"))

grid.text("C",
          x = unit(0.06, "npc"),
          y = unit(0.94, "npc"),
          gp = gpar(fontsize = 28, fontface = "bold"))
