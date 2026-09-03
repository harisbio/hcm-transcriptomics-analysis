# =============================================================================
# Install R / Bioconductor packages used in this repository.
# Analysis environment (as reported in the manuscript): R v4.5.3, DESeq2 v1.40.2.
# Run once:  source("install_packages.R")
# =============================================================================

# CRAN packages
cran_pkgs <- c("ggplot2", "ggrepel", "dplyr", "tidyr", "readr", "forcats",
               "pheatmap", "RColorBrewer", "circlize", "factoextra", "tibble")

new_cran <- cran_pkgs[!(cran_pkgs %in% installed.packages()[, "Package"])]
if (length(new_cran)) install.packages(new_cran, repos = "https://cloud.r-project.org")

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioc_pkgs <- c("DESeq2", "ComplexHeatmap")
new_bioc <- bioc_pkgs[!(bioc_pkgs %in% installed.packages()[, "Package"])]
if (length(new_bioc)) BiocManager::install(new_bioc, update = FALSE, ask = FALSE)

cat("Dependencies installed.\n")
# For a reproducible record of exact versions, run:  sessionInfo()
