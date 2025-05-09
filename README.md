# DESeq2_volcano
This repository contains an R script to perform differential gene expression analysis using the DESeq2 package. The pipeline is tailored for RNA-seq count data and includes result reporting, diagnostic plots, annotation integration, and Excel export.

## Features

- Automated DESeq2 pipeline
- Normalization and dispersion estimation
- Generation of key plots: heatmap, PCA, volcano, and correlogram
- Export of full DESeq2 results and filtered (FDR) results
- Integration of gene annotations (optional)
- Export of annotated results as `.tsv` and `.xlsx`
- Color-formatted Excel tables with frozen headers

## Requirements

This script uses the following R packages:

- `DESeq2`
- `ggplot2`, `ggrepel`, `RColorBrewer`, `gplots`, `corrplot`, `ggfortify`, `pheatmap`, `tidyverse`, `tibble`
- `DEGreport`
- `openxlsx`

Install them with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "DEGreport", "pheatmap"))

install.packages(c("ggplot2", "ggrepel", "RColorBrewer", "gplots", "corrplot", "ggfortify", "tibble", "openxlsx", "tidyverse"))

