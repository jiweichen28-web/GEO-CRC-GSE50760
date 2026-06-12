## ── 安装说明 ────────────────────────────────────────────────
## 在 R 里运行以下命令安装全部依赖，只需执行一次

# CRAN 包
install.packages(c(
  "tidyverse",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "RColorBrewer",
  "dplyr",
  "stringr",
  "openxlsx"
))

# Bioconductor 包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "GEOquery",
  "limma",
  "DESeq2",
  "edgeR",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "pathview",
  "fgsea",
  "DOSE"
))
