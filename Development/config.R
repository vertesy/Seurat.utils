# Configuration for the Package
# file.edit("~/GitHub/Packages/Seurat.utils/Development/config.R")

DESCRIPTION <- list(
  package.name = "Seurat.utils",
  version = "2.9.1",
  title = "Seurat.utils - utility functions for Seurat v5",
    description = "Seurat.utils is a collection of utility functions for Seurat single cell analysis.
      Functions allow 3D plotting, visualisation of statistics & QC,
      the automation / multiplexing of plotting, interaction with the Seurat object, etc.
      Some functionalities require functions from CodeAndRoll2 and MarkdownReports libraries.",
  depends = "tidyverse, Seurat, Stringendo, CodeAndRoll2, ggExpress, magrittr",
  imports = "cowplot, dplyr, ggcorrplot, ggpubr, ggrepel, HGNChelper, htmlwidgets,
    Matrix, matrixStats, princurve, pheatmap,
    R.utils, readr, reshape2, scales, SoupX, sparseMatrixStats, stringr, tibble, tictoc,
    ReadWriter, MarkdownHelpers, MarkdownReports,
    plotly, qs, foreach, harmony, EnhancedVolcano, rstudioapi,
    RColorBrewer, SeuratObject, checkmate, fs, future,
    ggplot2, gplots", # , grDevices
  suggests = "clusterProfiler, enrichplot, DatabaseLinke.R, vroom",
  author.given = "Abel",
  author.family = "Vertesy",
  author.email = "av@imba.oeaw.ac.at",
  github.user = "vertesy",
  license = "GPL-3 + file LICENSE"
)
