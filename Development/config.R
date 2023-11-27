# Configuration for the Package
DESCRIPTION <- list(
  package.name = "Seurat.utils",
  version = "2.3.0",
  title = "Seurat.utils - utility functions for Seurat",
  description = "Seurat.utils Is a collection of utility functions for Seurat v3.
    Functions allow the automation / multiplexing of plotting, 3D plotting, visualisation of statistics &
    QC, interaction with the Seurat object, etc. Some functionalities require functions from CodeAndRoll and MarkdownReports libraries.",
  depends = "ggplot2, Seurat, Stringendo, CodeAndRoll2, ggExpress, magrittr",
  imports = "cowplot, Seurat, ReadWriter, dplyr, ggcorrplot, ggpubr, ggrepel, grDevices, HGNChelper,
                          htmlwidgets, MarkdownHelpers, MarkdownReports, Matrix, matrixStats, princurve, pheatmap,
                          R.utils, readr, reshape2, scales, SoupX, sparseMatrixStats, stringr, tibble, tictoc,
                          vroom, job, qs, foreach, tidyverse",
  suggests = "SoupX, princurve, EnhancedVolcano",

  author.given = "Abel",
  author.family = "Vertesy",
  author.email = "av@imba.oeaw.ac.at",
  github.user = "vertesy",
  license = "GPL-3 + file LICENSE"
)
