% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.Visualization.R
\name{multiFeaturePlot.A4}
\alias{multiFeaturePlot.A4}
\title{multiFeaturePlot.A4}
\usage{
multiFeaturePlot.A4(
  list.of.genes,
  obj = combined.obj,
  subdir = TRUE,
  foldername = substitute(list.of.genes),
  plot.reduction = "umap",
  intersectionAssay = c("RNA", "integrated")[1],
  layout = c("tall", "wide", FALSE)[2],
  colors = c("grey", "red"),
  nr.Col = 2,
  nr.Row = 4,
  raster = if (ncol(obj) > 1e+05) TRUE else FALSE,
  cex = round(0.1/(nr.Col * nr.Row), digits = 2),
  cex.min = if (raster) TRUE else FALSE,
  gene.min.exp = "q01",
  gene.max.exp = "q99",
  prefix = NULL,
  suffix = NULL,
  background_col = "white",
  aspect.ratio = c(FALSE, 0.6)[2],
  saveGeneList = FALSE,
  w = 8.27,
  h = 11.69,
  scaling = 1,
  format = c("jpg", "pdf", "png")[1],
  ...
)
}
\arguments{
\item{list.of.genes}{List of gene names for which the plots are to be generated. No default.}

\item{obj}{Seurat object, Default: combined.obj}

\item{subdir}{Should plots be saved in a sub-directory? Default: TRUE}

\item{foldername}{Folder name to save the generated plots. Default: The name of the list of genes.}

\item{plot.reduction}{Dimension reduction technique to use for plots. Default: 'umap'}

\item{intersectionAssay}{The assay to intersect with, either 'RNA' or 'integrated'. Default: 'RNA'}

\item{layout}{Layout orientation of the plot. Default: 'wide'}

\item{colors}{Vector of colors to be used in the plot. Default: c("grey", "red")}

\item{nr.Col}{Number of columns in the plot grid. Default: 2}

\item{nr.Row}{Number of rows in the plot grid. Default: 4}

\item{cex}{Point size in the plot. Default: round(0.1/(nr.Col * nr.Row), digits = 2)}

\item{gene.min.exp}{Minimum gene expression level for plotting. Default: 'q01'}

\item{gene.max.exp}{Maximum gene expression level for plotting. Default: 'q99'}

\item{prefix}{Prefix for the plot filenames. Default: NULL}

\item{suffix}{Suffix for the plot filenames. Default: NULL}

\item{background_col}{Background color of the plots. Default: "white"}

\item{aspect.ratio}{Should the aspect ratio be fixed? Default: Yes, at 0.6}

\item{saveGeneList}{Should the list of genes be saved? Default: FALSE}

\item{w}{Width of the plot. Default: 8.27}

\item{h}{Height of the plot. Default: 11.69}

\item{scaling}{Scaling factor for plot size. Default: 1}

\item{format}{Format to save the plot file. Default: 'jpg'}

\item{...}{Pass any other parameter to the internally called functions (most of them should work).}
}
\description{
Save multiple FeaturePlots, as jpeg, on A4 for each gene, which are stored as a list of gene names.
}
\seealso{
\code{\link[tictoc]{tic}}
\code{\link[cowplot]{plot_grid}}
}
