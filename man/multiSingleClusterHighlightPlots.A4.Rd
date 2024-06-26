% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.Visualization.R
\name{multiSingleClusterHighlightPlots.A4}
\alias{multiSingleClusterHighlightPlots.A4}
\title{Generate Cluster Highlight UMAPs compiled into A4 pages}
\usage{
multiSingleClusterHighlightPlots.A4(
  obj = combined.obj,
  ident = GetClusteringRuns(obj)[1],
  foldername = ident,
  plot.reduction = "umap",
  intersectionAssay = DefaultAssay(combined.obj),
  layout = c("tall", "wide", FALSE)[2],
  colors = c("grey", "red"),
  nr.Col = 2,
  nr.Row = 4,
  cex = round(0.1/(nr.Col * nr.Row), digits = 2),
  subdir = TRUE,
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
\item{obj}{A Seurat object combining multiple datasets. Default: \code{combined.obj}.}

\item{ident}{The name of the metadata column in the Seurat object \code{obj} to use for identifying clusters.}

\item{foldername}{Name of the folder to save the plots in. Default: Value of \code{ident}.}

\item{plot.reduction}{The dimensionality reduction technique to use for the plots. Default: \code{"umap"}.}

\item{intersectionAssay}{The assay to use when calculating intersections. Default: \code{"RNA"}.}

\item{layout}{Plot layout, can be \code{"tall"}, \code{"wide"}, or \code{FALSE} for no specific layout. Default: \code{"wide"}.}

\item{colors}{A vector of colors to use for non-highlighted and highlighted clusters. Default: \code{c("grey", "red")}.}

\item{nr.Col}{Number of columns in the plot grid. Default: 2.}

\item{nr.Row}{Number of rows in the plot grid. Default: 4.}

\item{cex}{Size of the text in the plot, calculated based on the number of rows and columns. Default: Calculated value.}

\item{subdir}{Logical flag indicating whether to create a subdirectory for the plots. Default: \code{TRUE}.}

\item{prefix}{Optional prefix for the plot file names. Default: \code{NULL}.}

\item{suffix}{Optional suffix for the plot file names. Default: \code{NULL}.}

\item{background_col}{Background color of the plots. Default: \code{"white"}.}

\item{aspect.ratio}{Aspect ratio of the plots, can be \code{FALSE} for default ratio or a numeric value. Default: 0.6.}

\item{saveGeneList}{Logical flag indicating whether to save the list of genes used in the plots. Default: \code{FALSE}.}

\item{w}{Width of the plots, in inches. Default: \code{8.27}.}

\item{h}{Height of the plots, in inches. Default: \code{11.69}.}

\item{scaling}{Scaling factor for adjusting the size of the plots. Default: 1.}

\item{format}{Format to save the plots in, can be \code{"jpg"}, \code{"pdf"}, or \code{"png"}. Default: \code{"jpg"}.}

\item{...}{Additional arguments passed to lower-level plotting functions.}
}
\value{
Invisible. This function primarily saves plots to files.
}
\description{
This function generates and saves cluster highlight plots for both single and multiple
clusters using UMAP or other dimensionality reduction techniques. It supports saving plots in various
formats and allows customization of plot appearance and layout.
}
\examples{
multiSingleClusterHighlightPlots.A4(ident = "cluster_id", obj = yourSeuratObject)

}
