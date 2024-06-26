% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.R
\name{runDGEA}
\alias{runDGEA}
\title{Run Differential Gene Expression Analysis (DGEA)}
\usage{
runDGEA(
  obj,
  param.list = p,
  ident = NULL,
  res.analyzed.DE = if (is.null(ident)) c(0.1) else ident,
  reorder.clusters = if (is.null(ident)) TRUE else FALSE,
  reorder.dimension = 1,
  directory = OutDirOrig,
  subdirectory = "DGEA_res",
  add.combined.score = TRUE,
  save.obj = TRUE,
  calculate.DGEA = TRUE,
  plot.DGEA = TRUE,
  plot.av.enrichment.hist = TRUE,
  plot.log.top.gene.stats = TRUE,
  auto.cluster.naming = TRUE,
  clean.misc.slot = TRUE,
  clean.meta.data = TRUE,
  ...
)
}
\arguments{
\item{obj}{Seurat object, assumed to be pre-configured with necessary data.}

\item{param.list}{List of parameters for DE analysis. Default: p.}

\item{ident}{Use this to specify a non-standard cluster identity, such as named clusters.
\code{runDGEA} will use this ident explicitly for the DE analysis. Default: NULL.}

\item{res.analyzed.DE}{Vector of numeric values specifying the resolutions to analyze.
Default: c(0.1).}

\item{reorder.clusters}{Logical indicating whether to reorder clusters based on dimension.
Default: TRUE.}

\item{reorder.dimension}{Integer specifying the dimension for reordering (1 for x, -1 for y).
Default: 1.}

\item{directory}{Character string specifying the base directory for saving results.
Default: OutDirOrig.}

\item{subdirectory}{Character string specifying the subdirectory for saving outputs within
the base directory. Default: "DGEA".}

\item{add.combined.score}{Logical indicating whether to add a combined score to the markers.
Default: TRUE.}

\item{save.obj}{Logical indicating whether to save the modified Seurat object.
Default: TRUE.}

\item{calculate.DGEA}{Logical determining if the DE analysis should be calculated.
Default: TRUE.}

\item{plot.DGEA}{Logical determining if results should be plotted.
Default: TRUE.}

\item{plot.av.enrichment.hist}{Logical indicating whether to plot the average enrichment histogram.
Default: TRUE.}

\item{plot.log.top.gene.stats}{Logical indicating whether to plot the log top gene statistics.}

\item{auto.cluster.naming}{Logical indicating automatic labeling of clusters.
Default: TRUE.}

\item{clean.misc.slot}{Logical indicating whether to clean the misc slots of previous
clustering results. Default: TRUE.}

\item{clean.meta.data}{Logical indicating whether to clean the metadata slots of
previous clustering results. Default: TRUE.}
}
\value{
Modified Seurat object and markers list.
}
\description{
Runs a differential gene expression analysis based on specified parameters,
reorders clusters if needed, and optionally saves results. Supports output and plotting configurations.
}
\examples{
runDGEA(obj = mySeuratObject, param.list = myListParams, directory = "Results/MyAnalysis")

}
