% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.Visualization.R
\name{PlotTopGenesPerCluster}
\alias{PlotTopGenesPerCluster}
\title{Plot Top N Differentially Expressed Genes Per Cluster}
\usage{
PlotTopGenesPerCluster(
  obj = combined.obj,
  cl_res = GetClusteringRuns()[1],
  nrGenes = p$n.markers,
  order.by = c("combined.score", "avg_log2FC", "p_val_adj")[1],
  df_markers = obj@misc$df.markers[[paste0("res.", cl_res)]],
  ...
)
}
\arguments{
\item{obj}{Seurat object containing single-cell RNA-seq data and clustering information;
Default: \code{combined.obj}.}

\item{cl_res}{Cluster resolution used to identify distinct clusters for analysis; Default: \code{res}.}

\item{nrGenes}{Number of top DE genes to display for each cluster;
Default: GetClusteringRuns()\link{1}.}

\item{order.by}{Criteria for ranking DE genes within clusters; Default: \code{"combined.score"}.}

\item{df_markers}{Data frame or list of DE genes across clusters. If not provided,
attempts to retrieve from \code{obj@misc$df.markers[[paste0("res.", cl_res)]]};
Default: calculated based on \code{cl_res}.}
}
\description{
Visualizes the top N differentially expressed (DE) genes for each cluster within a
specified clustering resolution of a Seurat object, facilitating the exploration of gene
expression patterns across clusters.
}
\examples{
\dontrun{
if (interactive()) {
  PlotTopGenesPerCluster(obj = combined.obj, cl_res = 0.5, nrGenes = 10)
}
}

}
