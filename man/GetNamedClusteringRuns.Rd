% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.R
\name{GetNamedClusteringRuns}
\alias{GetNamedClusteringRuns}
\title{GetNamedClusteringRuns}
\usage{
GetNamedClusteringRuns(
  obj = combined.obj,
  res = list(FALSE, 0.5)[[1]],
  topgene = FALSE,
  pat = c("^cl.names.Known.*[0,1]\\\\.[0-9]$", "Name|name")[2]
)
}
\arguments{
\item{obj}{Seurat object, Default: combined.obj}

\item{res}{Clustering resoluton to use, Default: c(FALSE, 0.5)\link{1}}

\item{topgene}{Match clustering named after top expressed gene (see vertesy/Seurat.pipeline/~Diff gene expr.), Default: FALSE}

\item{pat}{Pettern to match, Default: '^cl.names.Known.*\link{0,1}\.\link{0-9}$'}
}
\description{
Get Clustering Runs: metadata column names
}
\examples{
\dontrun{
if (interactive()) {
  GetNamedClusteringRuns()
}
}
}
