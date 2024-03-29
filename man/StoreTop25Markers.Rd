% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.R
\name{StoreTop25Markers}
\alias{StoreTop25Markers}
\title{Save Top 25 Markers per Cluster}
\usage{
StoreTop25Markers(obj = combined.obj, df_markers = df.markers, res = 0.5)
}
\arguments{
\item{obj}{Seurat object to update with top 25 markers information; Default: \code{combined.obj}.}

\item{df_markers}{Data frame containing results from differential gene expression analysis
via \code{FindAllMarkers()}, specifying significant markers across clusters; Default: \code{df.markers}.}

\item{res}{Clustering resolution at which the markers were identified; Default: 0.5.}
}
\description{
Stores the top 25 markers for each cluster identified in a Seurat object, based on
the \code{avg_log2FC} from the output table of \code{FindAllMarkers()}. The result is saved under \verb{@misc$df.markers$res...},
rounding insignificant digits to three decimal places.
}
\examples{
\dontrun{
if (interactive()) {
  combined.obj <- StoreTop25Markers(obj = combined.obj, df_markers = df.markers, res = 0.5)
}
}

}
\seealso{
\code{\link[Seurat]{FindAllMarkers}}, \code{\link[dplyr]{top_n}}
}
