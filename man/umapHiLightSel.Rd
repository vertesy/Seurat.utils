% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.Visualization.R
\name{umapHiLightSel}
\alias{umapHiLightSel}
\title{Highlight Selected Clusters on UMAP}
\usage{
umapHiLightSel(
  obj = combined.obj,
  COI = c("0", "2", "4", "5", "11"),
  res.cl = "integrated_snn_res.0.3"
)
}
\arguments{
\item{obj}{Seurat object to be visualized; Default: \code{combined.obj}.}

\item{COI}{Vector of cluster IDs to highlight on the UMAP plot;
Default: \code{c("0", "2", "4", "5", "11")}.}

\item{res.cl}{Name of the metadata column containing cluster IDs;
Default: 'integrated_snn_res.0.3'.}
}
\value{
Saves a UMAP plot highlighting specified clusters to the current working directory.
The function itself does not return an object within R.
}
\description{
Generates a UMAP plot from a Seurat object with specified clusters highlighted.
It saves the resulting UMAP plot directly to the current working directory.
}
\examples{
\dontrun{
if (interactive()) {
  umapHiLightSel(obj = combined.obj, COI = c("0", "1"), res.cl = "resolution_0.8")
}
}

}
\seealso{
\code{\link[Seurat]{DimPlot}}
}
