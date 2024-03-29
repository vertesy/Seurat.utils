% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.Visualization.R
\name{SeuratColorVector}
\alias{SeuratColorVector}
\title{Regenerate Color Scheme for Clusters in Seurat Object as a vector}
\usage{
SeuratColorVector(
  ident = NULL,
  obj = combined.obj,
  plot.colors = FALSE,
  simple = FALSE
)
}
\arguments{
\item{ident}{Specific clustering identity to use for color extraction.
If NULL, the active identity in \code{obj} is used. Default: NULL.}

\item{obj}{Seurat object from which to extract cluster colors.
Default: \code{combined.obj}.}

\item{plot.colors}{If TRUE, visually displays the color scheme.
Default: FALSE.}

\item{simple}{If TRUE, returns only the unique set of colors used.
If FALSE, returns a named vector mapping cluster identities to colors.
Default: FALSE.}
}
\description{
Extracts and optionally displays the color scheme assigned to cluster identities
within a Seurat object, facilitating consistent color usage across visualizations. You can
check results in a barplot with \code{MarkdownHelpers::color_check()}.
}
\examples{
\dontrun{
# Display colors for the active identity
SeuratColorVector()
# Retrieve and plot colors for a specified clustering identity
SeuratColorVector(ident = "RNA_snn_res.1", plot.colors = TRUE)
}

}
\seealso{
\code{\link[scales]{hue_pal}}
}
