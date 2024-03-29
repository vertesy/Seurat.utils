% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.Visualization.R
\name{ww.check.quantile.cutoff.and.clip.outliers}
\alias{ww.check.quantile.cutoff.and.clip.outliers}
\title{Check Quantile Cutoff and Clip Outliers}
\usage{
ww.check.quantile.cutoff.and.clip.outliers(
  expr.vec = plotting.data[, gene],
  quantileCutoffX = quantileCutoff,
  min.cells.expressing = 10
)
}
\arguments{
\item{expr.vec}{A numeric vector representing gene expression data.}

\item{quantileCutoffX}{The quantile cutoff for clipping outliers.}

\item{min.cells.expressing}{The minimum number of cells that should remain expressing after clipping.}
}
\value{
The expression vector with outliers clipped, ensuring the minimum number of cells expressing.
}
\description{
Checks a specified quantile cutoff and clips outliers from an expression vector,
ensuring that a minimum number of cells expressing a gene remain.
}
\examples{
\dontrun{
expr.vec <- c(...)
quantileCutoff <- 0.99
min.cells.expressing <- 10
ww.check.quantile.cutoff.and.clip.outliers(expr.vec, quantileCutoff, min.cells.expressing)
}
}
