% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.R
\name{recall.all.genes}
\alias{recall.all.genes}
\title{Recall all.genes global variable from a Seurat object}
\usage{
recall.all.genes(obj = combined.obj, overwrite = FALSE)
}
\arguments{
\item{obj}{Seurat object, Default: combined.obj}
}
\description{
all.genes set by calc.q99.Expression.and.set.all.genes() #
}
\examples{
\dontrun{
if (interactive()) {
  recall.all.genes()
  all.genes
}
}
}
