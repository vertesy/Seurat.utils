% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.R
\name{.getNrScaledFeatures}
\alias{.getNrScaledFeatures}
\title{Get number of scaled features}
\usage{
.getNrScaledFeatures(obj, assay = Seurat::DefaultAssay(obj))
}
\arguments{
\item{obj}{A Seurat object containing scaled data in  \code{obj@assays$RNA@scale.data}.}
}
\value{
Integer representing the number of scaled features
}
\description{
Get number of scaled features
}
