% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.Metadata.R
\name{metaColnameExists}
\alias{metaColnameExists}
\title{Check if a Column Exists in the Metadata of an S4 Object}
\usage{
metaColnameExists(col_name, obj = combined.obj)
}
\arguments{
\item{col_name}{A character string specifying the name of the column.}

\item{obj}{A Seurat object.}
}
\value{
A logical value indicating whether the column exists (TRUE) or not (FALSE).
}
\description{
This function checks whether a given column exists in the meta.data of a Seurat object.
}
