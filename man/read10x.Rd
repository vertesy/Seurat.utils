% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.R
\name{read10x}
\alias{read10x}
\title{Load 10X Genomics Data as Seurat Object}
\usage{
read10x(dir)
}
\arguments{
\item{dir}{A character string specifying the path to the directory containing the 10X dataset files.
This directory should contain \code{matrix.mtx.gz}, \code{features.tsv.gz}, and \code{barcodes.tsv.gz} files.}
}
\value{
A Seurat object containing the single-cell RNA-seq data extracted from the provided 10X
Genomics dataset.
}
\description{
Reads 10X Genomics dataset files (gzipped) including matrix, features, and barcodes,
to a single expression matrix. This function handles the unzipping of these files, reads the data,
and re-compresses the files back to their original gzipped format.
}
\details{
This function facilitates the loading of 10X Genomics datasets into R for analysis with
the Seurat package. It specifically caters to gzipped versions of the \code{matrix.mtx}, \code{features.tsv},
and \code{barcodes.tsv} files, automating their decompression, reading, and subsequent recompression.
The function relies on Seurat's \code{Read10X} function for data reading and object construction.
}
\note{
Ensure that the specified directory contains the required gzipped files.
If the \code{features.tsv.gz} file is named differently (e.g., \code{genes.tsv.gz}), please rename it
accordingly before running this function.
}
\examples{
\dontrun{
if (interactive()) {
  # Replace `path_to_10x_data` with the path to your 10X data directory
  seuratObject <- read10x(dir = "path_to_10x_data")
  # `seuratObject` is now a Seurat object containing the loaded 10X data
}
}

}
\seealso{
\code{\link[Seurat]{Read10X}} for the underlying function used to read the 10X data.
}
