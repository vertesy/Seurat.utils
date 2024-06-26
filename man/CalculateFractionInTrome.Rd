% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.R
\name{CalculateFractionInTrome}
\alias{CalculateFractionInTrome}
\title{Calculate Fraction of Genes in Transcriptome}
\usage{
CalculateFractionInTrome(
  genesCalc.Cor.Seuratet = c("MALAT1"),
  obj = combined.obj,
  data.slot = c("counts", "data")[2]
)
}
\arguments{
\item{obj}{A Seurat object containing gene expression data; Default: \code{combined.obj}.
The function extracts gene expression data from this object to calculate fractions.}

\item{data.slot}{The data slot from which to extract expression data. This can be \code{"counts"}
for raw counts or \code{"data"} for normalized data; Default: second element (\code{"data"}).}

\item{geneset}{A character vector of gene symbols for which the fraction in the transcriptome will be calculated.
Default: \code{c("MALAT1")}. The function will check for the existence of these genes in the Seurat object.}
}
\value{
A numeric vector where each element represents the fraction of the specified geneset's expression
relative to the total transcriptome of a cell, expressed as a percentage. The names of the vector correspond to cell IDs.
}
\description{
Calculates the fraction of specified genes within the entire transcriptome of
each cell in a Seurat object.
This function is useful for assessing the relative abundance of a set of genes across cells,
such as identifying cells with high expression of marker genes.
}
\note{
This function calls \code{check.genes} to verify the existence of the specified genes within the Seurat object.
If genes are not found, it will return a warning.
}
\examples{
\dontrun{
if (interactive()) {
  # Assuming `combined.obj` is your Seurat object
  fractionInTranscriptome <- CalculateFractionInTranscriptome(geneset = c("MALAT1", "GAPDH"), obj = combined.obj)
  # This will return the fraction of MALAT1 and GAPDH in the transcriptome of each cell
}
}

}
\seealso{
\code{\link[Seurat]{GetAssayData}} for retrieving expression data from a Seurat object.
}
