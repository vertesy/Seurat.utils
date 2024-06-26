% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.Visualization.R
\name{filterGoEnrichment}
\alias{filterGoEnrichment}
\title{Filter GO Enrichment Results}
\usage{
filterGoEnrichment(
  df.enrichments,
  pvalueCutoff = NULL,
  qvalueCutoff = NULL,
  colname = "Description"
)
}
\arguments{
\item{df.enrichments}{An object of class \code{enrichResult} containing the GO enrichment results.}

\item{pvalueCutoff}{Numeric. The p-value cutoff for filtering the results. Default: NULL, meaning
that the default cutoff of the input object is used. It is stored in \code{df.enrichments@pvalueCutoff}.}

\item{qvalueCutoff}{Numeric. The q-value cutoff for filtering the results. Default: NULL,
meaning that the default cutoff of the input object is used. It is stored in \code{df.enrichments@qvalueCutoff}.}

\item{colname}{Character. The name of the column containing the GO-term names, or else.}
}
\value{
A character vector of descriptions of the filtered GO enrichment results.
}
\description{
This function filters GO enrichment results based on adjusted p-value and q-value
cutoffs, and retrieves the descriptions of the filtered results.
}
\examples{
# Assuming GO.Enriched.DL.Ctrl is an object of class `enrichResult` created by clusterprofiler or equivalent
descriptions <- filterGoEnrichment(GO.Enriched.DL.Ctrl)
print(descriptions)

}
