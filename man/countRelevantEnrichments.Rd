% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.Visualization.R
\name{countRelevantEnrichments}
\alias{countRelevantEnrichments}
\title{Count Relevant Enrichments}
\usage{
countRelevantEnrichments(
  df,
  pval_col = "p_val_adj",
  logfc_col = "avg_log2FC",
  pval_cutoff = 0.01,
  logfc_cutoff = 1
)
}
\arguments{
\item{df}{Data frame containing the gene expression data.}

\item{pval_col}{Character. Name of the column containing the adjusted p-values. Default: "p_val_adj".}

\item{logfc_col}{Character. Name of the column containing the log2 fold change values. Default: "avg_log2FC".}

\item{pval_cutoff}{Numeric. The maximum adjusted p-value to consider a gene relevant. Default: 1e-2.}

\item{logfc_cutoff}{Numeric. The minimum log2 fold change to consider a gene relevant. Default: 1.}
}
\value{
A list with the counts of enriched and depleted genes.
}
\description{
This function counts the number of relevantly expressed genes from a differential
gene expression table. It considers genes to be relevant if they fall under a maximum p-value
cutoff and are above a minimum log2 fold change cutoff. The function reports the number of
enriched and depleted genes.
}
\examples{
df <- data.frame(
  p_val_adj = c(0.001, 0.02, 0.03, 0.0001),
  avg_log2FC = c(1.5, -2, 0.5, 2)
)
countRelevantEnrichments(df)
}
