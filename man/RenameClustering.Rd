% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.R
\name{RenameClustering}
\alias{RenameClustering}
\title{RenameClustering}
\usage{
RenameClustering(
  namedVector = ManualNames,
  orig.ident = "RNA_snn_res.0.3",
  suffix.new.ident = "ManualNames",
  new.ident = ppp(orig.ident, suffix.new.ident),
  obj = combined.obj,
  suffix.plot = "",
  plot_umaps = TRUE,
  ...
)
}
\arguments{
\item{namedVector}{named vector, where values = new, names(vec) = old}

\item{orig.ident}{meta.data colname original}

\item{suffix.new.ident}{How to name (suffix) the new identity. Default: "ManualNames"}

\item{new.ident}{meta.data colname new}

\item{obj}{Seurat object}

\item{suffix.plot}{Suffix description (short string) to be added to the umap plots.}

\item{...}{Pass any other parameter to the internally called functions (most of them should work).}
}
\description{
Rename clustering in a Seurat object.
}
