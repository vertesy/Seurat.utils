% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.R
\name{downsampleListSeuObjsPercent}
\alias{downsampleListSeuObjsPercent}
\title{Downsample a List of Seurat Objects to a Fraction}
\usage{
downsampleListSeuObjsPercent(
  ls.obj = ls.Seurat,
  fraction = 0.1,
  seed = 1989,
  save_object = FALSE
)
}
\arguments{
\item{ls.obj}{List of Seurat objects to be downsampled; Default: \code{ls.Seurat}.}

\item{fraction}{Fraction of cells to retain in each Seurat object; Default: 0.1.}

\item{save_object}{Logical indicating whether to save the downsampled Seurat objects using
\code{isaveRDS} or return them; Default: FALSE.}
}
\description{
Downsampling a list of Seurat objects to a specified fraction of their original size.
This is useful for reducing dataset size for quicker processing or testing workflows.
}
\examples{
\dontrun{
if (interactive()) {
  downsampled_objs <- downsampleListSeuObjsPercent(ls.obj = yourListOfSeuratObjects, fraction = 0.1)
}
}

}
