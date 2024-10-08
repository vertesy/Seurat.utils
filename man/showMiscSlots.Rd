% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.R
\name{showMiscSlots}
\alias{showMiscSlots}
\title{Display Slots in the @misc of an Seurat Object}
\usage{
showMiscSlots(obj = combined.obj, max.level = 1, subslot = NULL, ...)
}
\arguments{
\item{obj}{An object whose \verb{@misc} slot needs to be examined. Default: combined.obj}

\item{max.level}{Max depth to dive into sub-elements.}

\item{subslot}{A subslot within \verb{@misc}.}

\item{...}{...}
}
\description{
See \code{showToolsSlots} for details. Prints the names of slots in the \verb{@misc} of a given object.
It specifically targets list elements, skipping over data frames and other non-list objects.
}
\examples{
showToolsSlots(obj)

}
