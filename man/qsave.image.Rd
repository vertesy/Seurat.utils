% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.R
\name{qsave.image}
\alias{qsave.image}
\title{qsave.image}
\usage{
qsave.image(..., showMemObject = T, options = c("--force", NULL)[1])
}
\arguments{
\item{...}{Pass any other parameter to the internally called functions (most of them should work).}

\item{showMemObject}{PARAM_DESCRIPTION, Default: T}

\item{options}{PARAM_DESCRIPTION, Default: c("--force", NULL)\link{1}}
}
\description{
Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression. #
}
\examples{
\dontrun{
if(interactive()){
 #EXAMPLE1
 }
}
}
\seealso{
\code{\link[Stringendo]{kollapse}}, \code{\link[function]{iprint}}
}