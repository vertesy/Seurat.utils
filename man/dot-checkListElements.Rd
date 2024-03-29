% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.R
\name{.checkListElements}
\alias{.checkListElements}
\title{Check List Elements}
\usage{
.checkListElements(param_list = NULL, elements = character(0))
}
\arguments{
\item{param_list}{A list containing variables to be checked. Default: \code{NULL}.}

\item{elements}{A character vector of element names in \code{param_list} to check.
Default: \code{character(0)}.}
}
\value{
A message for each element that is defined, and a warning for elements that are not.
}
\description{
Tests if list elements are defined and reports the value or issues a warning.
}
\examples{
param_list <- list(a = 1, b = NULL)
elements <- c("a", "b", "c")
.checkListElements(param_list, elements)
}
