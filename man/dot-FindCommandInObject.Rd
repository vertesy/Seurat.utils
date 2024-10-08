% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.R
\name{.FindCommandInObject}
\alias{.FindCommandInObject}
\title{Find Command in Seurat Object by Partial Match}
\usage{
.FindCommandInObject(obj, pattern, perl = TRUE)
}
\arguments{
\item{obj}{A Seurat object. \strong{Default:} None.}

\item{pattern}{A character string representing the pattern to match command names. \strong{Default:} None.}
}
\value{
If exactly one match is found, the function returns the content of the first match. If
multiple matches are found, it returns \code{NULL} after displaying the number of matches and their names.
}
\description{
This function searches for commands in a list within a Seurat object using a partial match
(e.g., pattern matching) on the command names. It returns the content of the first match if only
one match is found. If multiple matches are found, it outputs the number of hits and their names.
}
\examples{
# Assuming 'combined.obj' is your Seurat object
result <- FindCommandInObject(combined.obj, "^FindVariable")

}
