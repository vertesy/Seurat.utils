% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.Metadata.R
\name{transferMetadata}
\alias{transferMetadata}
\title{Transfer Multiple Metadata Columns Between Two Seurat Objects}
\usage{
transferMetadata(
  from,
  to,
  colname_from,
  colname_to = colname_from,
  overwrite = FALSE,
  plotUMAP = TRUE,
  strict = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{from}{The source Seurat object from which metadata will be transferred.}

\item{to}{The destination Seurat object to which metadata will be added.}

\item{colname_from}{Vector of names for the columns in the source object's metadata to transfer.}

\item{colname_to}{Vector of names for the columns in the destination object's metadata.
Defaults to the same names as \code{colname_from}. Must be the same length as \code{colname_from} unless
it is the same as \code{colname_from}.}

\item{overwrite}{Logical, indicating whether to overwrite the column in the destination object
if it already exists. Defaults to FALSE.}

\item{plotUMAP}{Logical, indicating whether to plot UMAPs of the destination object with
the new identity.}

\item{strict}{Logical, indicating whether to fail if the destination object have cells not found in the source object.}

\item{verbose}{Logical, indicating whether to print details about the transfer, including the
number and percentage of matching cells between objects, and unique cells in each object.}

\item{...}{Additional arguments to be passed to \code{transferMetadata}.}
}
\value{
Returns the destination Seurat object (\code{to}) with the new metadata columns added.
}
\description{
Transfers specified metadata columns from one Seurat object to another,
with options for verbose output and overwriting existing columns. Checks for cell overlap and
reports percentages of matching and unique cells.
}
\details{
This function is useful for merging related data from separate Seurat objects,
ensuring that relevant metadata is consistent across datasets. The function checks for
the existence of the specified columns in the source object and whether the columns
can be overwritten in the destination object. It also reports on cell overlap between
the two objects, which can be useful for understanding the relationship between datasets.
}
\examples{
# Assuming `object1` and `object2` are Seurat objects, and you want to transfer
# metadata columns named 'patientID' and 'treatmentGroup' from `object1` to `object2`:
object2 <- transferMetadata(
  from = object1, to = object2,
  colname_from = c("patientID", "treatmentGroup")
)

}
