% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.R
\name{LoadAllSeurats}
\alias{LoadAllSeurats}
\title{LoadAllSeurats}
\usage{
LoadAllSeurats(
  InputDir,
  file.pattern = "^filtered.+Rds$",
  string.remove1 = c(F, "filtered_feature_bc_matrix.", "raw_feature_bc_matrix.")[2],
  string.remove2 = c(F, ".min.cells.10.min.features.200.Rds")[2]
)
}
\arguments{
\item{InputDir}{Input directory}

\item{file.pattern}{PARAM_DESCRIPTION, Default: '^filtered.+Rds$'}

\item{string.remove1}{PARAM_DESCRIPTION, Default: c(F, "filtered_feature_bc_matrix.", "raw_feature_bc_matrix.")\link{2}}

\item{string.remove2}{PARAM_DESCRIPTION, Default: c(F, ".min.cells.10.min.features.200.Rds")\link{2}}
}
\description{
Load all Seurat objects found in a directory. Also works with symbolic links (but not with aliases). #
}
\examples{
\dontrun{
if(interactive()){
 ls.Seurat <- LoadAllSeurats(InputDir)
 }
}
}