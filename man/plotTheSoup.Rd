% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.R
\name{plotTheSoup}
\alias{plotTheSoup}
\title{plotTheSoup}
\usage{
plotTheSoup(
  CellRanger_outs_Dir = "~/Data/114593/114593",
  SeqRun = str_extract(CellRanger_outs_Dir, "[[:alnum:]_]+(?=/outs/)"),
  ls.Alpha = 1
)
}
\arguments{
\item{CellRanger_outs_Dir}{CellRanger 'outs' (output) directory, Default: '~/Data/114593/114593'}

\item{SeqRun}{Aka SampleName (the folder above 'outs;). Default: str_extract(CellRanger_outs_Dir, "[\link{:alnum:}_]+(?=/outs/)").}
}
\description{
Plot stats about the ambient RNA content in a 10X experiment.
}
\seealso{
\code{\link[Matrix]{colSums}}
\code{\link[tibble]{rownames}}
\code{\link[ggrepel]{geom_label_repel}}
}
