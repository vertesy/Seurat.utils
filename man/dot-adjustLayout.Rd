% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.Visualization.R
\name{.adjustLayout}
\alias{.adjustLayout}
\title{Adjust Layout Parameters for multi* plotting fucntions}
\usage{
.adjustLayout(layout, scaling, wA4, hA4, env)
}
\arguments{
\item{layout}{A string specifying the layout type. Can be either "tall" or "wide". Default: NULL.}

\item{scaling}{A numeric scaling factor to adjust the dimensions. Default: 1.}

\item{wA4}{The width of the A4 paper in inches. Default: 8.27.}

\item{hA4}{The height of the A4 paper in inches. Default: 11.69.}

\item{env}{The environment where the layout dimensions and properties should be assigned.
Default: parent.frame().}
}
\value{
Invisible NULL. The function operates by side effects, updating the \code{env} environment.
}
\description{
Adjusts layout dimensions and properties based on the specified layout type.
Updates the provided environment with new dimensions and layout configuration.
}
\examples{
env <- new.env()
.adjustLayout("tall", 1, 8.27, 11.69, env)
print(env$w) # Should print the width based on "tall" layout scaling.

}
