% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.R
\name{parallel.computing.by.future}
\alias{parallel.computing.by.future}
\title{parallel.computing.by.future}
\usage{
parallel.computing.by.future(cores = 4, maxMemSize = 4000 * 1024^2)
}
\arguments{
\item{cores}{Number of cores}

\item{maxMemSize}{memory limit}
}
\description{
Run gc(), load multi-session computing and extend memory limits.
}
