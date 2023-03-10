% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getDriverGenes.R
\name{getDriverGenes}
\alias{getDriverGenes}
\title{Get driver genes}
\usage{
getDriverGenes(drivingforcelist, cancerData, topgenes_ind = 5)
}
\arguments{
\item{drivingforcelist}{A list consists with driving force matrix of each patient.}

\item{cancerData}{A list consists with expression and mutation data after per-processing.}

\item{topgenes_ind}{Number of top genes that considered as individual driver genes.}
}
\value{
A list with individual driver list and ranked population driver list
}
\description{
Get individual-level and population-level driver genes from driving force matrix.
}
\examples{
\dontrun{
drivergenes <- getDriverGenes(drivingforcelist = genescore_individual,
                              cancerData = cancerData)
}
}
\author{
Yan Li
}
