% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getDriverRelatedPathway.R
\name{getDriverRelatedPathway}
\alias{getDriverRelatedPathway}
\title{Get driver related pathway/phenotype}
\usage{
getDriverRelatedPathway(
  drivingforcelist,
  enrichreactomeRes,
  drivergenes,
  quary
)
}
\arguments{
\item{drivingforcelist}{A list consists with driving force matrix of each patient.}

\item{enrichreactomeRes}{Abnormal pathways for each patient.}

\item{drivergenes}{A list with individual driver list and ranked population driver list.}

\item{quary}{indicating the query driver gene}
}
\value{
A dataframe contains the information of driver-associated pathways
}
\description{
Calculating driver-associated abnormal pathways for a giving driver gene
}
\examples{
\dontrun{
tp53relatedpathway <- getDriverRelatedPathway(drivingforcelist = genescore_individual,
                                              enrichreactomeRes = enrichreactomeRes,
                                              drivergenes = drivergenes,
                                              quary = 'TP53')
}
}
\author{
Yan Li
}
