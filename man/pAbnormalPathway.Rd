% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pAbormalPathway.R
\name{pAbnormalPathway}
\alias{pAbnormalPathway}
\title{Abnormal pathways discovering}
\usage{
pAbnormalPathway(
  diff,
  reactome,
  parallelworker,
  th.diff.type = "p.value",
  th.diff = 0.01
)
}
\arguments{
\item{diff}{a list consists with P-value (or P-adjust) matrix which indicates whether
a gene is differentially expressed.}

\item{reactome}{Reactome pathway information}

\item{parallelworker}{number of workers that parallel calculating differentally
expression genes of every patients.}

\item{th.diff.type}{indicating whether P-value or P-adjust used to ditect DEGs.}

\item{th.diff}{giving threshold of differentially expression pvalue.}
}
\value{
the results of enrichment analysis for every patient.
}
\description{
Discovering individual patient's clinically abnormal phenotypes (pathways)
}
\examples{
\dontrun{
abnormalPathway <- pAbnormalPathway(diff = diffexpgene,
                                    reactome = reactome,
                                    parallelworker = 4)
}
}
\author{
Yan Li
}
