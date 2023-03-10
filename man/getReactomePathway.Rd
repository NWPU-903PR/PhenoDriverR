% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getReactomePathway.R
\name{getReactomePathway}
\alias{getReactomePathway}
\title{Download Reactome Pathway file}
\usage{
getReactomePathway(
  directory,
  mingenenum = 10,
  maxgenenum = 200,
  cleanDownload = FALSE
)
}
\arguments{
\item{directory}{Directory/Folder where the data was downloaded.}

\item{mingenenum}{only pathways consist with larger than 'mingenenum' genes are considered.}

\item{maxgenenum}{only pathways consist with lower than 'maxgenenum' genes are considered.}

\item{cleanDownload}{Should downloaded data be deleted?}
}
\value{
A list of Reactome pathways information
}
\description{
Download Reactome Pathway file from Reactome web site
}
\examples{
\dontrun{
reactome <- getReactomePathway(directory = './',
                               cleanDownload = T)
}
}
\author{
Yan Li
}
