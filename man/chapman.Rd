\name{chapman}
\alias{chapman}
\title{Chapman Robson catch curve survival estimator}
\usage{
chapman(ages)
}
 
\arguments{
\item{ages}{vector of individual ages, coded such that the first fully-recruited age class is coded as zero.}
}
 
\description{Computes the Chapman-Robson survival estimator and it's standard error from age structure data}
 
\value{list with elemets S_hat and SE(S_hat)}
 
 
\author{Michael Scroggie}
 
\examples{
data(geocrinia)
chapman(geocrinia$age)

randata<-rgeom(1000, 1-0.3)
chapman(randata)
}
 
