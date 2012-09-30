\name{geomodel}
\alias{geomodel}
\title{Linear models for catch curve data, using the geometric probability distribution}
\description{Fits linear models to catch curve data with covariates using maximum likelihood}
\usage{
geomodel(formula, data)
}
 
\arguments{
\item{formula}{formula relating the age data to the covariates}
\item{data}{a dataframe containing the variables in the formula}
}

\value{A list with elements "formula", "coefficients", "vcov", "logLik", "df", "AICc"}
 
\author{Michael Scroggie}

\references{Scroggie, M.P. (2012) Survival of adult smooth froglets (\emph{Geocrinia laevis} complex, Anura, Myobatrachidae) in and around a hybrid zone. \emph{Herpetological Conservation and Biology} 152: 136-144.}
 
\examples{
data(geocrinia)
geomodel(age~sex, data=geocrinia)
geomodel(age~1, data=geocrinia)
}
