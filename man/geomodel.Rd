\name{geomodel}
\title{Linear models for catch curve data, using the geometric probability distribution}
\description{fits linear models to catch curve data with covariates using maximum likelihood}
\usage{
geomodel(formula, data)
}
 
\arguments{
\item{formula}{formula relating the age data to associated covariates}
}
 
 
\value{A list with elements "formula", "coefficients", "vcov", "logLik", "df", "AICc"}
 
 
\author{Michael Scroggie}
 
\examples{
data(geocrinia)
geomodel(age~sex)
geomodel(age~1)
}
