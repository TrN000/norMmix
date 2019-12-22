\name{dfnMm}
\alias{dfnMm}
\title{Number of Free Parameters of Multivariate Normal Mixture Models}
\usage{
dfnMm(k, p, model = c("EII","VII","EEI","VEI","EVI",
	               "VVI","EEE","VEE","EVV","VVV"))
}
\arguments{
  \item{k}{number of mixture components}
  \item{p}{dimension of data space, i.e., number of variables (aka \dQuote{features}).}
  \item{model}{a \code{\link{character}} string.  One of the 10 models
    above, see also \sQuote{Description}.}
}
\description{
  \code{npar()} returns an integer (vector, if \code{p} or \code{k} is)
  with the number of free parameters of the corresponding model, which is
  also the \code{\link{length}(.)} of the parameter vector in our
  parametrization, see \code{\link{nMm2par}()}.
}
\examples{
(m <- eval(formals(dfnMm)$model)) # list of 10 models w/ differing Sigma
# A nice table for a given 'p'  and all models, all k in 1:8
sapply(m, dfnMm, k=setNames(,1:8), p = 20)
}
