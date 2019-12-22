% Code in ../R/llnorMmix.R
\name{sllnorMmix}
\alias{sllnorMmix}
\title{Simple wrapper for Log-Likelihood Function or Multivariate Normal Mixture}
\usage{
sllnorMmix(x, obj, trafo=c("clr1", "logit"))
}
\arguments{
  \item{x}{data \code{\link{matrix}}.}
  \item{obj}{an \R object of class \code{"\link{norMmix}"}.}
  \item{trafo}{a \code{\link{character}} string specifying the transform to
    use for the weights, see  \code{\link{llnorMmix}}.}
}
\description{
  \code{sllnorMmix()} returns a number, the log-likelihood of the data
  \code{x}, given a normal mixture \code{obj}.
}
\details{
  Calculates log-likelihood of a dataset, x, given a normal mixture model;
  just a simplified wrapper for \code{\link{llnorMmix}}.
  Removes functionality in favor of ease of use.
}
\examples{
set.seed(2019)
x <- rnorMmix(400, MW27)
sllnorMmix(x, MW27) # -1986.315
}

