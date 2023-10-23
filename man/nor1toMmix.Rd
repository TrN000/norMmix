\name{nor1toMmix}
\alias{nor1toMmix}
\title{Cast nor1mix object as norMmix.}
\usage{
nor1toMmix(object)
}
\description{
  Cast nor1mix object as norMmix.
}
\arguments{
  \item{object}{A nor1mix mixture model to be coerced to norMmix.}
}
\value{
A `norMmix` object if the appropriate S3method has been implemented.
}
\details{
  This package was designed to extend the `nor1mix` package to the case of
  multivariate mixture models. Therefore we include a utility function to cast
  1-dimensional mixtures as defined in `nor1mix` to `norMmix`.
}
