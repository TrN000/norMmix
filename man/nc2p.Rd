% documentation in R/param.R
\name{nc2p}
\title{Wrapper function for \code{\link{nMm2par}}}
\alias{nc2p}
\description{
    \code{nc2p} returns same as \code{\link{nMm2par}}, using \code{object$model}
    and \code{trafo="clr1"} as defaults.
}
\usage{
    nc2p(object)
}
\arguments{
  \item{object}{an \code{"\link{norMmix}"} object.}
}
\value{
    the same as \code{\link{nMm2par}}. Real valued vector of length 1.

    \code{nMm2par(object, trafo="clr1", model=object$model)}
}
\examples{
    str(MW213)
    nc2p(MW213)
    #  [1]  0.0000000  0.0000000  0.0000000 30.0000000 30.0000000  0.3465736
    #  [7]  0.5493061  0.3465736 -0.5493061  3.0000000  2.0000000
    nMm2par(MW213, trafo="clr1", model=MW213$model)
    #  [1]  0.0000000  0.0000000  0.0000000 30.0000000 30.0000000  0.3465736
    #  [7]  0.5493061  0.3465736 -0.5493061  3.0000000  2.0000000
}
\keyword{misc}
