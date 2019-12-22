\name{nMm2par}
\alias{nMm2par}
\title{Multivariate Normal Mixture Model to parameter for MLE}
\description{
  From a \code{"\link{norMmix}"}(-like) object, return the numeric
  parameter vector in our MLE parametrization.
}
\usage{
nMm2par(obj, trafo = c("clr1", "logit"),
        model = c("EII", "VII", "EEI", "VEI", "EVI",
                  "VVI", "EEE", "VEE", "EVV", "VVV"),
        meanFUN = mean)
}
\arguments{
  \item{obj}{a \code{\link{list}} containing \describe{
      \item{\code{sig}:}{covariance matrix array,}
      \item{\code{mu}:}{mean vector matrix,}
      \item{\code{w}:}{= weights,}
      \item{\code{k}:}{= number of components,}
      \item{\code{p}:}{= dimension}
    }
  }
  \item{trafo}{a \code{\link{character}} string specifying the transform to
    use for the weights \code{w} (\eqn{= \pi_j}{= pi[j]}), currently either
    \code{"clr1"} or \code{"logit"}.}
  \item{model}{a \code{\link{character}} string specifying the (Sigma)
    model, one of those listed above.}
  \item{meanFUN}{a \code{\link{function}} to compute a mean (of variances typically).}
}
\details{
    This transformation forms a vector from the parameters of a normal
    mixture. These consist of weights, means and covariance matrices.
    Weights are transformed according to 'trafo' param; means are
    unchanged.

    Cov mats are given as D and L from the LDLt decomposition
}
\seealso{
    \code{\link{nc2p}} simplified wrapper function.
    
    the \emph{inverse} function of \code{nMm2par()} is \code{\link{par2nMm}()}.
}

\examples{
A <- MW24
nMm2par(A, trafo = "clr1", model = A$model)
# [1] -0.3465736  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
# [7] -2.3025851
}
