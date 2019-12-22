% code in  ../R/norMmix.R
\name{norMmix}
\alias{norMmix}
\title{Constructor for Multivariate Normal Mixture Objects}
\description{
  \code{norMmix()} creates a multivariate normal (aka Gaussian) mixture
  object, conceptually a mixture of \eqn{k} multivariate
  (\eqn{p}-dimensional) Gaussians
  \eqn{\mathcal{N}(\mu_j, \Sigma_j)}{N(mu_j, Sigma_j)}, for \eqn{j=1, \dots, k}.
}
\usage{
norMmix(mu, Sigma = NULL, weight = rep(1/k, k), name = NULL,
        model = c("EII", "VII", "EEI", "VEI", "EVI",
                  "VVI", "EEE", "VEE", "EVV", "VVV"))
}
\arguments{
    \item{mu}{matrix of means. should mu be a vector it will assume k=1
              to circumvent this behavoiur use as.matrix(mu) beforehand}
    \item{Sigma}{array of covariance matrices}
    \item{weight}{weights of mixture model components}
    \item{name}{gives the option of naming mixture}
    \item{model}{see desc}
}
\value{
  currently, a \code{\link{list}} of class \code{"norMmix"}, with
  a \code{name} attribute and components
  \item{\code{model}}{three-letter \code{\link{character}} string, specifying the
    Sigma-parametrization}
  \item{\code{mu}}{(p x k) matrix of component means \code{mu[,j]},
    \eqn{j=1,\dots,k}.}
  \item{\code{Sigma}}{(p x p x k) array of component Covariance matrices
    \code{Sigma[,,j]}.}
  \item{\code{weight}}{p-vector of mixture probability weights;
    non-negative, summing to one: \code{sum(weight) == 1}.}
  \item{k}{integer, the number of components}
  \item{dim}{integer, the dimension \eqn{p}.}
}
%% \details{
%% }
\references{
  __ TODO __
}
\author{Nicolas Trutmann}
\seealso{
  \code{\link{norMmixMLE}()} to fit such mixture models to data (an \eqn{n
    \times p}{n * p} matrix).

  \dQuote{Marron-Wand}-like examples (for testing, etc), such as
  \code{\link{MW21}}.
}
\examples{
# TODO
}
\keyword{distributions}
