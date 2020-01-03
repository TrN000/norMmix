\name{llmvtnorm}
\alias{llmvtnorm}
\title{Log-Likelihood of Multivariate Normal Mixture Relying on \code{mvtnorm::dmvnorm}}
\description{
  Compute the log-likelihood of a multivariate normal mixture, by calling
  \code{\link[mvtnorm]{dmvnorm}()} (from package \CRANpkg{mvtnorm}).
}
\usage{
llmvtnorm(par, x, k,
          model = c("EII", "VII", "EEI", "VEI", "EVI",
                    "VVI", "EEE", "VEE", "EVV", "VVV"))
}
\arguments{
  \item{par}{parameter vector as calculated by nMm2par}
  \item{x}{numeric data \code{\link{matrix}} (of dimension \eqn{n \times p}{n x p}).}
  \item{k}{number of mixture components.}
  \item{model}{assumed model of the distribution}
}
\value{
  returns the log-likelihood (a number) of the specified model for the data
  (\eqn{n} observations) \code{x}.
}
%% \references{
%% }
\seealso{
  \code{\link[mvtnorm]{dmvnorm}()} from package \CRANpkg{mvtnorm}.  Our own
  function, returning the same: \code{\link{llnorMmix}()}.
}
\examples{
set.seed(1); x <- rnorMmix(50, MW29)
para <- nMm2par(MW29, model=MW29$model)

llmvtnorm(para, x, 2, model=MW29$model)
# [1] -236.2295
}
