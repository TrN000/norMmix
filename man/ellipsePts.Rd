\name{ellipsePts}
\alias{ellipsePts}
\title{Compute Points on Bivariate Gaussian Confidence Ellipse}
\description{
  From 2-dimensional mean vector \code{mu}\eqn{= \mu} and 2x2 covariance matrix
  \code{sigma}\eqn{= \Sigma}, compute \code{npoints} equi-angular points on
  the \code{1-alpha}\eqn{ = 1-\alpha} confidence ellipse of bivariate
  Gaussian (normal) distribution \eqn{\mathcal{N}_2(\mu,\Sigma)}{N(\mu, \Sigma)}.
}
\usage{
ellipsePts(mu, sigma, npoints, alpha = 0.05, r = sqrt(qchisq(1 - alpha, df = 2)))
}
\arguments{
  \item{mu}{mean vector (\code{\link{numeric}} of length 2).}
  \item{sigma}{2x2 \code{\link{matrix}}, the covariance matrix.}
  \item{npoints}{integer specifying the number of points to be computed.}
  \item{alpha}{confidence level such that the ellipse should contain
    1-alpha of the mass.}
  \item{r}{radius of the ellipse, typically computed from \code{alpha}, via
    the default value.}
}
\value{
  a numeric matrix of dimension \code{npoints x 2}, containing the
  x-y-coordinates of the ellipse points.
}
\author{Martin Maechler}
\note{
  This has been inspired by package \CRANpkg{mixtools}'s \code{\link[mixtools]{ellipse}()}
  function.
}
%% \seealso{

%% }
\examples{
xy <- ellipsePts(c(10, 100), sigma = cbind(c(4, 7), c(7, 28)),  npoints = 20)
plot(xy, type = "b", col=2, cex=2,
     main="ellipsePts(mu = (10,100), sigma, npoints = 20)")
points(10, 100, col=3, cex=3, pch=3)
text  (10, 100, col=3, expression(mu == "mu"), adj=c(-.1, -.1))

stopifnot(is.matrix(xy), dim(xy) == c(20, 2))
}
\keyword{distribution}
