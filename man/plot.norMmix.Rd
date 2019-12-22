\name{plot.norMmix}
\alias{plot.norMmix}
\title{Plot Method for "norMmix" Objects}
\description{
  This is the S3 method for plotting \code{"\link{norMmix}"} objects.
}
\usage{
\S3method{plot}{norMmix}(x, y=NULL, ...)
}
\arguments{
  \item{x}{an \R object inheriting from \code{"\link{norMmix}"}.}
  \item{y}{further data matrix, first 2 columns will be plotted by \code{"\link{points}"}}
  \item{\dots}{further arguments to be passed to \code{"\link{plot}"}}
}
\value{
  \code{plot.norMmix} returns invisibly coordinates of bounding ellipses of distribution.
}
\examples{
plot(MW212) ## and add a finite sample realization:
points(rnorMmix(n=500, MW212))

## or:
x <- points(rnorMmix(n=500, MW212))
plot(MW212, x)
}
\keyword{hplot}
