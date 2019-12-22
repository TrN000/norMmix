\name{rnorMmix}
\alias{rnorMmix}
\title{Random Sample from Multivariate Normal Mixture Distribution}
\description{
  Draw \code{n} (p-dimensional) observations randomly from the multivariate normal
  mixture distribution specified by \code{obj}.
}
\usage{
rnorMmix(n, obj, index = FALSE, permute = TRUE)
}
\arguments{
  \item{n}{sample size, non-negative.}
  \item{obj}{a \code{"\link{norMmix}"} object}
  \item{index}{Logical, store the clustering information as first column}
  \item{permute}{Logical, indicating if the observations should be randomly
    permuted after creation \dQuote{cluster by cluster}.}
}
%% \details{
%% }
\value{
  n p-dimensional observations, as numeric \eqn{n \times p}{n * p} matrix.
}
%% \references{
%% }
\author{Nicolas Trutmann}
\seealso{
  \code{\link{rmultinom}}
}
\examples{
x <- rnorMmix(500, MW213)
plot(x)
x <- rnorMmix(500, MW213, index=TRUE)
plot(x[,-1], col=x[,1]) ## using index column to color components
}
