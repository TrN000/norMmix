\name{clr1}
\title{Centered Log Ratio Transformation and Inverse}
\alias{clr1}
\alias{clr1inv}
\description{% --> code in ../R/param.R
  The centered log ratio transformation is Maechler's solution to allowing
  unconstrained mixture weights optimization.

  It has been inspired by Aitchison's \bold{c}entered \bold{l}og \bold{r}atio,
  see also CRAN package \CRANpkg{compositions}' \code{\link[compositions]{clr}()}, and typically
  other references on modelling proportions.
}
\usage{
clr1(w)
clr1inv(lp)
}
\arguments{
  \item{w}{numeric vector of length \eqn{k}, say, of mixture weights, i.e.,
    non-negative and summing to one.}
  \item{lp}{numeric vector of length \eqn{k-1} clr-transformed weights.}
}
\details{
  Aitchison's clr transformation is slightly different, as it does
  \emph{not} drop one coordinate, as we do.  Hence the extra \sQuote{1} in
  the name of our version.
}
\value{
  a numeric vector of length \eqn{k-1} or \eqn{k}, see above.
}
\references{
  Aitchison, J., 1986.
  \emph{The Statistical Analysis of Compositional Data}
  Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK).

  More in the CRAN package \CRANpkg{compositions} vignette
  \file{UsingCompositions.pdf} % FIXME CRAN url
}
\author{Martin Maechler}
%% \note{
%% }
\seealso{
  The first implementation of these was in \CRANpkg{nor1mix}, June 2019, in its
  \code{\link[nor1mix]{par2norMix}()} and
  \code{\link[nor1mix]{nM2par}()} functions.% --> code  in ../../nor1mix/R/norMix.R
}
\examples{
## Apart from error checking and very large number cases, the R implementation is simply
..clr1 <- function (w) {
  ln <- log(w)
  ln[-1L] - mean(ln)
}

## and its inverse
..clr1inv <- function(lp) {
  p1 <- exp(c(-sum(lp), lp))
  p1/sum(p1)
}

lp <- clr1( (1:3)/6 )
clr1inv(lp)
stopifnot(all.equal(clr1inv(lp), (1:3)/6))

for(n in 1:100) {
   k <- 2 + rpois(1, 3) # #{components}
   lp <- rnorm(k-1) # arbitrary unconstrained
   ## clr1() and clr1inv() are inverses :
   stopifnot(all.equal(lp, clr1(clr1inv(lp))))
}

wM <- clr1inv(c(720,720,720))
w2 <- clr1inv(c(720,718,717))
stopifnot(is.finite(wM), all.equal(wM, c(0, 1/3, 1/3, 1/3))
        , is.finite(w2), all.equal(w2, c(0, 0.84379473, 0.1141952, 0.042010066))
         )
}
\keyword{math}
