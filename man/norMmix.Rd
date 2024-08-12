% code in  ../R/norMmix.R
\name{norMmix}
\alias{norMmix}
\title{Constructor for Multivariate Normal Mixture Objects}
\description{
  \code{norMmix} creates a multivariate normal (aka Gaussian) mixture
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
    \item{mu}{matrix of means, or a vector in which case \eqn{k = 1} is
      assumed.  Otherwise use \code{\link{as.matrix}(mu)}.}
    \item{Sigma}{NULL, number, numeric, vector (length = k), matrix (dim = p x k), or
        array (p x p x k). See details.
    }
    \item{weight}{weights of mixture model components}
    \item{name}{gives the option of naming mixture}
    \item{model}{see \sQuote{Details}}
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
    non-negative, summing to one: \code{\link{sum}(weight) == 1}.}
  \item{k}{integer, the number of components}
  \item{dim}{integer, the dimension \eqn{p}.}
}
\details{
  \code{model} must be specified by one of the (currently 10)
  \code{\link{character}} strings shown in the default.
  (In a future version, \code{model} may become \emph{optional}).% TODO - "smart"

  \code{norMmix} as a few nifty ways of constructing simpler matrices from
  smaller givens. This happens according to the dimension of the given
  value for the Sigma argument:
  \describe{
    \item{0.}{for a single value \code{d} or \code{NULL}, \code{norMmix()} assumes
      all covariance matrices to be diagonal with entries \code{d} or \code{1}, respectively. }
    \item{1.}{for a vector \code{v}, \code{norMmix} assumes all matrices to
      be diagonal with the i-th matrix having diagonal entries \code{v[i]}.}
    \item{2.}{for a matrix \code{m}, \code{norMmix} assumes all matrices to
      be diagonal with diagonal vector \code{m[,i]}, i.e., it goes by
      columns.}
    \item{3.}{an array is assumed to be the covariance matrices, given
      explicitly.}
  }


  FIXME ... give "all" the details ... (from Bachelor's thesis ???)
  %% ~/Betreute-Arbeiten/NicolasTrutmann/BSc_thesis+MM/chapter1.Rnw
  %% ~/Betreute-Arbeiten/NicolasTrutmann/GH_BachelorArbeit/Bachelorthesis/chapter1.Rnw


%%  The '\tabular' command takes two arguments.  The first gives for each
%% of the columns the required alignment ('l' for left-justification, 'r'
%% for right-justification or 'c' for centring.)  The second argument
%% consists of an arbitrary number of lines separated by '\cr', and with
%% fields separated by '\tab'.  For example:

%%        \tabular{rlll}{
%%          [,1] \tab Ozone   \tab numeric \tab Ozone (ppb)\cr
%%          [,2] \tab Solar.R \tab numeric \tab Solar R (lang)\cr
%%          [,3] \tab Wind    \tab numeric \tab Wind (mph)\cr
%%          [,4] \tab Temp    \tab numeric \tab Temperature (degrees F)\cr
%%          [,5] \tab Month   \tab numeric \tab Month (1--12)\cr
%%          [,6] \tab Day     \tab numeric \tab Day of month (1--31)
%%        }

%% There must be the same number of fields on each line as there are
%% alignments in the first argument, and they must be non-empty (but can
%% contain only spaces).  (There is no whitespace between '\tabular' and
%% the first argument, nor between the two arguments.)


}
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
## Some of the "MW" objects : % --> ../R/zmarrwandnMm.R

# very simple 2d:
M21 <- norMmix(mu = cbind(c(0,0)), # 2 x 1 ==> k=2, p=1
               Sigma = 1, model = "EII")
stopifnot(identical(M21, # even simpler, Sigma = default :
                    norMmix(mu = cbind(c(0,0)), model = "EII")))

m2.2 <- norMmix(mu = cbind(c(0, 0), c(5, 0)), Sigma = c(1, 10),
                weight = c(7,1)/8, model = "VEI")

m22 <- norMmix(
    name = "one component rotated",
    mu = cbind( c(0,0) ),
    Sigma = array(c(55,9, 9,3), dim = c(2,2, 1)),
    model = "EVV")
stopifnot( all.equal(MW22, m22) )

m213 <- norMmix(
    name = "#13 test VVV",
    weight = c(0.5, 0.5),
    mu = cbind( c(0,0), c(30,30) ),
    Sigma = array(c( 1,3,3,11, 3,6,6,13 ), dim=c(2,2, 2)),
    model = "VVV")
stopifnot( all.equal(MW213, m213) )
str(m213)

m34 <- norMmix(
    name = "#4 3d VEI",
    weight = c(0.1, 0.9),
    mu = matrix(rep(0,6), 3,2),
    Sigma = array(c(diag(1:3), 0.2*diag(3:1)), c(3,3, 2)),
    model = "VVI" )
stopifnot( all.equal(MW34, m34) )


}
\keyword{distribution}
