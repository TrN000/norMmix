% Code in ../R/norMmixMLE.R
\name{norMmixMLE}
\title{Maximum Likelihood Estimation for Multivariate Normal Mixture Models}
\alias{norMmixMLE}
\alias{ssClaraL}
\description{
  Direct Maximum Likelihood Estimation (MLE) for multivariate normal
  mixture models \code{"\link{norMmix}"}.  Starting from a
  \code{\link[cluster]{clara}} (package \CRANpkg{cluster}) clustering plus
  one M-step, or alternatively from the default start of (package)
  \CRANpkg{mclust}, perform direct likelihood maximization via \code{\link{optim}()}.
}
\usage{
norMmixMLE(x, k,
           model = c("EII", "VII", "EEI", "VEI", "EVI",
                     "VVI", "EEE", "VEE", "EVV", "VVV"),
           ini,
           ll = c("nmm", "mvt"),
           method = "BFGS", maxit = 100, trace = 2,
           optREPORT = 10, reltol = sqrt(.Machine$double.eps),
 	   \dots)

ssClaraL(n, k, p)
}
\arguments{
  \item{x}{numeric [n x p] matrix}
  \item{k}{positive number of components}
  \item{model}{a \code{\link{character}} string, specifying the model (for
    the \code{k} covariance matrices) to be assumed.}
  \item{ini}{a \code{\link{function}}, that takes arguments \code{x} and 
      \code{k} and returns a clustering index; a vector of length \code{ncol(x)},
      with entries in \code{1:k}.}
  \item{ll}{a string specifying the method to be used for the likelihood
    computation;  the default, \code{"nmm"} uses \code{\link{llnorMmix}()},
    whereas \code{"mvt"} uses \code{\link{llmvtnorm}()} which is based on
    the MV normal density from package \CRANpkg{mvtnorm}.}

  \item{method, maxit, trace, optREPORT, reltol, \dots}{arguments for tuning the
    optimizer \code{\link{optim}(*, method=method, control = list(...))}.}

  \item{n,p}{matrix dimensions \code{\link{nrow}(x)} and \code{\link{ncol}(x)}.}
}
\value{
  \code{norMmixMLE} returns an object of \code{\link{class}}
  \code{"norMmixMLE"} which is a \code{\link{list}} with components
  \item{norMmix}{the \code{"\link{norMmix}"} object corresponding to the
    specified model and the fitted (MLE) parameter vector.}
  \item{optr}{the [r]eturn value of \code{\link{optim}()}.}
  \item{npar}{the number of free parameter, a function of \eqn{(p, k, model)}.}
  \item{n}{the sample size, i.e., the number of observations or rows of \code{x}.}
  \item{cond}{the result of \code{parcond(..)}, that is the ratio of sample size 
      over parameter count.}% fixme
}
\details{
  Uses clara() and one M-step from EM-algorithm to initialize parameters
  after that uses general optimizer optim() to calculate ML.
}
\examples{
    str(MW214)
    set.seed(105)
    x <- rnorMmix(1000, MW214)
    ## Fitting assuming we know the true parametric model
    fm1 <- norMmixMLE(x, k = 6, model = "VII", ini=norMmix:::clarafunc)
    if(interactive()) ## Fitting "wrong" overparametrized model: typically need more iterations:
    fmW <- norMmixMLE(x, k = 7, model = "VVV", maxit = 200, ini=norMmix:::clarafunc)# default maxit=100 is often too small
}
