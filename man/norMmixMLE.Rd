% Code in ../R/norMmixMLE.R
\name{norMmixMLE}
\title{Maximum Likelihood Estimation for Multivariate Normal Mixtures}
\alias{norMmixMLE}
\alias{claraInit}
\alias{mclVVVinit}
\alias{ssClara2kL}
\description{
  Direct Maximum Likelihood Estimation (MLE) for multivariate normal
  mixture models \code{"\link{norMmix}"}.  Starting from a
  \code{\link[cluster]{clara}} (package \CRANpkg{cluster}) clustering plus
  one M-step by default, or alternatively from the default start of (package)
  \CRANpkg{mclust}, perform direct likelihood maximization via \code{\link{optim}()}.
}
\usage{
norMmixMLE(x, k,
           model = c("EII", "VII", "EEI", "VEI", "EVI",
                     "VVI", "EEE", "VEE", "EVV", "VVV"),
           initFUN = claraInit,
           ll = c("nmm", "mvt"),
           keep.optr = TRUE, keep.data = keep.optr,
           method = "BFGS", maxit = 100, trace = 2,
           optREPORT = 10, reltol = sqrt(.Machine$double.eps),
 	   \dots)

claraInit(x, k, samples = 128,
          sampsize = ssClara2kL, trace)
mclVVVinit(x, k, \dots)

ssClara2kL(n, k, p)
}
\arguments{
  \item{x}{numeric [n x p] matrix}
  \item{k}{positive number of components}
  \item{model}{a \code{\link{character}} string, specifying the model (for
    the \code{k} covariance matrices) to be assumed.}
  \item{initFUN}{a \code{\link{function}}, that takes arguments \code{x} and
    \code{k} and returns a clustering index; a vector of length \eqn{p = }\code{ncol(x)},
      with entries in \code{1:k}.}
  \item{ll}{a string specifying the method to be used for the likelihood
    computation;  the default, \code{"nmm"} uses \code{\link{llnorMmix}()},
    whereas \code{"mvt"} uses \code{\link{llmvtnorm}()} which is based on
    the MV normal density from package \CRANpkg{mvtnorm}.}

  \item{keep.optr, keep.data}{\code{\link{logical}}, each indicating of the
    optimization result (from \code{\link{optim}()}, currently), or the data
    \code{x} respectively, should be saved as part of the result (function
    \sQuote{value}, see also below).}

  \item{method, maxit, optREPORT, reltol}{arguments for tuning the
    optimizer \code{\link{optim}(*, method=method, control = list(...))}.}

  %% both norMmixMLE() and  <some>():
  \item{trace}{\describe{
      \item{in \code{norMmixMLE()}:}{passed to \code{optim(*, control=..)}, see above.}
      \item{in \code{claraInit()}: }{a non-negative integer indicating how
	much \code{\link{clara}()} calls should be traced.}
  }}
  \item{\dots}{\describe{
      \item{in \code{norMmixMLE()}:}{passed to \code{optim(*, control=..)}, see above.}
      \item{in \code{mclVVVinit()}:}{further arguments passed to (package
	\CRANpkg{mclust}) function \code{\link[mclust]{hcVVV}()}.}
  }}

  %% claraInit():
  \item{samples}{the number of subsamples to take in
    \code{\link[cluster]{clara}()}, package \CRANpkg{cluster}, see its help.}
  \item{sampsize}{the sample size to take in
    \code{\link[cluster]{clara}()}, package \CRANpkg{cluster}.
    Here, can be a positive integer \emph{or}, as by default, a
    \code{\link{function}} with arguments \code{(n,k,p)}.}

  %% ssClara2kL() :
  \item{n,p}{matrix dimensions \code{\link{nrow}(x)} and \code{\link{ncol}(x)}.}
}
\value{
  \code{norMmixMLE} returns an object of \code{\link{class}}
  \code{"norMmixMLE"} which is a \code{\link{list}} with components

  \item{norMmix}{the \code{"\link{norMmix}"} object corresponding to the
    specified model and the fitted (MLE) parameter vector.}
  \item{optr}{(if \code{keep.optr} is true:) the [r]eturn value of
    optimization, i.e., currently, \code{\link{optim}()}.}
  \item{npar}{the number of free parameters, a function of \eqn{(p, k, model)}.}
  \item{n}{the sample size, i.e., the number of observations or rows of \code{x}.}
  \item{cond}{the result of \code{parcond(..)}, that is the ratio of sample size
    over parameter count.}% fixme
  \item{x}{(if \code{keep.optr} is true:) the \eqn{n \times p}{n x p} data matrix.}
}
\details{
  By default, \code{initFUN=claraInit}, uses \code{\link[cluster]{clara}()}
  and one M-step from EM-algorithm to initialize parameters
  after that uses general optimizer \code{\link{optim}()} to calculate the MLE.
}
\examples{
MW214
set.seed(105)
x <- rnorMmix(1000, MW214)

\donttest{## Fitting, assuming to know the true model (k=6, "VII")
fm1  <- norMmixMLE(x, k = 6, model = "VII", initFUN=claraInit)
fm1 # {using print.norMmixMLE() method}
fm1M <- norMmixMLE(x, k = 6, model = "VII", initFUN=mclVVVinit)

## Fitting "wrong" overparametrized model: typically need more iterations:
fmW <- norMmixMLE(x, k = 7, model = "VVV", maxit = 200, initFUN=claraInit)
## default maxit=100 is often too small    ^^^^^^^^^^^
}% dont.. %

x <- rnorMmix(2^12, MW51)
fM5 <- norMmixMLE(x, k = 4) # k = 3 is sufficient
fM5
c(logLik = logLik(fM5), AIC = AIC(fM5), BIC = BIC(fM5))
plot(fM5, show.x=FALSE)
plot(fM5, lwd=3, pch.data=".")

\donttest{# this takes several seconds
 fM5big <- norMmixMLE(x, model = "VVV", k = 4, maxit = 300) # k = 3 is sufficient
 summary(warnings())
 fM5big ; c(logLik = logLik(fM5big), AIC = AIC(fM5big), BIC = BIC(fM5big))
 plot(fM5big, show.x=FALSE)
}%dont.. %
}
