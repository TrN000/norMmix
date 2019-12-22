% Code in ../R/llnorMmix.R
\name{llnorMmix}
\title{Log-likelihood of parameter vector given data}
\alias{llnorMmix}
\description{
  Calculates log-likelihood of a dataset, tx, given a normal mixture model as
  specified by a parameter vector.
    A parameter vector can be obtained by applying \code{\link{nMm2par}} to a
    \code{\link{norMmix}} object.
}
\usage{
llnorMmix(par, tx, k,
          trafo = c("clr1", "logit"),
          model = c("EII", "VII", "EEI", "VEI", "EVI",
                    "VVI", "EEE", "VEE", "EVV", "VVV"))
}
\arguments{
  \item{par}{parameter vector}
  \item{tx}{\emph{Transposed} numeric data matrix, i.e. \code{tx := t(x)}
    is of dimension \eqn{p \times n}{p x n}; its rows are variables and
    columns are observations.}
  \item{k}{number of mixture components.}
  \item{trafo}{a \code{\link{character}} string specifying the transform to
    use for the weights \code{w} (\eqn{= \pi_j}{= pi[j]}), currently either
    \code{"clr1"} or \code{"logit"}.}
  \item{model}{assumed distribution model of normal mixture}
}
\value{
  returns the log-likelihood (a number) of the specified model for the data
  (\eqn{n} observations) \code{x}.
}
\seealso{
  Our alternative function \code{\link{llmvtnorm}()} (which is based on
  \code{dmvnorm()} from package \CRANpkg{mvtnorm}).
}
\examples{
set.seed(1); tx <- t(rnorMmix(50, MW29))
para <- nMm2par(MW29, model=MW29$model)

llnorMmix(para, tx, 2, model=MW29$model)
# [1] -236.2295
}

