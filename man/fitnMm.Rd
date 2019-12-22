% Code in ../R/fit.R
\name{fitnMm}
\alias{fitnMm}
\title{Fit Several Normal Mixture Models to a Dataset}
\description{
  \code{fitnMm()} fits several multivariate normal mixture models
  (\code{norMmix}) to the data set \code{x}, and returns (among other info)
  a \code{\link{list}} of fitted \code{\link{norMmix}} objects.
}
\usage{
fitnMm(x, k, models = 1:10,
            trafo = c("clr1", "logit"),
            ll = c("nmm", "mvt"), 
            savdir = NULL, name = NULL,
         \dots)
}
\arguments{
  \item{x}{data matrix, rows are observations and columns are variables}
  \item{k}{vector of positive integers, indicating the \emph{number} of
    mixture components to fit.}
  \item{models}{vector of integers from \code{1:10}, indexing models from\preformatted{
      c("EII", "VII", "EEI", "VEI", "EVI",
        "VVI", "EEE", "VEE", "EVV", "VVV")}
    which are to be fit.}
  \item{trafo}{a \code{\link{character}} string specifying the transform to
    use for the weights \code{w} (\eqn{= \pi_j}{= pi[j]}), currently either
    \code{"clr1"} or \code{"logit"}.}
  \item{ll, \dots}{further arguments passed to function \code{\link{norMmixMLE}()}.}
  \item{savdir}{valid file path to save directory}
  \item{name}{name to be used when saving return value as RDS file}
}
\value{
  \code{fitnMm()} returns a \code{\link{list}} with components
  \item{nMm}{a \code{\link{list}} containing all fitted models}
  \item{models}{\code{\link{character}} vector of model (names) that were
    fitted.}
  \item{n}{number of observations of x}
  \item{p}{number of variables of x}
}
%% \details{
%% }
\note{
    Given an object \code{r} of class \code{fittednorMmix},
    for use with \code{\link{massbic}}, do:

    saveRDS(list(fit=r))
}
\seealso{
  \code{\link{norMmixMLE}()} which is called \code{length(k) *
    length(models)} times.
}
\examples{
x <- rnorMmix(500, MW21)
fitnMm(x, 1:2, models=1:4) ## will fit models 1:4 with 1:3 components
}
