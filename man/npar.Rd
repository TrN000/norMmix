\name{npar}
\title{Degrees of freedom of (Fitted) Multivariate Normal Mixtures}
\alias{npar}
\description{
    This function is generic; method functions can be written to
    handle specific classes of objects. The following classes have
    methods written for them:

    norMmix

    norMmixMLE

    fittednorMmix
}
\usage{
npar(object, \dots)
}
\arguments{
  \item{object}{any \R object from the list in the \sQuote{Description}.}
  \item{\dots}{potentially further arguments for methods;
    Currently, none of the methods for the listed classes do have such.}
}
%% \details{
%%     % TODO: description of behaviour.
%% }
\value{
  Depending on \code{object} : %% FIXME? other *syntax* to list value for different *methods*
  \item{norMmix }{integer number.}
  \item{norMmixMLE }{integer number.}
  \item{fittednorMmix }{integer\code{\link{matrix}} with
    \code{\link{dimnames}} set to k and models.}
}
\author{Nicolas Trutmann}
\seealso{
  \code{\link{norMmix}}, \code{\link{norMmixMLE}}.
}
\examples{
methods(npar) # list available methods

npar(MW213)
}
