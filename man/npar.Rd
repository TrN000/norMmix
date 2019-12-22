\name{npar}
\alias{npar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Extract degrees of freedom from objects of the \code{\link{norMmix}} package
}
\description{
    This function is generic; method functions can be written to 
    handle specific classes of objects. The following classes have
    methods written for them:

    norMmix

    norMmixMLE

    fittednorMmix
}
\usage{
    npar(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
      Any object from the list in the Description.
    }
  \item{...}{
      In place to leave room for further arguments.
      None of the methods for the listed classes take arguments beyond 
      \code{object}.
    }
}
\details{
    % TODO: description of behaviour.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    \item{norMmix}{Integer vector of length 1.}
    \item{norMmixMLE}{Integer vector of length 1.}
    \item{fittednorMmix}{Integer matrix with dimnames set to k and models.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
    Nicolas Trutmann
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(fSMI.12)
npar(fSMI.12)

npar(MW213)
}
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
