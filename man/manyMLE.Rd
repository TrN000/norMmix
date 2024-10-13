\name{manyMLE}
\alias{manyMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
manyMLE
}
\description{
Run \code{\link{norMmixMLE}} over various component sizes and model types.
}
\usage{
manyMLE(x, k, models = 1:10, trafo = c("clr1", "logit"), ll = c("nmm", "mvt"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Data set to fit models to. Expects the same shape as
    \code{\link{norMmixMLE}}, that is, a \code{[n x p]} numeric matrix.
  }
  \item{k}{
    Integer vector. list of cluster sizes to fit to.
  }
  \item{models}{
    Integer vector. This vector indexes the list of models.
  }
  \item{ll}{
    Log-likelihood function. Options are \code{"nmm"} for our own
    \code{\link{llnorMmix}} or \code{"mvt"} for
    \code{\link[mvtnorm]{dmvnorm}} from the \CRANpkg{mvtnorm} package.
  }
  \item{\dots}{
    Additional arguments to be passed to \code{\link{norMmixMLE}}.
    All arguments from there are valid, except for
    \code{x},\code{k},\code{model} and \code{ll}.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  List of class \code{manyMLE}:

  \item{nMm}{Doubly indexed list of return values of
    \code{norMmixMLE}. First index is along (as in
    \code{\link{seq_along}}) \code{k}; second index is along
    \code{models}.
  }
  \item{nMmtime}{
    timings from \code{\link{system.time}}. Indexed like \code{nMm}.
  }
  \item{k}{The \code{k} given as argument.}
  \item{m}{Vector of models used.}
  \item{n}{Number of datapoints.}
  \item{p}{Dimension of dataset.}
  \item{x}{Dataset.}
}
\author{
  Nicolas Trutmann
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{norMmixMLE}}
}
\examples{
# generate sample data
x <- rnorMmix(100, MW21)

# run fitting
fitted <- manyMLE(x, 1:2, 1:2)

# > fitted$nMm
#   EII          VII         
# 1 norMmixMLE,6 norMmixMLE,6
# 2 norMmixMLE,6 norMmixMLE,6

# custom print method which includes best fit and time.
print(fitted)
}
