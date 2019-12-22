\name{massbicm}
\alias{massbicm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Do mclust along .rds files from fitnMm
%%  ~~function to do ... ~~
}
\description{
    \code{massbicm} applies \code{\link{Mclust}} along an existing
    \code{fittednorMmix} object, assumed to be in an RDS file in a list,
    named 'fit'.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
massbicm(string, DIR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{string}{
%%     ~~Describe \code{string} here~~
}
  \item{DIR}{
%%     ~~Describe \code{DIR} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (string, DIR) 
{
    nm <- readRDS(file.path(DIR, string[1]))
    cl <- nm$fit$k
    mo <- nm$fit$models
    valm <- array(0, lengths(list(cl, mo, string)))
    dims <- vector(mode = "integer", length = length(string))
    for (i in 1:length(string)) {
        nm <- readRDS(file.path(DIR, string[i]))
        x <- nm$fit$x
        valm[, , i] <- mclust::Mclust(x, G = cl, modelNames = mo)$BIC
        dims[i] <- nm$fit$p
    }
    dimnames(valm) <- list(components = cl, models = mo, files = string)
    attr(valm, "dims") <- dims
    -valm
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
