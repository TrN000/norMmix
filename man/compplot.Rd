\name{compplot}
\alias{compplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    composition plot
}
\description{
    Creates a \code{\link{massplot}} like plot. Takes two \code{\link{massbic}}
    arrays and overlays their plots.
}
\usage{
compplot(f, g, h=NULL, main = NULL, adj = 1/dim(f)[3],
         col = nMmcols[1:3], mar = 0.1 + c(1.4, 2, 3, 1), 
         compnames = c("clara", "mclVVV", "Mclust"), 
         oma=c(7,0,2.8,0), pars=TRUE,
      ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{f}{
        Array as from \code{\link{massbic}} or \code{\link{massbicm}}
    }
    \item{g}{
        Array as from \code{\link{massbic}} or \code{\link{massbicm}}
    }
    \item{h}{
        Optional third array.
    }
    \item{main}{
        Character string to be used as title of the plot.
    }
    \item{adj}{
        Adjustment for alpha value of the plot color.
    }
    \item{col}{
        Color of plot.
    }
    \item{mar}{
        Margin space adjustment
    }
    \item{compnames}{
        Names to be used for legend. Set to NULL to supress.
    }
    \item{oma}{
        Margin parameter. See par.
    }
    \item{pars}{
        Logical, whether to display number of parameters on top of plot.
    }
    \item{...}{
        further arguments passed on to \code{\link{matplot}}
    }
}
\details{
    The intended use for this function is to run \code{\link{massbic}} and 
    \code{\link{massbicm}} with the same arguments, \code{string} and \code{DIR}
    and feeding the results into \code{compplot}, producing line plots, that 
    can be compared against each other. Currently uses \code{rainbow(20)} to 
    generate plot colors. Red, orange and green correspond to \code{f} and 
    cyan, blue and purple to \code{g}.
}
\value{
    No return value. Intended side effect is the generation of a plot.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
    Nicolas Trutmann
}
\note{
    While possible to be used as standalone, the function \code{\link{epfl}} 
    allows for more complex generation of PDFs and is the intended use for 
    simulations.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{massplot}},
    \code{\link{epfl}},
    \code{\link{massbic}},
    \code{\link{massbicm}}
}
\examples{
## TODO: add example
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
