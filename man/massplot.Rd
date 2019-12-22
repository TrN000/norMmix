\name{massplot}
\alias{massplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    plot from massbic
%%  ~~function to do ... ~~
}
\description{
    plots the result of \code{\link{massbic}} and \code{\link{massbicm}}
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
massplot(f, main = "unnamed", adj = exp(-0.002 * size), 
         col = nMmcols[1], mar = 0.1 + c(1.4, 2, 3, 1), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{f}{
        Result of \code{\link{massbic}}
    }
    \item{main}{
        Character string, title of plot.
    }
    \item{adj}{
        Alpha adjustment for plot color.
    }
    \item{col}{
        Plot color.
    }
    \item{mar}{
        A numerical vector of the form ‘c(bottom, left, top, right)’.
    }
    \item{...}{
        Further parameters, passed to \code{\link{matplot}}
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (f, main = "unnamed") 
{
    ran <- extendrange(f)
    size <- dim(f)[3]
    cl <- as.numeric(dimnames(f)$components)
    p <- attr(f, "dims")
    adj <- exp(-0.002 * size)
    models <- mods()
    op <- sfsmisc::mult.fig(mfrow = c(4, 5), main = main, mar = 0.1 + 
        c(2, 4, 4, 1))
    for (i in 1:10) {
        if (!is.null(p)) {
            matplot(f[, i, ], lty = 1, col = adjustcolor(rainbow(10)[i], 
                adj), type = "l", ylim = ran, main = models[i])
            axis(3, at = seq_along(cl), labels = npar(cl, p[1], 
                models[i]))
        }
        else {
            matplot(f[, i, ], lty = 1, col = adjustcolor(rainbow(10)[i], 
                adj), main = models[i], type = "l", ylim = ran)
        }
    }
    for (i in 1:10) {
        boxplot(t(f[, i, ]), lty = 1, col = adjustcolor(rainbow(10)[i], 
            0.4), main = models[i], type = "l", ylim = ran)
    }
    par(op$old.par)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
