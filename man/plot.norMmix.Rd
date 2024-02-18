\name{plot.norMmix}
\title{Plot Method for "norMmix" Objects}
\alias{plot.norMmix}
\alias{plot.norMmixMLE}
\alias{plot.fittednorMmix}
\alias{plot2d}
\alias{plotnd}
\alias{Trubetskoy10}
\description{
  This is the S3 method for plotting \code{"\link{norMmix}"} objects.
}
\usage{
\S3method{plot}{norMmix}(x, y=NULL, ...) %%-- TODO all arguments
\S3method{plot}{norMmixMLE}(x, y = NULL,
     show.x = TRUE,
     main = sprintf(
         "norMmixMLE(*, model=\"\%s\") fit to n=\%d observations in \%d dim.",
         nm$model, x$nobs, nm$dim
     ),
     sub = paste0(
         sprintf("log likelihood: \%g; npar=\%d", x$logLik, x$npar), 
         if (!is.null(opt <- x$optr)) paste("; optim() counts:", named2char(opt$counts))
     ),
     cex.data = par("cex") / 4, pch.data = 4,
     ...)
 
\S3method{plot}{fittednorMmix}(x, main = "unnamed", plotbest = FALSE, ...)

plot2d (nMm, data = NULL,
        add = FALSE,
        main = NULL,
        sub = NULL,
        type = "l", lty = 2, lwd = if (!is.null(data)) 2 else 1,
        xlim = NULL, ylim = NULL, f.lim = 0.05,
        npoints = 250, lab = FALSE,
        col = Trubetskoy10[1],
        col.data = adjustcolor(par("col"), 1/2),
        cex.data = par("cex"), pch.data = par("pch"),
        fill = TRUE, fillcolor = col, border = NA,
        ...)
plotnd(nMm, data = NULL,
       main = NULL,
       diag.panel = NULL,
       ...)
Trubetskoy10
}
\arguments{
  \item{x, nMm}{an \R object inheriting from \code{"\link{norMmix}"}.}
  \item{y}{further data matrix, first 2 columns will be plotted by
    \code{"\link{points}"}
  }
  \item{\dots}{further arguments to be passed to another plotting function.}
  \item{show.x}{
    Option for `plot.norMmixMLE`. Plot data points along with estimated model.
    Defaults to `TRUE`.
  }
  \item{data}{
    Data points to plot. %% TODO: what format.
  }
  \item{add}{
    This argument is used in the internal function, `plot2d`, to control whether
    to create a new plot or add to an existing one. Should not be set by the user.
    Defaults to `FALSE`
  }
  \item{main}{
    Set main title. See `Usage` section for default values.
  }
  \item{sub}{
    Set subtitle. See `Usage` section for default values.
  }
  \item{type}{
    Graphing type for ellipses border. Defaults to "l".
  }
  \item{lty}{
    Line type to go with the `type`. See \code{"\link{par}"}.
  }
  \item{lwd}{
    Line width as in `lty`.
  }
  \item{xlim}{
    Set explicit x limits for 2d plots.
  }
  \item{ylim}{
    As `xlim`.
  }
  \item{f.lim}{
    Percentage value for how much to extend `xlim` and `ylim`.
    As in the `f` argument to \code{"\link{extendrange}"}.
  }
  \item{npoints}{
    How many points to use in the drawn ellipses. Larger values make them prettier
    but might affect plot times.
  }
  \item{lab}{
    Whether to print labels for mixture components.
    Will print "comp %s" with s in 1 through k.
  }
  \item{col}{
    Fill color for ellipses. Default is "#4363d8".
  }
  \item{col.data}{
    Color to be used for data points.
  }
  \item{cex.data}{
    See \code{"\link{par}"}.
  }
  \item{pch.data}{
    See \code{"\link{par}"}.
  }
  \item{fill}{
    Leave ellipses blank with outline or fill them in.
  }
  \item{fillcolor}{
    Color for infill of ellipses.
  }
  \item{border}{
    %% TODO: I think this should be removed. collides with `type` argument over drawing borders or not. also using lines in any case for border might be sketchy.
  }
  \item{diag.panel}{
    Function to plot 2d projections of a higher-dimensional mixture model.
    Used by `plotnd`.
    Requires function with signature `function(x, y, data = NULL, ...)`
    Should not be set by the user.
  }
  \item{plotbest}{
    Used by `fittednorMmix`. Plot best fitting model using `plot.norMmix`.
  }
}
\value{
  \code{plot.norMmix} In the 2 dimensional case, returns invisibly coordinates 
  of bounding ellipses of distribution.
}
\details{
  The plot method calls one of two auxiliary functions, one for dim=2,
  another for higher dimensions. The method for 2 dimensional plots also takes a
  \code{add} parameter (\code{FALSE} by default), which allows for the ellipses
  to be drawn over an existing plot.

  The higher dimensional plot method relies on the \code{pairs.default} function
  to draw a lattice plot, where the panels are built using the 2 dimensional method.

  \code{Trubetskoy10}: A vector of colors for these plots,
  chosen to be distinguishable and accessible for the colorblind, according to
  \url{https://sashamaps.net/2017/01/11/list-of-20-simple-distinct-colors/},
  slightly rearranged, so that the first five colors stand out well on white background.
}
\examples{% see also  example( norMmix )
plot(MW212) ## and add a finite sample realization:
points(rnorMmix(n=500, MW212))

## or:
x <- points(rnorMmix(n=500, MW212))
plot(MW212, x)


## Example of dim. = p > 2 :
plot(MW34)
}
\keyword{hplot}
