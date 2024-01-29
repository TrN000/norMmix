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
                                if (!is.null(opt <- x$optr)) { paste("; optim() counts:", named2char(opt$counts)) }
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
        col.data = adjustcolor(par("col"), 1 / 2),
        cex.data = par("cex"), pch.data = par("pch"),
        fill = TRUE, fillcolor = col, border = NA,
        ...)
plotnd(nMm, data = NULL,
       main = NULL,
       diag.panel = NULL,
       ...)
Trubetskoy10 %% <<< FIXME: *rename* e.g. 'Trubetskoy10' # because we only use 10 of his
}
\arguments{
  \item{x, nMm}{an \R object inheriting from \code{"\link{norMmix}"}.}
  \item{y}{further data matrix, first 2 columns will be plotted by \code{"\link{points}"}}
  \item{\dots}{further arguments to be passed to another plotting function.}
  \item{show.x}{
  }
  \item{data}{
  }
  \item{add}{
  }
  \item{main}{
  }
  \item{sub}{
  }
  \item{type}{
  }
  \item{lty}{
  }
  \item{lwd}{
  }
  \item{xlim}{
  }
  \item{ylim}{
  }
  \item{f.lim}{
  }
  \item{npoints}{
  }
  \item{lab}{
  }
  \item{col}{
  }
  \item{col.data}{
  }
  \item{cex.data}{
  }
  \item{pch.data}{
  }
  \item{fill}{
  }
  \item{fillcolor}{
  }
  \item{border}{
  }
  \item{diag.panel}{
  }
  \item{plotbest}{
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
