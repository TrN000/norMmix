\name{plot.norMmix}
\title{Plot Method for "norMmix" Objects}
\alias{plot.norMmix}
\alias{plot.norMmixMLE}
\alias{plot.fittednorMmix}
\alias{nMmcols}
\description{
  This is the S3 method for plotting \code{"\link{norMmix}"} objects.
}
\usage{
\S3method{plot}{norMmix}(x, y=NULL, ...) %%-- TODO all arguments
\S3method{plot}{norMmixMLE}(x, y=NULL, ...)
\S3method{plot}{fittednorMmix}(x, y=NULL, ...)
%% Exportieren und Dokumentieren ?
        %% plot2d(x, ...)
        %% plotnd(x, ...)
nMmcols %% <<< FIXME: *rename* e.g. 'Trubetskoy10' # because we only use 10 of his
}
\arguments{
  \item{x}{an \R object inheriting from \code{"\link{norMmix}"}.}
  \item{y}{further data matrix, first 2 columns will be plotted by \code{"\link{points}"}}
  \item{\dots}{further arguments to be passed to \code{"\link{plot}"}}
}
\value{
  \code{plot.norMmix} returns invisibly coordinates of bounding ellipses of distribution.
}
\details{
  The plot method calls one of two auxiliary functions, one for dim=2,
  another for higher dimensions.% ... FIXME: add more

  \code{nMmcols}: A vector of colors for these plots,
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
