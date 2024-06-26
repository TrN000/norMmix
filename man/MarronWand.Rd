\name{MarronWand}
\title{Marron-Wand-like Specific Multivariate Normal Mixture 'norMmix' Objects}
%% Note: These are ``datasets'' , but cannot be handled strictly
%%       as such since they *require* the package.
\alias{MarronWand}% to use as "link"
%% --- 2-dim ---
\alias{MW21}
\alias{MW22}
\alias{MW23}
\alias{MW24}
\alias{MW25}
\alias{MW26}
\alias{MW27}
\alias{MW28}
\alias{MW29}
\alias{MW210}
\alias{MW211}
\alias{MW212}
\alias{MW213}
\alias{MW214}
\alias{MW215}
%% 3 d:
\alias{MW31}
\alias{MW32}
\alias{MW33}
\alias{MW34}
%% 5 d:
\alias{MW51}
\description{
  Nicolas Trutmann constructed multivariate versions from most of the
  univariate (i.e., one-dimensional) "Marron-Wand" densities as defined in
  CRAN package \CRANpkg{nor1mix}, see \code{\link[nor1mix]{MarronWand}} (in
  that package).
}
\usage{
## 2-dim examples:
MW21   # Gaussian
MW22   # Skewed
MW23   # Str Skew
MW24   # Kurtotic
MW25   # Outlier
MW26   # Bimodal
MW27   # Separated (bimodal)
MW28   # Asymmetric Bimodal
MW29   # Trimodal
MW210  # Claw
MW211  # Double Claw
MW212  # Asymmetric Claw
MW213  # Asymm. Double Claw
MW214  # Smooth   Comb
MW215  # Trimodal

## 3-dim :
MW31
MW32
MW33
MW34

## 5 - dim:
MW51    # Gaussian
}
\value{
  A normal mixture model. The first digit of the number in the variable name encodes the dimension
  of the mixture; the following digits merely enumerate models, with some correlation to the
  complexity of the model.
}
%% \source{
%% }
\author{Martin Maechler for 1D; Nicolas Trutmann for 2-D, 3-D and 5-D.}
\references{
  Marron, S. and Wand, M. (1992)
  Exact Mean Integrated Squared Error;
  \emph{Annals of Statistcs} \bold{20}, 712--736;
  \doi{10.1214/aos/1176348653}.
}
\examples{
MW210
plot(MW214, main = "plot( MW214 )")% FIXME: show attr(*, "name")  in title

plot(MW51, main = paste("plot( MW51 );  name:", attr(MW51, "name")))
}
\keyword{datasets}
\keyword{distribution}
