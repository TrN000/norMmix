\name{plot.fittednorMmix}
\alias{plot.fittednorMmix}
\title{Plot method for the class \code{fittednorMmix}}
\description{This is the S3 method for plotting the results of \code{\link{fitnMm}}}
\usage{
    \S3method{plot}{fittednorMmix}(x, main="unnamed", plotbest=FALSE, \dots) 
}
\arguments{
    \item{x}{object of class \code{"fittednorMmix"}}
    \item{main}{plot title}
    \item{plotbest}{logical, determines whether to plot BIC values or best fitted model.
                    See Details.}
    \item{\dots}{further arguments to be passed to \code{\link{plot}} if plotbest=TRUE,
                 and \code{\link{matplot}} if FALSE}
}
\details{
    This plot method has two main capabilities, selected by the argument \code{plotbest}.
    If \code{plotbest} is TRUE, then the model will be plotted using the 
    \code{\link{plot.norMmix}} method with added points of the fitted data.
    And if \code{plotbest} is FALSE, then the \code{\link{BIC}} values will be plotted
    using \code{\link{matplot}}
}
\seealso{
    \code{\link{fitnMm}}, \code{\link{norMmix}}, \code{\link{plot.norMmix}}
}
