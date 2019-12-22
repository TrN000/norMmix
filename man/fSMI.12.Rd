\name{fSMI.12}
\alias{fSMI.12}
\docType{data}
\title{
    Model selection of the SMI.12 dataset from the \code{\link[copula]{SMI.12}} package.
%%   ~~ data name/kind ... ~~
}
\description{
    Result of \code{fitnMm(SMI.12, k=1:8, ini="clara", maxit=1e4, optREPORT=1e4)}
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("fSMI.12")}
\format{
 The format is:
     List of 7
     $ nMm    :List of 80
      ..- attr(*, "dim")= int [1:2] 8 10
      ..- attr(*, "dimnames")=List of 2
     $ nMmtime:List of 80
      ..- attr(*, "dim")= int [1:2] 8 10
      ..- attr(*, "dimnames")=List of 2
     $ k      : int [1:8] 1 2 3 4 5 6 7 8
     $ models : chr [1:10] "EII" "VII" "EEI" "VEI" ...
     $ n      : int 141
     $ p      : int 20
     $ x      : num [1:141, 1:20] 16.1 15.7 15.7 16.1 16.6 ...
      ..- attr(*, "dimnames")=List of 2
     - attr(*, "class")= chr "fittednorMmix"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(fSMI.12)
## maybe str(fSMI.12) ; plot(fSMI.12) ...
}
\keyword{datasets}
