\name{EII}
\name{VII}
\name{EEI}
\name{VEI}
\name{EVI}
\name{VVI}
\name{EEE}
\name{VEE}
\name{EVV}
\name{VVV}
\alias{EII}
\alias{VII}
\alias{EEI}
\alias{VEI}
\alias{EVI}
\alias{VVI}
\alias{EEE}
\alias{VEE}
\alias{EVV}
\alias{VVV}
\docType{data}
\title{
  Normal Mixture Model Types
}
\description{
  This package uses the Cholesky decomposition to parametrize the
  covariance matrices of normal mixture models.

  This object lists them roughly, but not precisely, in increasing order
  of complexity. When a model is or can be refered to by an integer it
  is to be understood as an index into this list. See details section.
}
\usage{data("EII")}
\usage{data("VII")}
\usage{data("EEI")}
\usage{data("VEI")}
\usage{data("EVI")}
\usage{data("VVI")}
\usage{data("EEE")}
\usage{data("VEE")}
\usage{data("EVV")}
\usage{data("VVV")}
\format{
  The format is:
 chr [1:10] "EII" "VII" "EEI" "VEI" "EVI" "VVI" "EEE" "VEE" "EVV" "VVV"
}

\details{
  A symmetric positive definite matrix can be decomposed into a constant
  factor, a lower triagonal, and a diagonal matrix; known as the
  Cholesky decomposition.
  The model names correspond to how those components may vary.

  The first component is a multiplicative factor, corresponding to the
  determinant of the matrix. Across multiple matrices they may either be
  'E'qual or 'V'arying.

  The second component describes the diagonal components. These are one
  of the following: the 'I'dentity diagonal matrix, 'E'qual across
  matrices, or 'V'arying.

  The third corresponds to the lower triagonal matrix. Which also can be
  the 'I'dentity matrix, making the whole matrix diagonal, 'Equal', or
  'V'arying.
}
\keyword{datasets}
