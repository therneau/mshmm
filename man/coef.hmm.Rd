\name{coef.hmm}
\alias{coef.hmm}
\title{Extract coefficients from an hmm object}
\description{Return the coefficients of an hmm object, either as
  a vector or a matrix.
}
\usage{
\method{coef}{hmm}(object, type = c("matrix", "raw"), \dots)
}
\arguments{
  \item{object}{an hmm object}
  \item{type}{return a matrix of coefficients, or the raw vector?}
  \item{\dots}{currently ignored}
}
\details{
  The default form is to return a matrix with one row for each
  transition and a column for each coefficient.
  (This is the \code{coef2} component of the object.)

  The raw form is the vector of coefficients as seen by the internal
  maximizer.  Some transitions may share coefficients, and some
  coefficients may be fixed; thus this vector will normally be much
  shorter than the number of non-zero values from the matrix form of
  the print.  See the \code{cmap} component of the object for the
  mapping between transitions and coefficients.
}
\value{a matrix or vector of coefficients}
\seealso{\code{\link{hmm}}, \code{\link{printSpMatrix}}
}
