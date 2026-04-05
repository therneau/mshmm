\name{hmmeigen}
\alias{hmmeigen}
\title{Spectral decomposition of a matrix}
\description{
  Compute the eigenvalues along with the left and right eigenvectors of
  a square matrix
}
\usage{
hmmeigen(x)
}
\arguments{
  \item{x}{a square matrix, normally the transition rate matrix R for a
    multistate hazard model}
}
\details{
This function is a mimic of the R \code{eigen} function, and uses essentially
the same underlying LAPACK call and code.  The only difference is that it
both the left and right eigenvectors, whereas \code{eigen} only returns
the right eigenvectors.  (It is common practice that the term
`eigenvector', unqualified, refers to right eigenvectors).
Both left and right eigenvectors are used by the \code{survexpm}
routine, this routine is exposed mostly as an aid to testing that
function.
}
\value{The spectral decomposition of \code{x} is returned as a list with
  components
  \item{values}{a vector containing the p eigenvalues of \code{x},
    sorted in decreasing order according to
    \code{Mod(values)}. Eigenvalues may be real or complex, complex
    values, if present, will occur as complex conjugate pairs.}
  \item{left}{a square matrix whose rows contain the left eigenvectors,
    each normalized to unit length}
  \item{right}{a square matrix whose columns contain the right
    eigenvectors, each normalized to unit length}
}
\author{Terry Therneau}
\seealso{\code{\link{survexpm}}
}
