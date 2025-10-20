\name{pade}
\alias{pade}
\title{Pade method for matrix exponential}
\description{
  Compute a matrix exponential, and optionally the derivative,
  using a Pade approximation.
}
\usage{
pade(A, deriv)
}
\arguments{
  \item{A}{a square matrix}
  \item{deriv}{a list of matrices for which derivatives are required}
}
\details{
  When the eigenvalues of \code{A} are unique the \code{hmm} routine
  uses a much faster method based on matrix decompostion.
  This routine is called when this is not applicable.
}
\value{
  the exponential exp(A), and optionally deriviatives of exp(A)
}
\references{
  The expansion is based on a routine in the expm package.  
}
\keyword{ models }

