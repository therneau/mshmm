\name{survexpm}
\alias{survexpm}
\title{
Matrix exponential for transition matrices}
\description{
  Compute the matrix exponential and its derivatives
}
\usage{
survexpm(R, time = 1, deriv = FALSE, tol=1e-8, method=c("eigen", "pade"))
}
\arguments{
  \item{R}{the transition rate matrix R for a multistate hazard model}
  \item{time}{the result will be exp(time * R)}
  \item{deriv}{compute derivatives}
  \item{tol}{the tolerance for singularity in the eigenvalues}
  \item{method}{use the eigenvalue approach (if possible) or a Pade
  approximation}
}
\details{
  For a Markov process with k states, the k by k transition rate matrix
  R will have off diagonal elements that are 0
  for transitions that cannot occur, >0 for those that can occur, and
  a diagonal such that all row sums are 0. Then P= exp(R*time), the matrix
  exponential, will be a transition probability matrix:
  all elements of P are >=0 and all rows sum to 1. It estimates the
  transition probabilities over a time interval \code{time}.

  The derivative of P with respect to each non-zero off-diagonal element
  is computed (for R[i,j] elements that are zero the derivative is 0).
  Changing any element of the R matrix can change all the elements of
  P= exp(R). The derivatives are thus a set of matrices, each of the
  same shape as R, one per non-zero off diagonal element of R.

  If the eigenvectors corresponding to the non-zero eigenvalues are of
  full rank, then an eigenvalue based decompostion method is used for
  the matrix exponential and its derivatives, see Kalbfleisch and Lawless.
  Otherwise a Pade approximation is used that is parallel to the Ward
  method found in \code{expm:expm}.
  See the code vignette for mathematical details.

  The \code{Matrix::expm} and \code{expm:expm} functions compute the
  matrix exponential in the more general case, and if the fitting functions
  in this package did not require derivatives there would have been
  little point in creating our own.  The eigenvectors and eigenvalues are
  computed using the lapack dgeev routine, the same function used by
  the R \code{eigen} function.

  The matrix exponential has the reputation of being a difficult
  computation due to multiple edge cases, but transition rate matrices
  are one of the more benign arguments to the function.  Our Pade code
  skips preconditioning steps found in expm which are important for more
  difficult cases, but which would also make derivatives much harder.
  A very rough check using a simple 5 state model with 7 transitions
  showed relative run times of approximately 1, 1.4 and 2.6
  for the eigenvalue based, Pade, and expm functions, resepectively;
  2.9 and 9.7 for the eigenvalue and Pade versions with derivatives. 
  A more complex rate matrix with 17 states and 40 transitions
  increased these relative times to 2.2, 2.9, 5.6, 29.1, and 148,
  respectfully.
}
\value{
  If \code{deriv=FALSE} the matrix exponential P= exp(R),
  otherwise a list with elements
  \item{P}{the matrix exponential exp(R)}
  \item{dmat}{an array with dimension(k, k, m) where k is the number of
    states and m the number of non-zero off diagonal elements in R.
    The [,,j] slice of the array contains the derivatives of each
    element of P with respect to the jth non-zero off diagonal.}
  \item{method}{the method used for the computation}
}
\references{
J Kalbfleisch and J Lawless, The analysis of panel data under a Markov
assumption, J Am Stat Assoc, 1985.

C Moler and C van Loan, Nineteen dubious ways to calculate the
exponential of a matrix, twenty-five years later. Siam Review, 2003.
}
\author{
Terry Therneau
}
\examples{
rmat <- rbind(c(-.3,  .2,   0,   0,  .1),
              c(  0, -.2,   0, .15, .05),
              c(  0,   0, -.1,   0,  .1),
              c(  0,   0,   .3, -.4, .1), 0)
survexpm(rmat, .3)
}
\keyword{array}

