\name{survexpm}
\alias{survexpm}
\alias{survexpinit}
\title{
Matrix exponential for transition matrices}
\description{
  Compute the matrix exponential and its derivatives, for the special
  case of a transition matrix
}
\usage{
survexpminit(R)
survexpm(R, time = 1, setup, deriv = 0, all=FALSE, eps = 1e-06)
}
\arguments{
  \item{R}{the transition rate matrix R for a multistate hazard model}
  \item{time}{the result will be exp(time * R)}
  \item{setup}{result of a call to \code{survexpminit}, if missing the
    Pade approach is used}
  \item{deriv}{0: no derivatives, 1= derivative with respect to eta,
    2= derivative with respect to exp(eta)}
  \item{all}{if TRUE, also return intermediate results}
  \item{eps}{criteria for a tied eigenvalue}
}
\details{
  For a Markov process with k states, the k by k transition rate matrix
  R will have off diagonal elements that are 0
  for transitions that cannot occur, >0 for those that can occur, and
  a diagonal such that all row sums are 0. Then P= exp(R), the matrix
  exponential, will be a transition probability matrix:
  all elements of P are >=0 and all rows sum to 1.

  The derivative of P with respect to each element off-diagonal element
  of R can be computed, however, for R[i,j] elements that are zero the
  derivative is 0, and is not returned in the output. Within the icmsh model,
  R[j,k] = exp(eta[j,k]) where each eta is a linear predictor, a
  separate one for each non-zero rate in R. This ensures that rates are
  positive, the same constraint as a Poisson or a Cox model.
  If deriv=1 the routine
  returns derivatives wrt eta, if deriv=2 with respect to exp(eta),
  i.e., the individual non-zero element of R.  Note that changing any
  element of the R matrix can change all the elements of P.

  If the rate matrix is acyclical (no loops) and has unique eigenvalues,
  then an eigenvalue based decompostion method due to Kalbfleisch and
  Lawless is used for
  the matrix exponential and its derivatives, otherwise a Pade
  approximation is used that is parallel to that found in \code{expm:expm}.
  See the code vignette for mathematical details.
  The \code{survexpminit} routine verifies that \code{R} is a
  transition rate matrix, and checks whether it can be reordered into an
  upper triangular form. The result of \code{survexpminit} becomes the
  \code{setup} argument of \code{survexpm}.

  The \code{Matrix::expm} and \code{expm:expm} functions compute the
  matrix exponential in the more general case, and if the fitting functions
  in this package did not require derivatives there would be no point in
  creating our own. The internal routines need to call the setup
  function only once, but compute the matrix exponential hundreds or
  thousands of times (it consumes a majority of the execution
  time), hence the split into two routines.
  The \code{survexpm} function exposes our interal functions, largely
  for testing purposes.

  The matrix exponential has the reputation of being a difficult
  computation due to multiple edge cases, but transition rate matrices
  are one of the more benign arguments to the function.
  A very rough check using a simple 5 state model
  showed relative run times of about 1, 1 and 2.6
  for the eigenvalue based, Pade, and expm functions, resepectively,
  versus 2.5 and 9.2 for the eigenvalue and Pade versions with
  derivatives. See the code in tests/survexpm.R.
  (Though our Pade code is based on that in the expm package, our guess
  is that the speedup is due to the benign nature of rate matrices,
  which allows skipping several checks and special cases.) 
}
\value{
  If \code{deriv=FALSE} and \code{all=FALSE} the matrix exponential
  P= exp(R), otherwise a list with some subset of
  \item{P}{the matrix exponential exp(R)}
  \item{deriv}{an array with dimension(k, k, m) where k is the number of
    states and m the number of non-zero off diagonal elements in R.
    The [,,j] slice of the array contains the derivatives of each
    element of P with respect to the jth non-zero off diagonal.}
  \item{d}{the diagonal of the eigenvalue decomposition
    \eqn{R=ADA^{-1}}{R = A D A-inverse}, D is a diagonal matrix containing
    the eigenvalues}
  \item{A}{the A matrix of the eigenvalue decomposition}
  \item{Ainv}{the inverse of the A matrix}
  \item{nterm}{the number of terms in the Pade approximation}
}
\references{
J Kalbfleisch and J Lawless, The analysis of panel data under a Markov
assumption, J Am Stat Assoc, 1985.

C Moler and C van Loan, Nineteen dubious ways to calculate the
exponential of a matrix, Siam Review, 1978.
}
\author{
Terry Therneau
}
\examples{
rmat <- rbind(c(-.3,  .2,   0,   0,  .1),
              c(  0, -.2,   0, .15, .05),
              c(  0,   0, -.1,   0,  .1),
              c(  0,   0,   .3, -.4, .1), 0)
s1 <- survexpminit(rmat) 
survexpm(rmat, .3, s1)
}
\keyword{array}

