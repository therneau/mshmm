\name{hmminit}
\alias{hmminit}
\title{Default initial state function.
}
\description{A helper function for hmm, that allows the initial
  state to depend on covariates.
}
\usage{
hmminit(nstate, eta, gradient = FALSE)
}
\arguments{
  \item{nstate}{number of states}
  \item{eta}{vector of linear predictors}
  \item{gradient}{if TRUE return a gradient also}
}
\details{
  Say there were 6 states and \code{eta} was a matrix with two columns.  Then
  the probabilities for each state are assumed to be proportional to
  the vector (1, exp(eta[,1]), exp(eta[,2]), 0, 0, 0).
  These values are then rescaled so as to sum to 1.
  If the matrix \code{eta} has only a single row (a common case) then
  all subjects share the same initial value.

 Since this routine is called by the \code{hmm} function, the details of
the arguments and response are taken care of behind the scenes.
}
\value{
  a matrix of initial values with \code{nrow(eta)} rows and
  \code{nstate} columns.
  Optionally the gradient is also included, which will be an array
  whose dimensionse are \code{(nrow(eta), nstate, ncol(eta))} that
  contains the derivatives with respect to each column of \code{eta}.
}
\seealso{
\code{\link{hmm}}
}
