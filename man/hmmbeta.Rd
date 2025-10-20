\name{hmmbeta}
\alias{hmmbeta}
\alias{hmmexp}
\alias{hmmgamma}
\alias{hmmgauss}
\title{Continuous response functions}
\description{
These are continuous response functions for the hidden Markov model
(hmm) function.
}
\usage{
hmmbeta(y, nstate, eta, gradient = FALSE, statemap, weight=1)
hmmexp(y, nstate, eta, gradient = FALSE, statemap, weight =1)
hmmgamma(y, nstate, eta, gradient = FALSE, statemap, weight =1)
hmmgauss(y, nstate, eta, gradient = FALSE, statemap, weight =1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{the vector of observed responses
  }
  \item{nstate}{number of states
  }
  \item{eta}{matrix of linear predictors
  }
  \item{gradient}{shoud a gradient also be computed?
  }
  \item{statemap}{the mapping between true states and y
  }
  \item{weight}{give this the response a weight, e.g., the loglik
    term will be weight*log(f).}  
}
\details{
  These are used as building blocks by users, to create response functions
  for an \code{hmm} call.
  In this case the first four arguments will be passed through; the
  \code{statemap} argument and optionally \code{weight}
  will be supplied by the user.  As an example
  assume that we had 5 states of A-/N-, A+/N-, A-/N+, A+/N+ and death in
  that order, and that the outcome in question is a particular
  laboratory test which is related to A but not to N.
  The outcome is not measured at death.  Then \code{statemap} would be
  (1, 2, 1, 2, 0): it maps the 5 model states to the two relevant groups
  for this outcome; death corresponds to neither group.
  The \code{statemap} vector must be an integer vector with length equal
  to the number of states, which maps them into ngroup distinct sets.
  
  For the beta distribution \code{eta} will have 2*ngroup columns,
  exp(eta) will contain the 2 shape parameters for group 1, then those
  for group 2, etc.  The gamma distribution is similar, with but with
  shape and rate parameters for each group.

  For the exponential distribution \code{eta} with have ngroup columns
  with exp(eta) the rate for each distribution.

  For the gaussian distribution eta can have either ngroup +1 or
  2*ngroup columns.  In the second case columns 1 to ngroup are the means
  followed by ngroup colums containing the log(standard deviations).
  If there ngroup+1 columns then a common
  standard deviation is assumed whose logarithm appears in the last
  column.

  If two or more responses map to the same state(s), then the weight
  argument can be used to a priori give one of them less influence.
  This is experimental.
}
\seealso{\code{\link{hmm}}
}
\examples{
# a gaussian response matching the example given above
mygauss <- function(y, nstate, eta, gradient) {
    hmmgauss(y, nstate, eta, gradient, statemap=c(1,2,1,2,0))
}
}