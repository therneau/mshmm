\name{hmmtest}
\alias{hmmtest}
\title{Dummy maximization function}
\description{
  This is used as the maximization function in \code{hmm} to do a fit with
  no iterations, i.e., return the log-likelihood at the intial values.
}
\usage{
hmmtest(par, fn, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{the parameters}
  \item{fn}{internal maximization function}
  \item{\dots}{other arguments}
}
\details{
  When used as the \code{mfun} argument of \code{hmm}, it results in
  a fitted value and derivatives at the intial parameters.
  The function is useful in the test suite, or to check that initial
  values have been set correctly.
}
\seealso{\code{\link{hmm}}}
