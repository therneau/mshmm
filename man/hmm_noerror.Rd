\name{hmm_noerror}
\alias{hmm_noerror}
\title{Error function for perfect observation.}
\description{The error function that corresponds to perfect
  observation.  This is used to non-HMM models with the hmm function.}

\usage{
hmm_noerror(y, nstate, ...)
}
\arguments{
  \item{y}{a vector of observed response values, containing integers
    from 1 to the number of true states}
  \item{nstate}{the number of true states}
  \item{\dots}{arguments used by other response functions}
}
\details{
  This function allows the \code{hmm} function to fit non-hidden models.
  It is used mostly in the test suite, to compare results to \code{msm}.
}
\keyword{ models }

