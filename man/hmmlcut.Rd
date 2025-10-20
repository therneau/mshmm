\name{hmmlcut}
\alias{hmmlcut}
\alias{hmmncut}
\title{Response functions}
\description{
Continuous response functions}
\usage{
hmmncut(y, nstate, eta, gradient = FALSE, cuts, statemap)
hmmlcut(y, nstate, eta, gradient = FALSE, cuts, statemap)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{a vector of observed outcomes}
  \item{nstate}{number of states in the model}
  \item{eta}{matrix of linear predictors, one row per subject}
  \item{gradient}{should a gradient be returned}
  \item{cuts}{cutpoints for the uniform distribution}
  \item{statemap}{the mapping from overall states to this substate}
}
\details{
  These functions are used to build a problem and outcome specific
  response function.  The outcome is assumed to come from a density
  which is the convolution of uniform and gaussian
  (\code{hmmncut}) or the convolution of a logistic and uniform
  (\code{hmmlcut}).

  The statemap vector is used to collapse states.  Say a model has 5 underlying
  states that correspond to amyloid no/yes and neurodegeneration no/yes of:
  A-N-, A+N-, A-N+, A+N+, death.
  If a particular response variable is an (imprecise) measure amyloid burden,
  then statemap would be (1,2,1,2,0): there is one distribution for the
  A- subjects (independent of N), another for A+, and "does not apply"
  for death. 

  For the \code{hmmlcut} and \code{hmmncut}
  the distribution is a blurred uniform.  The two column \code{cuts}
  matrix will contain the range for each, one row per outcome, while
  eta has a single column log(std) of a mean zero convolution.

  For a user, all the arguments except \code{cut} and \code{statemap}
  will simply be copied forward.
}
\value{ an internal structure used by \code{hmm} }
\author{
Terry Therneau
}
\seealso{\code{\link{hmm}}, \code{\link{hmmemat}}
}
\examples{
# This is based on the discussion above
myresp <- function(y, nstate, eta, gradient)
   hmmgauss(y, nstate, eta, gradient, statemap=c(1,2,1,2,0))
}

\keyword{ models }% use one of  RShowDoc("KEYWORDS")
