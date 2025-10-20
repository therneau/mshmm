\name{hmmemat}
\alias{hmmemat}
\title{Set up an error matrix}
\description{
  Error function for a categorical response.
}
\usage{
hmmemat(y, nstate, eta, gradient = FALSE, statemap, setup, weight)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{a vector of observed outcomes}
  \item{nstate}{number of states in the model}
  \item{eta}{matrix of linear predictors, one row per subject}
  \item{gradient}{should a gradient be returned}
  \item{statemap}{the mapping from overall states to a smaller set,
    redundant to a component of \code{setup}}
  \item{setup}{setup information}
  \item{weight}{give this the response a weight, e.g., the loglik
    term will be weight*log(f).}  
}
\details{
  This is currently packaged into a user written function, called by the
  hmm routine.  The first 4 arguments are supplied by \code{hmm},
  the remaining as defaults set by the user.
  
  An error matrix has a row for each true state and a column for each
  observed outcome.  The matrix need not be square.
(The fact that users often use the same names for
these two  isjoint  concepts can often lead to confusion.)
The set of underlying states is always categorical, this function is
appropriate when the observed outcome is also discrete.

The normal usage is to first call the \code{hmmesetup} routine with a design
matrix; see the documentation for that routine.

Each row of the error matrix produced by \code{hmmemat} corresponds to a
true underlying state, and gives the probability of each outcome if that
state were true.
}
\value{
  an internal structure used by \code{hmm}
}
\seealso{
  \code{\link{hmmesetup}, \link{hmmlcut}}, \code{\link{hmm}}
}
\examples{
# An example with 3 underlying progression states + death, 
#  Those truely in P0 are never misdiagnosed as moderate or high.
mypattern <- rbind(c(-1,  1,  0,  0, 0),
                   c( 2, -1,  3,  0, 0),
                   c( 0,  3, -1,  4, 0),
                   c( 0,  0,  0,  0, -1))
states <- c("P0", "P1", "P2", "death") # generic labels for the hidden state
ylab   <- c("none", "low", "moderate", "high", "death")  # study observations
dimnames(mypattern) <- list(states, ylab)
mset <- hmmesetup(mypattern, states, ylab)
myfun <- function(y, nstate, eta, gradient) {
    hmmemat(y, nstate, eta, gradient, statemap=c(1:4, 0), mset)
}
}
\keyword{ models }
