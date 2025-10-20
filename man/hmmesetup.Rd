\name{hmmesetup}
\alias{hmmesetup}
\title{Set up an error matrix}
\description{
  Accept user input for an error matrix, check, and  produce a
  summary that will be used by \code{hmmemat}
}
\usage{
hmmesetup(emat, states, response, statemap, partial=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{emat}{user defined matrix for the pattern of the errors}
  \item{states}{vector containing the names of the states, optional}
  \item{response}{vector containing the levels of the response variable,
  optional}
  \item{statemap}{the mapping from overall states to a smaller set}
  \item{partial}{partially observed states, if any}
}
\details{
  For an hmm, there is a mapping between the true states and each
  biomarker of those states, this function applies for discrete outcome
  markers.

  The first argument is a pattern matrix with a row for each true state,
  and a column for each possible observed result $y$ for a given marker.
  The is a row for each true state, and a column for each level of the
  marker (which will be a factor).
  Elements of the matrix are p[i,j] =0: =  outcome j will not occur when
  the true state is i; p[i,j] = -1: outcome j is the reference cell for
  this state; p[i,j] = 1, 2, 3, ... the probability of this outcome
  is based on linear predictor 1, 2, 3, etc.

  Arguments 2 and 3 are vectors of the true states the levels of the
  outcome, respectively.  If present the dimnames of \code{emat} are
  checked for consistency with these.  If missing, the dimnames of
  \code{emat} are assumed to be correct as is.  (Older versions of the
  routine lacked these checks.)

  The optional statemap argument is useful if some patterns repeat,
  e.g., if \code{statemap <- c(1,1,1,2,2,3)} there are 6 states, the
  error matrix is identical for the first three, and corresponds to
  the first row of \code{emat}.

  The last argument is for partial observed states, also referred to
  as \emph{censored}.  Say that the levels of the response were 
  or D, but there are some visits where it is only known that
  y= A or B, you know the subject is in category A or B but not which
  one, and that is coded as factor level "A/B".  Then the
  option \code{partial = list("A/B" = c("A", "B"))} codes for this. This
  is also referred to as a censored outcome, in parallel to left
  censoring for a continuous variable, e.g., a laboratory assay with a
  lower threshold of detection.
}   
\value{
  a list, used as input to the \code{hmmemat} function}
\seealso{
  \code{\link{hmmemat}},\code{\link{hmmlcut}}, \code{\link{hmm}}
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
