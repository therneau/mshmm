\name{qplus}
\alias{qplus}
\title{Create a Q matrix for a hidden Markov model}
\description{
 A helper function to create the legal transition matrix for an HMM,
  when progression can only be forward
}
\usage{
qplus(data, death = NULL)
}
\arguments{
  \item{data}{A data frame with one row per state, and at least one
    variable \code{states}.  Other variables decribe aspects of the
    state that can "only go forwards".  The order of the states will
    be the order in the Q matrix.}

\item{death}{The label of the death state, if any.  All other states
  are assumped to be able to transition to death. }
}
\details{
This is useful for creating the matrix for a composite state, as in the
example below.
}
\value{A transition matrix with elements of 0 (not transition possible)
  or 1.
}
\seealso{
  \code{\link{cmsh}}
}
\examples{
# An example with states of amyloid burden (A0, A1, A2) and cognitive state
# (C0, C1, C2).
sdata <- data.frame(A= c(rep(0:2, 3), NA), 
                    C= c(rep(0:2, each=3), NA),
                    state= c(outer(c("A0", "A1", "A2"),
                                   c("C0", "C1", "C2"), paste0),
                             "death"))
qmat <- qplus(sdata, death="death")
# qmat
# We see a transtion from A1C0 to A1C1, but not to A0C0 (backing up)
#  or to A1C2 (jumping)
}