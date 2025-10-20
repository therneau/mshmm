\name{ptplot}
\alias{ptplot}
\title{
Plot the parse tree of a formula object
}
\description{
The \code{hmm} function uses a list of modified formulas for the
transition rates from state to state, and another such list for
mapping markers of each state to the state. To do this is "walks" the
parse tree of the formula.
This function was very helpful in creating and debugging that part of
the code.
}
\usage{
ptplot(x, depth = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a formula}
  \item{depth}{the maximum depth to plot, 0= no maximum}
  \item{\dots}{further arguments to \code{text} and \code{segments}
    functions that are used to create the plot}
}
\value{
 an invisible list containing the structure of the tree
}
\author{Terry Therneau}
\seealso{\code{\link{hmm}}}
\examples{
\dontrun{
ptplot(A0:A1 + A2:(A3 + death) ~ age + sex / init(.21, .3))
}
}

\keyword{ graphics }
