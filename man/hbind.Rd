\name{hbind}
\alias{hbind}
\title{create a response matrix}

\description{The left hand side of the model equation for an hmm model
  will have a time variable followed by multiple responses.
  The \code{hbind} function binds these together as a numberic
  matrix, properly handling characters and factors.
}
\usage{
hbind(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{a time variable followed by one or more response variables}
}
\details{
 This is essentiall a variant of cbind which deals with factors in a
 more proactive manner.  If all the arguments are numeric the result is
 no different than cbind.
}
\value{a matrix, with possible attributes for the factor levels}
\seealso{\code{\link{hmm}}}
