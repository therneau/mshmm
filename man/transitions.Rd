\name{transitions}
\alias{transitions}
\title{
Count state transtions
}
\description{
This creates a table of transitions between states, a subject with k
rows of data will contribute k-1 transitions to the counts.
Row for a subject are summarized in time order.
}
\usage{
transitions(id, time, state)
}
\arguments{
  \item{id}{the subject identifier}
  \item{time}{time of the observation}
  \item{state}{state at that time point}
}
\value{a simple tabulation, NA is included}
}
\keyword{ utilities }
