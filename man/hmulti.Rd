\name{hmulti}
\alias{hmulti}
\alias{mlogit}
\title{create a multinomial response}
\description{This is a helper function, used to create a multinomial
  response function for the \code{hmm} routine.
}
\usage{
hmulti(y, nstate, eta, gradient = FALSE, statemap, weight)
mlogit(eta, gradient=FALSE)
}

\arguments{
  \item{y}{the observed response value, numeric}
  \item{nstate}{the number of states in the model}
  \item{eta}{a matrix with one row per observation and one
    column for each linear predictor}
  \item{gradient}{logical: should the gradient be computed}
  \item{statemap}{a matrix with one row per true state and one
    column for each possible response, that maps outcomes
    to probabilities.}
  \item{weight}{give this the response a weight, e.g., the loglik
    term will be weight*log(f).}  
}
\details{
  Both of these are helper functions, used to construct response
  functions for an \code{hmm} call.
  
  For each observation \eqn{i} there is a vector of predictors
  \eqn{\eta_{ij}}{eta[i,]}.  If there are \eqn{k} predictors
    this defines a set of \eqn{k+1} probabilities:
    \deqn{p_{i1}= 1/\sum_j \eta_{ij}}{p[i,1] = 1/sum(eta[i,])}
    \deqn{p_{i2}= \eta_{i1}/\sum_j \eta_{ij}}{p[i,1] = eta[i,1]/sum(eta[i,])}
    \deqn{p_{i3}= \eta_{i2}/\sum_j \eta_{ij}}{p[i,1] = eta[i,2]/sum(eta[i,])}
    \deqn{\ldots}{\ldots}

    Any (true state, observed y) pair corresponds to one of these
    values. This latter mapping is defined by the \code{statemap}
    matrix, which should have one row per possible true state, one
    column for possible y, and contain the index 1 to k+1 of the
    corresponding probability.

    As an example suppose that we had four true states of
    1= slow/short, 2=slow/tall, 3=fast/short and 4=fast/tall,
    and that y was 1 if the person finished a mile run in under 6
    minutes and 2 otherwise.
    Assume a single column for eta with values <0 (so p[i,1] > p[i,2]).
    A reasonable setup for \code{statemap} would be
    \code{cbind(c(2,2,1,1), c(1,1,2,2))}, which states that if a subject
    is truely in the slow/short group we expect a small probability of
    \code{y==1} and a larger one for \code{y==2}, and etc.
}
\value{for \code{mlogit} the matrix of probabilities \eqn{p}, with or
  without an added gradient attribute.  For \code{hmulti} this has been
  reformated as a response function for \code{hmm}.
}

\seealso{\code{\link{hmm}}}

\examples{
# Create a response function for the example above, which would be
#  used in the \code{rfun} argument of an hmm call.
slowfun <- function(y, nstate, eta, gradient) \{
   map <- cbind(c(2,2,1,1), c(1,1,2,2))
   hmulti(y, nstate, eta, gradient, map)
\}
}
