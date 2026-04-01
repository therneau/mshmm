\name{icmsh}
\alias{icmsh}
\title{Interval censored multi-state hazard models}

\description{Fit a multi-state hazard model to data that may be interval
censored.  When all states are only paritially observed (latent or
hidden) this becomes a hidden Markov model (HMM).  This function allows
for mix of hidden and directly observed states.}

\usage{
icmsh(formula, data, subset, weights, 
    id, qmatrix, markers, iprob, init, fixed,  
    penalty,  constraint, statedata,
    iter=30, exact= "death", mfun=hmmscore, mpar=list(), mfattr, 
    mc.cores=getOption("mc.cores", 2L), control= icmsh.control(), ...)
}
\arguments{
  \item{formula}{a formula for the model.  The left hand side will be a
    numeric vector or a call to Surv.}
  \item{data}{dataframe containing the data to be fit}
  \item{subset}{subset of the data to be used}
  \item{weights}{not currently implemented}
  \item{id}{the variable name that contains the subject identifier.
    This will normally be in the data set pointed to by \code{data}.}
  \item{qmatrix}{a square matrix whose row and column labels are the
    states, positive elements indicate valid transitions}
  \item{markers}{an optional list of formulas defining state:marker
      relationships, a key component of hidden Markov models (HMM)}
  \item{iprob}{optional, the initial distribution of states. Either a
    matrix with one row per unique id and one column per state, a
    vector with one element per state, or a formula}
  \item{init}{optional, initial values for the parameters}
  \item{fixed}{optional, initial values that are fixed}
  \item{penalty}{an optional data frame descibing penalties.}
  \item{constraint}{an optional data frame descibing constraints.}
  \item{statedata}{a data frame which provides aliases for the states}
  \item{iter}{maximum number of iterations}
  \item{exact}{states for which the time of entry to the state
    is known exactly}
  \item{mfun}{the optimization function to use}
  \item{mpar}{any parameters for the optimization function}
  \item{mfattr}{optional information about the optimization function}
  \item{mc.cores}{how many parallel cores to use for computing the
    per-subject log-likelihood contributions}
  \item{control}{a list of control parameters, normally created with
      the \code{icmsh.control} function.}
  \item{...}{optional parameters, any names that match control arguments are
      used as input to \code{icmsh.control}}
}
\value{an icmsh fit object}
\details{
If the response is \code{Surv(time, status)}, then the status variable
should be a factor whose first level corresponds to censored, and
other levels are a subset of the states in \code{qmatrix}.  Any observed
outcomes are assumed to me measured without error, which would be the
case for simple interval censored data.  For an HMM, an observation will
normally be marked as censored and the true state is (imperfectly) inferred
from marker variables.
In the particular case that the only completely observed state is death,
then a 0/1 status code is allowed, with the assumption that 1=death.
  
It is possible to have exact states other than death, one example
might be a particular surgery.  Formally, a censored observation is only
known to not be any state that is both exact and absorbing, of which death is
the most common example.
(Absorbing states are inferred from \code{qmatrix}).

The \code{icoef} argument can be a prior model fit, often a smaller
model with fewer parameters.
The argument may also be a matrix with transitions as the columns and
predictors as the rows, i.e., in the same form as the \code{cmap}
component or as provided by \code{coef(fit, matrix=TRUE)}.
Last it can be a named vector with names similar to \code{coef(fit)}.
In each case not all transitions or 
covariates need be present, the remainder are assumed to be 0. 
This makes it easy to use a smaller fit as intial values for a larger one.

Missing values are dealt with in a way peculiar to the icmsh function.
Any \code{id} value with a missing value for the response, weight, or intial
rate covariates is removed, i.e., all rows for that id. 
This is because retaining only some rows
is likely to lead to an invalid timeline for that subject, e.g., one with gaps.
Rate covariates are
those that appear on the right hand side of the primary formula.
Missing rate variables at intermediate times are filled in using last value
carried forward. 
Missing markers are ignored, as episodic observation of marker variables is
expected.
}
\author{Terry Therneau}