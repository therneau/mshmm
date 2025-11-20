\name{hmm}
\alias{hmm}
\title{Hidden Markov Models}

\description{Fit an hmm model to a population data set}

\usage{
hmm(formula, data, subset, weights, na.action=na.pass, id,
    istate, otype, qmatrix, qcoef,
    rfun, rcoef, pfun=hmminit, pcoef, entry, iprob, death,
    mfun=hmmtest,  mpar=list(), mc.cores=getOption("mc.cores", 2L),
    icoef, scale=c(FALSE,FALSE), penalty, constraint, cutoff=0,
    debug=0, makefork=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula for the HMM.  The left hand side will be a
    numeric vector or a call to Surv.
  \item{data}{dataframe containing the data to be fit}
  \item{subset}{subset of the data to be used}
  \item{weights}{not currently implemented}
  \item{na.action}{any special actions for missingness}
  \item{id}{the variable name that contains the subject identifier.
    This will normally be in the data set pointed to by \code{data}.}
  \item{qmatrix}{a square matrix whose row and column labels are the
    states, positive elements indicate valid transitions}
  \item{markers}{a list of formulas defining state:marker relationships}
  \item{iprob}{optional, the initial distribution of states. Either a
    matrix with one row per unique id and one column per state, or a
    vector with one element per state.}
  \item{pfun}{a function that provides the initial probabilities for
    each state.}
  \item{pcoef}{a data frame which describes which model terms are
    associated with the initial probability function.}
  \item{entry}{a vector giving the probability of enrollment for each
    state, i.e., 0 for the death state.
    This requires data preprocessing, see below.}
  \item{death}{optional, the state which corresponds to death}
  \item{exact}{optional, states for which the time of entry to the state
  is known exactly}
  \item{mfun}{the optimization function to use}
  \item{mpar}{any parameters for the optimization function}
  \item{mc.cores}{how many parallel cores to use for computing the
    per-subject log-likelihood contributions}
  \item{icoef}{initial coefficients.  An alternative to inital values in
    the \code{qcoef}, \code{rcoef} and \code{pcoef} arguments.  (Useful
    for restarting a model.)}

  \item{scale}{two logical values where the first determines if the
    covariates are recentered and the second if they should also be
    scaled.  A single logical is replicated.
    Alternatively, a list containing the centering and scaling vectors.}

  \item{penalty}{an optional data frame descibing penalties.}

  \item{constraint}{an optional data frame descibing constraints.}

  \item{cutoff}{the threshold value for constraints}

  \item{debug}{higher values lead to selected debuging printout.}
  \item{makefork}{controls which portion of the parallel library is
    used: the fork option is not available in MS-Windows.}
}
\value{an hmm fit object}
\details{
  If the response is \code{Surv(time, status)}, then the status variable
  should be a factor whose first level corresponds to censored, and
  other levels are a subset of the states in \code{qmatrix}.  Any states
  in said subset are assumed to me measured without error. Normally, the
  only state of that type, in a hidden Markov model (HMM), will be death,
  and in that particular case a 0/1 status code is allowed.
  
  It is possible to have exact states other than death, one example
  might be a particular surgery.  The difference arises with censored
  observations: at that point we know the subject is not dead, even for
  an out of study follow-up, but we would not know the status of other
  states.

  
  Entry into a trial is one important issue for HMM models, we allow three
  different likelihoods: a common entry probability vector for all, per-
  subject entry probabilities, and computed probabilities.  The first
  two can depend on parameters and covariates, and depend on the
  \code{iprob} and \code{rcoef} arguments.

The last is applicable when subjects were sampled at random from a
population, and makes use of the \code{entry} argument.  Assume that
subjects are entered at various ages, but we are willling to assume
that all come from a single distribution at younger ages, e.g., no one
is demented before age 50.
Then the data can be set up so that everyone ``enters'' at age 50 in the
known starting state, and the actual entry observation marked as \code{obtype=3}.
At that point the current state probality estimates will be rescaled
using the \code{entry} values, e.g., no one enrolls after death.
The necessary preprocessing, i.e., adding the extra rows to the data
set, is not feasable if there are time-dependent covariates for the hazard.

The \code{icoef} argument can be a matrix in the form of the returned
\code{beta} value, or a vector as in the returned \code{coefficient}
vector.  It is useful for restarting a model at a prior solution.
If the former, it is matched to the current model's \code{beta} matrix
by row and column names, replacing those rows/columns.
This allows one to use coefficients from a prior model that was smaller.
}
\author{Wiste, Weigand, and Therneau}
