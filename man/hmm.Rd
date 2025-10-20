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
  \item{formula}{a formula for the HMM.  The left hand side will
    normally be a call to hbind or Surv. The first column of y is the
    time of observation, remaining columns are the responses observed at
   that time point.}
  \item{data}{dataframe containing the data to be fit}
  \item{subset}{subset of the data to be used}
  \item{weights}{not currently implemented}
  \item{na.action}{any special actions for missingness}
  \item{id}{the variable name that contains the subject identifier.
    This will normally be in the data set pointed to by \code{data}.}
  \item{istate}{the initial distribution of states.  Only required for
    survival type data.}
  \item{otype}{the observation type. Normally a variable in the
      \code{data} argument.  Values are 0=censored, 1=usual, 2=exact,
      3=entry.  For most data sets most observations will be of type 1.}
  \item{qmatrix}{the square matrix which shows which transitions are
    possible, along with initial estimates of their rate.}
  \item{qcoef}{a data frame which describes which of the terms are
    associated with which transitions, along with contraints on those terms}
  \item{rfun}{a response function or list of response functions, one for
    each response.}
  \item{rcoef}{a data frame which describes which of the terms are associated
    with which response function, along with contraints on those terms}
  \item{pfun}{a function that provides the initial probabilities for
    each state.}
  \item{pcoef}{a data frame which describes which model terms are
    associated with the initial probability function.}
  \item{entry}{a vector giving the probability of enrollment for each
    state, i.e., 0 for the death state.
    This requires data preprocessing, see below.}
  \item{iprob}{optional vector or matrix containing probabilty of each
    state, at entry.
    If a matrix, there is one row for each unique subject, with subject
    as a row name.}
  \item{death}{the state which corresponds to death}
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

  \item{penalty}{an optional array descibing penalties.}

  \item{constraint}{an optional array descibing constraints.}

  \item{cutoff}{the threshold value for constraints}

  \item{debug}{higher values lead to selected debuging printout.}
  \item{makefork}{controls which subportion of the parallel library is
    used: makefork or mcapply.
    The default option of FALSE appears to be faster.}
}
\value{an hmm fit object}
\details{
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
