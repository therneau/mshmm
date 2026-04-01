\name{icmsh.object}
\alias{icmsh.object}
\title{Multi-state model regression object}
\description{This class of objects is returned by the \code{icmsh}
  function.
  Objects of this class have methods for \code{print}, \code{predict}
  and \code{coef}.
}

\section{Components}{
The following components should be included in a legitimate \code{icmsh}
object.
}
\arguments{
  \item{coefficients}{the vector fitted coefficients. }
  \item{loglik}{a vector of length 2 giving the log-likelihood at the
    initial and final iteration.}
  \item{beta}{the coefficients of the model arranged as a matrix.  This
    includes both fixed and estimated values.}
  \item{n}{number of observations in the data}
  \item{qmatrix}{The matrix of transtions between states. The [i,j]
    value is positive if
    transitions are possible directly from state i to state j.
    The row names of the matrix provides the state names.}
  \item{cmap}{the mapping matrix from states to coefficients.}
  \item{scale}{used to internally scale the X matrix.}
  \item{fit}{the object returned by the underlying fitting routine}
  \item{eval}{the number of likelihood evaluations, and the number of
    these that could not use the fast matrix factorization.  A single
    iteration will have one evaluation per subject.}
  \item{time}{the amount of compute time used by various parts of the
    code}
  \item{nlp}{the number of linear predictor columns for the rates matrix (q),
    response functions (m) and initial state (p). This is also marks
    which columns of \code{cmap} correspond to each portion.}
  \item{na.action}{the na.action attribute, if any, that was returned by
    the \code{\link{model.frame}} routine, containing information on
    any deleted observations.}
%% It's OK to have text outside the '\item' entries
The object will also contain the following, for documentation see
the \code{lm} object: \code{terms}, \code{assign}, \code{formula},
\code{call}, and, optionally, \code{x}, \code{y}, 
and/or \code{frame}.
}


\seealso{\code{\link{icmsh}}}

  
