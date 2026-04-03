\name{hmmscore}
\alias{hmmscore}
\title{Fisher scoring method}
\description{
Find a solution to the hmm likelihood using the method of Fisher
scoring,
along with Levenberg-Marquardt (LM) stabilization of the iteration.
}
\usage{
hmmscore(par, fn, gr, iter = 30, scale = 2, shrink = 1, debug = FALSE, constraint, eps = 1e-6, hessian = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{initial parameters}
  \item{fn}{function to evaluate the likelihood}
  \item{gr}{function to evaluate the likelihood and gradient}
  \item{iter}{maximum number of iterations}
  \item{scale}{reduction factor for the LM penalty.  Larger values cause
    the routine to more quickly reduce the penalty.}

  \item{shrink}{step shrinkage.  The actual step taken will be shrink
    times the proposed "best" step from the underlying LM algorithm.
    Shrinkage increases the run time but can make the solution path more
    stable.}
  
  \item{debug}{print extra information during iteration.  A value of 1
    causes some tracing, higher values may invoke the browser function
    at various stages.}

  \item{constraint}{an optional matrix of constraints}

  \item{eps}{convergence criterion for the fit.  Iteration continues
    until the relative change in the log-likelihood is < eps.}
  
  \item{hessian}{if TRUE compute the Hessian matrix at the final
    solution.
    A numeric second derivative is computed by incrementing each
    parameter, separately, by +- delta.
    This calculation requires 2 * (number of parameters)
    evaluations of the derivative, which is often more evaluations than
    the solution.
    The \code{hessian} is a numeric value that value is used for delta,
    otherwise the default value for delta is .001 if
    \code{hessian=TRUE}.
  }
  \item{\dots}{arguments used by other maximization functions}
}
\details{
  The \code{par}, \code{fn} and \code{constraint} arguments are filled
  in by the parent \code{hmm} routine.  Other arguments will be set by the
  user using the \code{mpar} argument of \code{hmm}.
  This list must include \code{gr = "hmmboth"} as one element.
  
  The overall log-likelihood is a sum over the unique subjects
  \eqn{\sum L_i}{sum(L[i])} where
  \eqn{L_i}{L[i]} is the contribuion of subject \eqn{i}.
  Let \eqn{d_{ij}}{d[i,j]} be the derivative of \eqn{L_i}{L[i]} with
  respect to parameter \eqn{j}.
  Then the score matrix \eqn{S} has elements
  \eqn{s_{jk} = \sum_i d_{ij}d_{ik}}{S[i,j] = sum(d[i,j] * d[i,k])};
  at the solution it is an emprical variance of the first derivative,
  which has mean 0 at the solution by definition.
  Fisher scoring uses a Newton-Raphson update, but with S in place of
  the hessian H.  S is much faster to compute.
  
  In Newton-Raphson iteration each trial value is the solution to a
  quadratic equation.  If there are constraints the \code{quadprog}
  library is used to find the best solution under the constraints, which
  is then used as the trial value.  The \code{hmmsore} function is thus
  an interior point method.
}
\value{
\item{loglik}{the log-likelihood at the solution}
\item{deriv}{vector of derivatives at the solution}
\item{S}{the S matrix per above}
\item{S2}{a centered S matrix, i.e. var(D) where D is the matrix of
  derivatives.  This was entertained at one point as an alternative for
  the Fisher scoring.}
\item{coef}{parameter vector at the solution}
\item{iter}{number of iterations used}
\item{converged}{FALSE if the number of iterations was exceeded}
\item{trace}{tracing information. A list with matrices loglik, coef, and
  active, each of which has a row per iteration.}
}
\seealso{\code{\link{cmsh}}}
