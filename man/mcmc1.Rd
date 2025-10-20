\name{mcmc1}
\alias{mcmc1}
\alias{mcmc0}
\title{Markov Chain Monte Carlo solver}
\description{
  Solve an hmm likelihood using Markov Chain Monte Carlo (MCMC).
}
\usage{
mcmc0(par, fn, iter=0, burnin = 1000, prior.std=4, prop.std=.1,
reset=300, trace=0)
mcmc1(par, fn, iter, burnin, prior.means, prior.std, group, pvar,
      reset= 100, trace = 0, constraint, bfit = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{vector of starting parameters for the iteration}
  \item{fn}{function to compute the log-likelihood}
  \item{iter}{number of MCMC iterations}
  \item{burnin}{number of burnin iterations prior to the MCMC, used to
    tune the proposal parameters}

  \item{prior.means, prior.std}{the mean and std of the Gaussian prior
    distributions that are assumed for each parameter.  Each can be
    single value or a vector of values.}

  \item{group}{optional vector to group the parameters.  If it is
    present this must be a vector of integers.  This causes the MCMC
    updates to be grouped: all of the parameters in group 1, then all
    those in group 2, etc.  The effect is ngroup*iter MCMC steps, of
    which \code{iter} are saved.
  }

  \item{pvar}{a variance matrix for the proposal distribution.  If there
    was a prior fit to the data, for instance, the variance of the
    parameters will lead to a more efficient proposals.}
  \item{prop.std}{standard deviations of the proposal}

  \item{reset}{during the burn in process, how often to recalculate the
    tuning parameters}

  \item{trace}{print a progess line every trace burn in or iteration steps}

  \item{constraint}{optional matrix of constraints}

  \item{bfit}{if TRUE, return results for the burnin process as well as
    iteration results}
}
\details{
The mcmc0 routine is a simple function for testing ideas and will likely
    be removed at a future date.
  
  The \code{par}, \code{fn} and \code{constraint} parameters are set by
the \code{hmm} routine, all the other parameters are the responsibility
of the user and are set via the \code{mpar} argument.
}
\author{Terry Therneau}

\seealso{\code{\link{hmm}}, \code{\link{hmmscore}}}
