\name{hmm.control}
\alias{hmm.control}
\title{Control parameters for the hmm function
}
\description{
  These are miscellaneous control parameters used by the solution
  functions underlying the \code{hmm} function.
}
\usage{
hmm.control(smallpos = 0.001, debug = FALSE, scale=TRUE,
             makecluster = .Platform$OS.type == "windows")
}
\arguments{
  \item{smallpos}{due to round off error the internal computations may
    sometimes create a rate that is less than 0.  When this occurs it is
    replaced by this value.  Usually such conditions are transient, and
    thus the actual value used here is not critical
}
  \item{debug}{print optional debugging information. Larger integers
    lead to more verbose printout.
  }
\item{scale}{internally scale the covariates to have mean 0 and std 1} 
\item{makecluster}{the form of parallel computing that will be used, via
  the \code{parallel} library, either \code{mclapply}, which is based on
  the fork system primitive, or \code{parLappy} which is based on
  clustering. We slihtly prefer the first, but it is not available on
  Windows systems.
}}
\value{a list containing the options}
\seealso{\code{\link{hmm}}
}
\keyword{ models }
