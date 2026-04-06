\name{cmsh.control}
\alias{cmsh.control}
\title{Control parameters for the cmsh function
}
\description{
  These are miscellaneous control parameters used by the solution
  functions underlying the \code{cmsh} function.
}
\usage{
cmsh.control(smallpos = 0.001, debug = FALSE, center=TRUE, scale=FALSE,
             makecluster = .Platform$OS.type == "windows", detail=FALSE)
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
  \item{center}{internally center covariates}
  \item{scale}{internally scale the covariates to have std 1} 
\item{makecluster}{the form of parallel computing that will be used, via
  the \code{parallel} library, either \code{mclapply}, which is based on
  the fork system primitive, or \code{parLappy} which is based on
  clustering. We slihtly prefer the first, but it is not available on
  Windows systems.
}
\item{detail}{If TRUE, no iteration is done, and the contributions of
  each individual id to the estimate and derivatives is returned}
}
\value{a list containing the options}
\seealso{\code{\link{cmsh}}
}
\keyword{ models }
