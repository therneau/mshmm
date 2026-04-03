\name{hmmpenalty}
\alias{hmmpenalty}
\alias{hmmcontrast}
\title{Create or test contrast matrices in the coefficients}
\description{Both penalties used during an hmm fit and tests of
  coefficient combinations after the fit depend on the creation
  of contrast matrices.
  These two routines for those tasks have a similar form.
}
\usage{
hmmpenalty(fit,  data, lp, contrast)
hmmcontrast(fit, data, lp, contrast)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{the result of an \code{hmm} fit}

  \item{data}{a data set of new subjects}

  \item{lp}{a vector of length \code{nrow(data)}.  It should be
    numeric or match the column names of \code{coef(fit)}.
    Together \code{data} and \code{lp} describe a set of
    \code{nrow(data)} linear predictors.
}
\item{contrast}{a contrast matrix in the linear predictors.
  This can either be a full contrast matrix with k rows
  and \code{nrow(data)} columns, or can be a two column matrix
  containing values from 1 to \code{nrow(data)}.  In the
  latter case all contrasts are pairwise.
}
}
\value{
  For \code{hmmpenalty} a contrast array whose first two dimension
  match those of \code{coef(fit)}, and whose third dimension is the
  number of contrasts.
  For \code{hmmconstrast} a list with elements of \code{estimate},
  a vector giving the value of each linear predictor, \code{var}
  the variance matrix of the linear predictors, and
  \code{test}, a matrix with columns of \code{test}, \code{std},
  \code{p} giving the value of each contast, its estimated standard
  error, and a two sided P value.
}
\seealso{\code{\link{cmsh}}, \code{\link{predict.cmsh}}}
\author{Terry Therneau}
\keyword{ model }% use one of  RShowDoc("KEYWORDS")

