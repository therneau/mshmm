\name{coef.hmm}
\alias{coef.hmm}
\title{Extract coefficients from an hmm object}
\description{Return the coefficients of an hmm object, either as
  a vector or a matrix.
}
\usage{
\method{coef}{hmm}(object, matrix= FALSE, fixed= TRUE, 
                   matrix1= matrix, matrix2= matrix, matrix3= matrix, ...)
}
\arguments{
  \item{object}{an hmm object}
  \item{matrix}{if TRUE, return the coefficients as a matrix, otherwise
    as a vector}
  \item{fixed}{include fixed coefficients in the result}
  \item{matrix1}{if TRUE, include group 1 coefficients in the matrix}
  \item{matrix2}{if TRUE, include group 2 coefficients in the matrix}
  \item{matrix3}{if TRUE, include group 3 coefficients in the matrix}
  \item{\dots}{currently ignored}
}
\details{
  The matrix form will have a row for each covariate used in the model
  and a column for each linear predictor. The structure is parallel to
  the \code{cmap} component; the matrix will have a 0 wherever
  \code{cmap} has a zero, i.e., for a covariate that is not used in the
  given linear predictor. 
  The coefficients can be broken into three groups: 1= those associated
  with the transition probabilities, 2= those associated with markers,
  3 = those associated with initial probability.
  (For some models,  there may be no coefficients in group 2 and/or group 3.)
  To return a matrix containing only coefficients for group 1, for
  instance, one can use \code{matrix1 = TRUE}.
  
  If none of the \code{matrix} arguments are TRUE and
  \code{fixed =FALSE}, the returned vector of coefficients is identical
  to the coefficient vector as seen by the maximization function.
  Some transitions may share coefficients, and some
  coefficients may be fixed; thus this vector will often be much
  shorter than the number of non-zero values from the matrix form.
  For \code{fixed=TRUE} it is the \code{coefficient} component of the
  \code{hmm} object.

  The \code{cmap} component of a fit has one row for each covariate that
  appears and one column for each linear predictor; each transition for
  instance is the exponent of a single linear predictor. If
  \code{cmap[i,j] = 0} then variable i is not part of linear predictor
  j. If \code{cmap[i,j] =k} then linear predictor j has X[i,]*
  coefficient[k] as one of its terms.  If k is negative this indicates
  that the coefficient in question was fixed.
}
\examples{
# assume hfit is the result of an hmm fit
# The following will give the matrix of estimated transition rates
#  for a new data set, one row per observation in new, one column
#  per transition
\dontrun{
 exp(model.matrix(hfit, data=new) %*% coef(hfit, matrix1=TRUE))
}
\value{a matrix or vector of coefficients}
\seealso{\code{\link{hmm}}}
}
