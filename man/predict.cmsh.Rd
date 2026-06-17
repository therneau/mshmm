\name{predict.cmsh}
\alias{predict.cmsh}
\title{Predictions for an cmsh model}
\description{Predictions for an cmsh model}
\usage{
\method{predict}{cmsh}(object, newdata, type = c("link", "rate",
"pstate"), se.fit = FALSE, time, p0, absorb, \dots)
}
\arguments{
  \item{object}{the result of an \code{cmsh} fit}
  \item{newdata}{a data set for prediction}
  \item{type}{the type of predicted value}
  \item{se.fit}{should standard errors of the predicted values be given}
  \item{time}{a vector, list, or data frame containing time values,
    required for \code{type='pstate'}.  If times is
    a data frame all columns must have the same spacing, and the first
    column is used to label the resulting rows.required for type =
    pstate}
  \item{p0}{starting distribution of states, for pstate curves}
  \item{absorb}{treat selected states as absorbing}
  \item{\dots}{required of a predict method}
}
\details{
  Predictions of type link are the linear predictors for the model. Rates
  are exp(linear predictor); in this case only rate parameters are
  returned.

  For the probability in state curves, the time, newdata and p0 arguments
  are required. Starting with p0, the updates to the estimates
  probability in state are p(t) = p(t-) exp(R(t-)) where R is the
  transition rate matrix, whose elements are exp() of the appropriate
  linear predictor for each state.  The newdata argument will contain
  fixed (not time dependent) predictors, a separate set of curvec will
  be created for each row of newdata.
}
\value{
  For type prevalence an array whose dimensions are time, state, and the
  rows of newdata; if newdata has a single row it will be a matrix.
  The [i,j,k] element contains the estimated probability of being in
  state j at time i, for covariate set k.

  For types link and rate, and with 
  \code{se.fit=FALSE}, a matrix of linear predictions with one row per
  subject and one column per parameter,
  otherwise a list with component fit and se.fit.
}
\author{Terry Therneau}
\seealso{\code{\link{cmsh}}, \code{\link{hmmcontrast}}}
