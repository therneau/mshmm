\name{predict.cmsh}
\alias{predict.cmsh}
\title{Predictions for an cmsh model}
\description{Predictions for an cmsh model}
\usage{
\method{predict}{cmsh}(object, newdata, type = c("link", "rate",
"pstate"), se.fit = FALSE, times, absorb, \dots)
}
\arguments{
  \item{object}{the result of an \code{cmsh} fit}
  \item{newdata}{a data set for prediction}
  \item{type}{the type of predicted value}
  \item{se.fit}{should standard errors of the predicted values be given}
  \item{times}{a vector of times, required for type = pstate}
  \item{absorb}{treat selected states as absorbing}
  \item{\dots}{required of a predict method}
}
\details{
  Predictions of type link are the linear predictors for the model. Rates
  are exp(linear predictor); in this case only rate parameters are
  returned.
  Type prevalence creates a prevalence curve, in which case the new data
  set will need to have the time variable.
}
\value{
  For type prevalence a matrix containing the time points along with the
  prevalence of each state.
  For types link and rate, and with 
  \code{se.fit=FALSE}, a matrix of linear predictions with one row per
  subject and one column per parameter,
  otherwise a list with component fit and se.fit.
}
\author{Terry Therneau}
\seealso{\code{\link{cmsh}}, \code{\link{hmmcontrast}}}
