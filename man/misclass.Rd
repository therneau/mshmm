\name{misclass}
\alias{misclass}
\title{Extract categorical marker emission probabilities}
\description{
Extract the estimated or fixed matrix of probabilities for a categorical
marker in a fitted \code{cmsh} model.
}
\usage{
misclass(fit, marker,
         what = c("matrix", "se", "lower", "upper",
                  "coefficients", "all"))
}
\arguments{
  \item{fit}{A fitted \code{cmsh} object.}
  \item{marker}{Marker name or position.  If omitted, the first marker is
  returned; with multiple markers some requests return a named list.}
  \item{what}{The probability matrix, its standard errors or confidence
  limits, the linear-predictor coefficient table, or the complete stored
  marker result.}
}
\details{
Rows correspond to latent marker levels and columns to observed factor levels.
The matrix entries are \eqn{P(\mathrm{observed}\mid\mathrm{latent\ level})}.
}
\value{A matrix, coefficient table, marker-result list, or named list.}
\seealso{\code{\link{cmsh}}, \code{\link{coef.cmsh}}}

