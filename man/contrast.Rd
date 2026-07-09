\name{contrast}
\alias{contrast}
\title{
Contrasts between selected linear predictors
}
\description{
Given a model and a data set containing the covariates X for a
hypothetical set of subjects, compute the linear predictor for each
subject, and return contrasts between those linear predictors along with
their standard errors.
}
\usage{
contrast(fit, data, global = FALSE, weight, transition, 
         tol = sqrt(.Machine$double.eps))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{a model fit, for any model that has coef and vcov methods,
    e.g., lm, glm, coxph, mshmm, etc.
}
  \item{data}{a data frame containing the covariates of hypothetical
    subjects, one per row
}
  \item{global}{do the global Wald test for \eqn{X\beta=0}{X beta =0},
    where \eqn{\beta}{beta} is the vector of coefficients and \eqn{X}{X}
    the X matrix corresponding to the data rows.
}
  \item{weight}{do a test of each row versus a weighted average of the rows
}
\item{transition}{for a multistate model, restrict the computation to
  selected transitions
}
  \item{tol}{for the global test, a tolerance for the detection of
    singularity in the coefficients, columns deemed to be redundant will
    be skipped
}
}
\details{
  For simple additive models this function is superfluous. If for
  instance the covariates in the model were x1, x2, x3, \ldots, then
  a comparison of two data rows that differ by 1 in x1, but no
  differences in any other, will, the result will be the first
  coefficient and its standard deviation, independent of what values
  were chosen for the others.  That is, it simply reprises a portion of
  the standard printout for the model.

  For a more complex models with for instance non-linear effect or interactions,
  the impact of a given variable cannot be easily read off from the
  standard printout of coefficients and their standard errors, and this
  function is useful. By making the input be in terms of the data,
  issues with respect to how dummy variable were coded become moot.
  
  If neither \code{global} nor \code{weight} are given, then each row is
  contrasted with the first, the first row is the 'reference'.
  For models with interactions, it is often useful to compare each
  group to the overall mean, rather than choosing one combination to be
  a reference.
  The user can decide whether to use weights of 1 for this overall mean
  (average value across the chosen set) or proportional to the
  prevalence of each combination within the larger data set. 
}
\value{a matrix with columns of estimate, std.err, and z}
\examples{
#fit <- glm(is.na(centiloid) ~ age + male * global, binomial, dementia)
#tdata <- expand.grid(age= 70, male=0:1, global=c(-1,0,1))
#cfit <- contrast(fit, tdata, weight=rep(1,6))
}
\keyword{ htest }

