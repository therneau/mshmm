\name{test1}
\alias{test1}
\docType{data}
\title{Small test data set}
\description{
A small subset (4 subjects) of data from an imaging study.
}
\usage{data("test1")}
\format{
  A data frame with 17 observations on the following 12 variables.
  \describe{
    \item{\code{id}}{subject identifier}
    \item{\code{years}}{years from enrollment}
    \item{\code{age}}{age in years}
    \item{\code{male}}{1=male, 0=female}
    \item{\code{educ}}{years of education}
    \item{\code{apoepos}}{positive APOE genotype}
    \item{\code{state}}{current state}
    \item{\code{pib}}{PIB binding (a measure of amyloid)}
    \item{\code{fdg}}{mean FDG image level}
    \item{\code{istate}}{integer version of state variable}
  }
}
\details{
Not useful for an actual fit, but used in several of the tests.
}
\keyword{datasets}
