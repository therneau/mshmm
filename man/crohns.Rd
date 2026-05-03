\name{crohns}
\alias{crohns}
\docType{data}
\title{
A Capitalized Title for the Data Set
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("crohns")}
\format{
  A data frame with 1173 observations on the following 6 variables.
  \describe{
    \item{\code{id}}{a numeric vector}
    \item{\code{time}}{a numeric vector}
    \item{\code{state}}{a factor with levels \code{censor} \code{remission} \code{mild} \code{drug responsive} \code{drug dependent} \code{refractory} \code{surgery} \code{death}}
    \item{\code{dxyr}}{a numeric vector}
    \item{\code{dxage}}{a numeric vector}
    \item{\code{male}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(crohns)
## maybe str(crohns) ; plot(crohns) ...
}
\keyword{datasets}
