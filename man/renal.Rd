\name{renal}
\alias{renal}
\docType{data}
\title{
Renal transplantation data
}
\description{
This is a data set of 407 patients suffering from chronic kidney disease
who underwent a primary renal transplantation with a graft from a deceased
or living donor in the University Hospital of the Catholic University of
Leuven (Belgium) between 21 January 1983 and 16 August 2000.
}
\usage{renal}
\format{
  A data frame with 73091 observations on the following 8 variables.
  \describe{
    \item{\code{id}}{subject identifier, 1 to 416
    \item{\code{day}}{day of observation}
    \item{\code{fail}}{1= graft failure, 0= censor}
    \item{\code{age}}{age at renal transplant}
    \item{\code{gender}}{a factor with levels \code{female} \code{male}}
    \item{\code{prot}}{presence of proteinuria: 0= no, 1= yes}
    \item{\code{haem}}{haemocrit: 100 times the ratio of the
    volume of red blood cells to the total volume of blood}
    \item{\code{gfr}}{glomular filtration rate, ml/min/m^2
  }
}
\details{
This is a data set from the joineRML package, placed into timeline
format. Interesting features are that the GFR and haematocrit trajectories of
individual patients are quite variable in their pattern, time-dependent
GFR and haematocrit are strongly predictive of graft failure, but
baseline values are not.
}
\source{
The joineRML package
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(t4)
## maybe str(t4) ; plot(t4) ...
}
\keyword{datasets}
