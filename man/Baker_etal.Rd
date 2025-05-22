\name{Baker_etal}
\alias{Baker_etal}
\docType{data}
\title{
Auditory detection of a gap in noise
}
\description{
A 300-ms noise burst containing a gap of 2–8 ms duration or no gap was presented to one ear of a subject. The symbols in the figure below show the proportion of responses "gap" as a function of gap duration. There were 12 trials with each gap duration and 84 trials with no gap. See \url{https://personalpages.manchester.ac.uk/staff/d.h.foster/software-modelfree/latest/examples/Baker_etalR.html}
}
\usage{data("Baker_etal")}
\format{
 A data frame with 8 rows and 3 columns.
  \describe{
    \item{\code{x}}{stimulus level}
    \item{\code{r}}{number of successes}
    \item{\code{m}}{number of trials}
  }
}
\references{
Baker, R. J., Jayewardene, D., Sayle, C., & Saeed, S. “Failure to find asymmetry in auditory gap detection”, Laterality: Asymmetries of Body, Brain and Cognition, 13, 1-21, 2008.
}
\examples{
data("Baker_etal")
x = Schofield_etal$x
r = Schofield_etal$r
m = Schofield_etal$m
plot( x, r / m, xlim = c( 0.16, 7.83 ), ylim = c( -0.01, 1.01 ), type = "p", pch="*" )
}
\keyword{datasets}
