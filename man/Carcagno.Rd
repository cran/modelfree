\name{Carcagno}
\alias{Carcagno}
\docType{data}
\title{
Discrimination of pitch
}
\description{
The subject had to identify the interval containing a tone whose fundamental frequency was different from that in the other two intervals. The symbols in the figure below show the proportion of correct responses as the difference between the tones varied. There were 3–49 trials at each stimulus level. See https://personalpages.manchester.ac.uk/staff/david.foster/software-modelfree/latest/examples/example03.html
}
\usage{data("Carcagno")}
\format{
  A data frame with 8 rows and 3 columns.
  \describe{
    \item{\code{x}}{stimulus level}
    \item{\code{r}}{number of successes}
    \item{\code{m}}{number of trials}
  }
}
\references{
Unpublished data from S. Carcagno, Lancaster University, July 2008.
}
\examples{
data("Carcagno")
x = Carcagno$x
r = Carcagno$r
m = Carcagno$m
plot( x, r / m, xlim = c( 1.95, 4.35 ), ylim = c( 0.24, 0.99 ), type = "p", pch="*" )
}
\keyword{datasets}

