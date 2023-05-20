\name{Levi_Tripathy}
\alias{Levi_Tripathy}
\docType{data}
\title{
Visual detection of path deviation
}
\description{
The subject was presented with the image of a dot moving rightwards on a linear path until it reached the midline of the display, when it changed direction either upwards or downwards. The subject had to indicate the direction. The symbols in the figure below show the proportion of correct responses in 30 trials as the deviation varied from –3 to 3 units. See \url{https://personalpages.manchester.ac.uk/staff/david.foster/software-modelfree/latest/examples/example02.html}
}
\usage{data("Levi_Tripathy")}
\format{
  A data frame with 7 rows and 3 columns.
  \describe{
    \item{\code{x}}{stimulus level}
    \item{\code{r}}{number of successes}
    \item{\code{m}}{number of trials}
  }
}
\references{
Levi, D. M. & Tripathy, S. P. “Is the ability to identify deviations in multiple trajectories compromised by amblyopia?”, Journal of Vision, 6(12), 1367-1379, 2006.
}
\examples{
data("Levi_Tripathy")
x = Levi_Tripathy$x
r = Levi_Tripathy$r
m = Levi_Tripathy$m
plot( x, r / m, xlim = c( -2.87, 2.87 ), ylim = c( 0.03, 0.97 ), type = "p", pch="*" )
}
\keyword{datasets}
