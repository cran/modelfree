\name{Xie_Griffin}
\alias{Xie_Griffin}

\docType{data}
\title{
Discrimination of "porthole" views of natural scenes
}
\description{
The subject was presented with a display split into two parts, one containing a pair of patches from the same image, the other a pair from different images, and the subject had to judge which pair came from the same image. The symbols in the figure below show the proportion of correct responses in 200 trials as a function of patch separation. See https://personalpages.manchester.ac.uk/staff/david.foster/software-modelfree/latest/examples/example04R.html
}
\usage{data("Xie_Griffin")}
\format{
  A data frame with 10 rows and 3 columns.
  \describe{
    \item{\code{x}}{stimulus level}
    \item{\code{r}}{number of successes}
    \item{\code{m}}{number of trials}
  }
}
\references{
Xie, Y. & Griffin, L. D. “A 'portholes' experiment for probing perception of small patches of natural images”, Perception, 36, 315, 2007.
}
\examples{
data("Xie_Griffin")
x = Xie_Griffin$x
r = Xie_Griffin$r
m = Xie_Griffin$m
plot( x, r / m, xlim = c( 0.25, 8.76 ), ylim = c( 0.52, 0.99 ), type = "p", pch="*" )
}
\keyword{datasets}

