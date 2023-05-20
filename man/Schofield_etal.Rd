\name{Schofield_etal}
\alias{Schofield_etal}
\docType{data}
\title{
Induction of a visual motion aftereffect
}
\description{
The subject was presented with a moving adaptation stimulus, followed by a test stimulus. The symbols in the figure below show the proportion of responses in which the subject indicated motion of the test stimulus in the same direction as the adapting stimulus, either up or down, as a function of relative modulation depth. There were 10 trials at each stimulus level. See \url{https://personalpages.manchester.ac.uk/staff/david.foster/software-modelfree/latest/examples/example05R.html}
}
\usage{data("Schofield_etal")}
\format{
 A data frame with 7 rows and 3 columns.
  \describe{
    \item{\code{x}}{stimulus level}
    \item{\code{r}}{number of successes}
    \item{\code{m}}{number of trials}
  }
}
\references{
Schofield, A. J., Ledgeway, T., & Hutchinson, C. V. “Asymmetric transfer of the dynamic motion aftereffect between first- and second-order cues and among different second-order cues”, Journal of Vision, 7(8), 1-12, 2007.
}
\examples{
data("Schofield_etal")
x = Schofield_etal$x
r = Schofield_etal$r
m = Schofield_etal$m
plot( x, r / m, xlim = c( 2, 98.2 ), ylim = c( 0.02, 0.98 ), type = "p", pch="*" )
}
\keyword{datasets}
