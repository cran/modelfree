\name{Nascimento_etal}
\alias{Nascimento_etal}
\docType{data}
\title{
Discrimination of image approximations
}
\description{
The subject was shown an image of a natural scene and an approximation of this image based on principal component analysis. The task was to distinguish between the images. The symbols in the figure below show the proportion of correct responses as a function of number of components in the approximation. There were 200 trials at each level pooled over a range of natural scenes. See \url{https://personalpages.manchester.ac.uk/staff/d.h.foster/software-modelfree/latest/examples/Nascimento_etalR.html}
}
\usage{data("Nascimento_etal")}
\format{
A data frame with 8 rows and 3 columns.
  \describe{
    \item{\code{x}}{stimulus level}
    \item{\code{r}}{number of successes}
    \item{\code{m}}{number of trials}
  }
}
\references{
Nascimento, S.M.C., Foster, D.H., & Amano, K. “Psychophysical estimates of the number of spectral-reflectance basis functions needed to reproduce natural scenes”, Journal of the Optical Society of America A-Optics Image Science and Vision, 22 (6), 1017-1022, 2005.
}
\examples{
data("Nascimento_etal")
x = Schofield_etal$x
r = Schofield_etal$r
m = Schofield_etal$m
plot( x, r / m, xlim = c( 1.17, 7.8 ), ylim = c( 0.47, 1.03 ), type = "p", pch="*" )
}
\keyword{datasets}
