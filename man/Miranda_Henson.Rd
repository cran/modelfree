\name{Miranda_Henson}
\docType{data}
\alias{Miranda_Henson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detection of a flash of light
}
\description{
A flash of light of variable intensity was presented repeatedly at a fixed location in the visual field of a subject who reported whether the flash was visible. There were 3–20 trials at each stimulus level. See \url{https://personalpages.manchester.ac.uk/staff/d.h.foster/software-modelfree/latest/examples/Miranda_HensonR.html}
}
\usage{
data("Miranda_Henson")
}
\format{A data frame with 10 rows and 3 columns.
  \describe{
    \item{\code{x}}{stimulus level}
    \item{\code{r}}{number of successes}
    \item{\code{m}}{number of trials}
  }
}
\references{
  Miranda, M. A. & Henson, D. B. “Perimetric sensitivity and response variability in glaucoma with single-stimulus automated perimetry and multiple-stimulus perimetry with verbal feedback”, Acta Ophthalmologica, 86, 202-206, 2008.
}
\examples{
data("Miranda_Henson")
x = Miranda_Henson$x
r = Miranda_Henson$r
m = Miranda_Henson$m
plot( x, r / m, xlim = c( 0.1, 1.302 ), ylim = c( 0.0165, 0.965 ), type = "p", pch="*" )
}
\keyword{datasets}

