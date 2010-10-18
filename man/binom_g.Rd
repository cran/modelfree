\name{binom_g}
\alias{binom_g}
\title{Psychometric function with guessing rate}
\description{THIS IS AN INTERNAL FUNCTION: USE \code{BINOM_LIMS} FOR BEST RESULTS. Maximum likelihood estimates of the parameters of psychometric function with guessing rate (GLM). The estimated parameters for the linear part are in vector b and the estimated guessing rate is guess (GLM).}

\usage{
binom_g( r, m, x, link, p, K, initval )
}

\arguments{
\item{r}{number of successes in points \code{x}}
\item{m}{number of trials in points \code{x}}
\item{x}{stimulus levels}
\item{link}{link function}
\item{p}{degree of the polynomial to be fitted on the linear scale}
\item{K}{Power parameter in Weibull and reverse Weibull models}
\item{initval}{initial value for guessing rate}
}
\value{
\item{value}{Object with 3 components: 
guessing: estimated guessing rate
b: vector of estimated coefficients for the linear part
fit: glm object to be used in evaluation of fitted values}
}

\examples{
data( "01_Miranda" )
value <-binom_g( example01$r, example01$m, example01$x, "logit", 1, 2, 0.01 )
}

\keyword{nonparametric}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}