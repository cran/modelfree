\name{binom_gl}
\alias{binom_gl}
\title{Psychometric function with guessing and lapsing rates}
\description{THIS IS AN INTERNAL FUNCTION: USE \code{BINOM_LIMS} FOR BEST RESULTS. Maximum likelihood estimates of the parameters of psychometric function with guessing and lapsing rates (GLM) or only guessing rate. The estimated parameters for the linear part are in vector b and the estimated limits are in lims.}

\usage{
binom_gl( r, m, x, link, p, K, initval )
}

\arguments{
\item{r}{number of successes in points \code{x}}
\item{m}{number of trials in points \code{x}}
\item{x}{stimulus levels}
\item{link}{link function}
\item{p}{degree of the polynomial to be fitted on the linear scale}
\item{K}{Power parameter in Weibull and reverse Weibull models}
\item{initval}{initial value for guessing and lapsing rates}
}

\value{
\item{value}{Object with 4 components: 
guessing: estimated guessing rate
lapsing: estimated lapsing rate
b: vector of estimated coefficients for the linear part
fit: glm object to be used in evaluation of fitted values}
}

\examples{
data( "01_Miranda" );
value <-binom_gl( example01$r, example01$m, example01$x, "logit", 1, 2, c( 0.01, 0.01 ) );
}

\keyword{nonparametric}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}