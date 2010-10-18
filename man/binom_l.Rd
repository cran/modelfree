\name{binom_l}
\alias{binom_l}

\title{Psychometric function with lapsing rate}

\description{
THIS IS AN INTERNAL FUNCTION: USE \code{BINOM_LIMS} FOR BEST RESULTS. Maximum likelihood estimates of the parameters of psychometric function with lapsing rate (GLM). The estimated parameters for the linear part are in vector b and the estimated guessing rate is guess (GLM).}

\usage{
value<-binom_l( r, m, x, link, p, K, initval );
}
\arguments{
\item{r}{number of successes in points \code{x}}
\item{m}{number of trials in points \code{x}}
\item{x}{stimulus levels}
\item{link}{link function}
\item{p}{degree of the polynomial to be fitted on the linear scale}
\item{K}{Power parameter in Weibull and reverse Weibull models}
\item{initval}{initial value for lapsing rate}
}
\value{
\item{value}{Object with 4 components: 
lapsing: estimated lapsing rate
b: vector of estimated coefficients for the linear part
fit: glm object to be used in evaluation of fitted values}
}

\examples{
data( "01_Miranda" );
value <-binom_l( example01$r, example01$m, example01$x, "logit", 1, 2, 0.01 );
}

\keyword{nonparametric}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}