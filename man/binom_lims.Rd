\name{binom_lims}
\alias{binom_lims}
\title{Psychometric function with guessing and lapsing rates}
\description{Maximum likelihood estimates of the parameters of psychometric function with guessing and lapsing rates (GLM) or only guessing rate. The estimated parameters for the linear part are in vector b and the estimated limits are in lims.}

\usage{
fit<-binom_lims( r, m, x, gl = "both", link = "logit", p = 1, K = 2, initval = NULL );
}
\arguments{
\item{r}{number of successes in points x}
\item{m}{number of trials in points x }
\item{x}{stimulus levels}
\item{gl}{indicator, calulate only guessing if "guessing", only lapsing if "lapsing" and both guessing and lapsing if "both"; default is "both"}
\item{link}{link function; default is "logit"}
\item{p}{degree of the polynomial to be fitted on the linear scale; default is 1}
\item{K}{Power parameter in Weibull and reverse Weibull models; default is 2}
\item{initval}{initial value for guessing and lapsing; default is c( 0.01, 0.01) if guessing and lapsing rates are estimated, and 0.01 if only guessing or only lapsing rate is estimated}
}
\value{
\item{value}{The object returned by internal functions binom_g, binom_l, or binom_gl, depending on only guessing, only lapsing, or guessing and lapsing are estimated}
}

\examples{
data( "01_Miranda" );
value <-binom_lims( example01$r, example01$m, example01$x );
}

\keyword{nonparametric}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}
