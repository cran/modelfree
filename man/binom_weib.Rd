\name{binom_weib}
\alias{binom_weib}
\title{Psychometric function fitting for Weibull link function}

\description{Maximum likelihood estimates of the parameters of the Weibull model (GLM). The estimated parameters for the linear part are in vector \code{b} and the estimated exponent is \code{K}.
}

\usage{
value<-binom_weib( r, m, x, p = 1, initK = 2, guessing = 0, lapsing = 0 )
}

\arguments{
\item{r}{number of successes in points x}
\item{m}{number of trials in points x }
\item{x}{stimulus levels}
\item{p}{degree of the polynomial to be fitted on the linear scale; default is 1}
\item{initK}{Power parameter in reverse Weibull model; default is 2}
\item{guessing}{guessing rate; default is 0}
\item{lapsing}{lapsing rate; default is 0}
}

\value{
\item{value}{Object with 3 components: 
b: vector of estimated coefficients for the linear part
K: estiamte of the power parameter in the reverse Weibull model
fit: glm object to be used in evaluation of fitted values}
}

\examples{
data( "01_Miranda" );
value <- binom_revweib( example01$r, example01$m, example01$x )
}

\keyword{nonparametric}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}
