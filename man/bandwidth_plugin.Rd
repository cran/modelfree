\name{bandwidth_plugin}
\alias{bandwidth_plugin}
\title{Plug in estimation of Bandwidth}
\description{Calculates an estimate of the AMISE optimal bandwidth for the local polynomial estimator of a psychometric function.}

\usage{
bandwidth_plugin( r, m, x, link = c( "logit" ), guessing = 0, lapsing = 0, K = 2, p = 1, ker = c( "dnorm" ) )
}

\arguments{
\item{r}{number of successes in points \code{x}}
\item{m}{number of trials in points \code{x}}
\item{x}{stimulus levels}
\item{link}{link function; default is "logit"}
\item{guessing}{guessing rate; default is 0}
\item{lapsing}{lapsing rate; default is 0}
\item{K}{power parameter for Weibull and reverse Weibull link; default is 2}
\item{p}{order of the polynomial; default is 1}
\item{ker}{kernel function for weights; default "dnorm"}
}

\value{
\item{h}{plug-in bandwidth (on eta-scale)}
}

\examples{
data("01_Miranda")
h<-bandwidth_plugin( example01$r, example01$m, example01$x )
}

\keyword{nonparametric}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}