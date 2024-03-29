% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binomfit_lims.R
\name{binomfit_lims}
\alias{binomfit_lims}
\title{Parameter estimates for a psychometric function with FIXED guessing and lapsing rates}
\usage{
binomfit_lims( r, m, x, p = 1, link = "logit", guessing = 0, lapsing = 0, K = 2 )
}
\arguments{
\item{r}{number of successes at points x}

\item{m}{number of trials at points x}

\item{x}{stimulus levels}

\item{p}{(optional) degree of the polynomial; default is p = 1}

\item{link}{(optional) name of the link function; default is "logit"}

\item{guessing}{(optional) guessing rate; default is 0}

\item{lapsing}{(optional) lapsing rate; default is 0}

\item{K}{(optional) power parameter for Weibull and reverse Weibull link; default is 2}
}
\value{
\verb{b    } vector of estiamted coefficients for the linear part

\verb{fit  } glm object to be used in evaluation of fitted values
}
\description{
This function fits a binomial generalised linear model with fixed guessing and lapsing rates.
}
\examples{
data("Carcagno")
x = Carcagno$x
r = Carcagno$r
m = Carcagno$m
plot( x, r / m, xlim = c( 1.95, 4.35 ), ylim = c( 0.24, 0.99 ), type = "p", pch="*" )
guess = 1/3; # guessing rate
laps = 0; # lapsing rate
val <- binomfit_lims( r, m, x, link = "probit", guessing = guess, lapsing = laps )
numxfit <- 199 # Number of new points to be generated minus 1
xfit <- (max(x)-min(x)) * (0:numxfit) / numxfit + min(x)
# Plot the fitted curve
pfit<-predict( val$fit, data.frame( x = xfit ), type = "response" )
lines(xfit, pfit )

}
