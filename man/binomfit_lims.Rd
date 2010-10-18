\name{binomfit_lims}
\alias{binomfit_lims}

\title{Generalized linear model fit with guessing and lapsing rates}

\description{The function fits a binomial generalised liner model with fixed guessing and lapsing rates}

\usage{
value <- binomfit_lims( r, m, x, p = 1, link = "logit", guessing = 0, lapsing = 0, K = 2 )
}

\arguments{
\item{r}{number of successes in points x}
\item{m}{number of trials in points x }
\item{x}{stimulus levels}
\item{p}{degree of the polynomial to be fitted on the linear scale; default is 1}
\item{link}{link function; default is "logit"}
\item{guessing}{guessing rate; default is 0}
\item{lapsing}{lapsing rate; default is 0}
}

\value{
\item{value}{Object with 2 components: 
b: vector of estimated coefficients for the linear part
fit: glm object to be used in evaluation of fitted values}
}

\examples{
data( "01_Miranda" );
value <- binomfit_lims( example01$r, example01$m, example01$x )
}

\keyword{nonparametric}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}