\name{locglmfit_sparse_private}
\alias{locglmfit_sparse_private}
\title{Local generalized linear fitting with sparse matrices}

\description{THIS IS AN INTERNAL FUNCTION: USE \code{LOCGLMFIT} FOR BEST RESULTS. Fisher scoring method for local polynomial estimator of a psychometric function (PF).}

\usage{
locglmfit_sparse_private( xfit, r, m, x, h, returnH, link, guessing, lapsing, K, p, ker, maxiter, tol )
}

\arguments{
\item{xfit}{points in which to calculate the estimate}
\item{r}{number of successes in points \code{x}}
\item{m}{number of trials in points \code{x}}
\item{x}{stimulus values}
\item{h}{bandwidths}
\item{returnH}{Boolean; Return or not the hat matrix H? default is \code{TRUE}}
\item{link}{name of the link function to be used; default is \code{"logit"}}
\item{guessing}{guessing rate; default is 0}
\item{lapsing}{lapsing rate; default is 0}
\item{K}{power parameter for Weibull and reverse Weibull link; default is 2}
\item{p}{degree of the polynomial; default \code{p = 1}}
\item{ker}{kernel function for weights; default \code{"dnorm"}}
\item{maxiter}{maximum number of iterations in Fisher scoring; default is 50}
\item{tol}{tolerance level at which to stop Fisher scoring; default is 1e-6}
}

\value{
\item{value}{Object with 2 or 3 components: 
pfit: value of the local polynomial estimate at points xfit
etafit: estimate of eta (link of pfit)
H: hat matrix (OPTIONAL)}
}

\examples{
data( "01_Miranda" )
xnew = 1.2 * (0:99)/99+0.1
h <- 0.2959
fit <- locglmfit_sparse_private( xnew,  example01$r,  example01$m, example01$x, h, FALSE, "logit_link", 0, 0, 2, 1, "dnorm", 50, 1e-6)
}

\keyword{nonparametric}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}