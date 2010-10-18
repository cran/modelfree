\name{locglmfit}
\alias{locglmfit}

\title{Local generalized linear fitting}

\description{Local polynomial estimator for the psychometric function (PF) and eta function (PF transformed by link) for binomial data; also returns the Hat matrix. Actual calculations are done in \code{LOCGLMFIT_PRIVATE} or \code{LOCGLMFIT_SPARSE_PRIVATE} depending on the size of the data set. Here the data are split into several parts to speed up the calculations.}

\usage{
value<-locglmfit( xfit, r, m, x, h, returnH = FALSE, link = "logit", guessing = 0, lapsing = 0, K = 2, p = 1, ker = "dnorm", maxiter = 50, tol = 1e-6 )
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
data( "01_Miranda" );
xnew = 1.2 * (0:99)/99+0.1;
h <- 0.2959;
fit <- locglmfit( xnew, example01$r, example01$m, example01$x, h );
}

\keyword{nonparametric}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}