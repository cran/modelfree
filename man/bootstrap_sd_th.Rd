\name{bootstrap_sd_th}
\alias{bootstrap_sd_th}

\title{Bootstrap estimate the standard deviation of threshold estimation}

\description{Finds bootstrap estimate of the standard deviation of the estimated thresholdfor the local polynomial estimation of psychometric function (PF) with guessing and lapsing rates as specified}

\usage{
value<-bootstrap_sd_th( TH, r, m, x, N, h0, X = (max(x)-min(x))*(0:999)/999+min(x), link = "logit", guessing = 0, lapsing = 0, K = 2, p = 1, ker = "dnorm", maxiter = 50, tol = 1e-6 );
}

\arguments{
\item{TH}{required threshold level}
\item{r}{number of successes in points x}
\item{m}{number of trials in points x }
\item{x}{stimulus levels}
\item{N}{number of bootstrap replications}
\item{h0}{pilot bandwidth; if not specified, then the scaled plug-in bandwidth is used}
\item{X}{set of value for which to calculate the estimates of PF for the thresholdestimation; if not given 1000 equally spaced points from min to max of xdes are used}
\item{link}{name of the link function to be used; default is "logit"}
\item{guessing}{guessing rate; default is 0}
\item{lapsing}{lapsing rate; default is 0}
\item{K}{power parameter for Weibull and reverse Weibull link; default is 2}
\item{p}{order of the polynomial; default is 1}
\item{ker}{kernel function for weights; default "dnorm"}
\item{maxiter}{maximum number of iterations in Fisher scoring; default is 50}
\item{tol}{tolerance level at which to stop Fisher scoring; default is 1e-6}
}

\value{
\item{value}{Object with 2 components: 
sd: bootstrap estimate of the standard deviation of the threshold estimate
th0: threshold estimate}
}

\examples{
data( "01_Miranda" )
bwd <- 0.2959;
value <- bootstrap_sd_th( 0.5, example01$r, example01$m, example01$x, 10, bwd );
}

\keyword{nonparametric}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}