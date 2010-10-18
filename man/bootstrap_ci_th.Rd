\name{bootstrap_ci_th}
\alias{bootstrap_ci_th}

\title{Bootstrap estimate of confidence interval for threshold estimation}

\description{Finds bootstrap estimate of a confidence interval at a significant level alpha for the estimated threshold for the local polynomial estimation of psychometric function (PF) with guessing and lapsing rates specified in lims. Confidence interval is based on bootstrap percentiles}

\usage{
bootstrap_ci_th( TH, r, m, x, N, h0, alpha = 0.05, X = (max(x)-min(x))*(0:999)/999+min(x), link = c( "logit" ), guessing = 0, lapsing = 0, K = 2, p = 1, ker = c( "dnorm" ), maxiter = 50, tol = 1e-6 )
}

\arguments{
\item{TH}{required threshold level}
\item{r}{number of successes in points x}
\item{m}{number of trials in points x }
\item{x}{stimulus levels}
\item{N}{number of bootstrap replications}
\item{h0}{pilot bandwidth; if not specified, then the scaled plug-in bandwidth is used}
\item{alpha}{significance level of the confidence interval}
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
ci: confidence interval based on bootstrap percentiles
th0: threshold estimate}
}

\examples{
data( "01_Miranda" )
bwd <- 0.2959;
value <- bootstrap_ci_th( 0.5, example01$r, example01$m, example01$x, 10, bwd );
}

\keyword{nonparametric}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}