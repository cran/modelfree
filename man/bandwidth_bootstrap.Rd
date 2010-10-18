\name{bandwidth_bootstrap}
\alias{bandwidth_bootstrap}
\title{Bootstrap estimate of bandwidth}
\description{Finds bootstrap estimate of the optimal bandwidth \code{h} for binomial data in local polynomial estimation of psychometric function (PF) with guessing and lapsing rates specified in lims.}

\usage{
bandwidth_bootstrap( r, m, x, H, N, h0 = NULL, link = c( "logit" ), guessing = 0, lapsing = 0, K = 2, p = 1, ker = c( "dnorm" ), maxiter = 50, tol = 1e-6, method = c( "all" ) )
}
\arguments{
\item{r}{number of successes in points \code{x}}
\item{m}{number of trials in points \code{x}}
\item{x}{stimulus levels}
\item{H}{minimum and maximum values of bandwidth to be considered}
\item{N}{number of bootstrap replications}
\item{h0}{pilot bandwidth; if not specified, then the scaled plug-in bandwidth is used}
\item{link}{name of the link function to be used; default is "logit"}
\item{guessing}{guessing rate; default is 0}
\item{lapsing}{lapsing rate; default is 0}
\item{K}{power parameter for Weibull and reverse Weibull link; default is 2}
\item{p}{order of the polynomial; default is 1}
\item{ker}{kernel function for weights; default "dnorm"}
\item{maxiter}{maximum number of iterations in Fisher scoring; default is 50}
\item{tol}{tolerance level at which to stop Fisher scoring; default is 1e-6}
\item{method}{loss function to be used in bandwidth: choose from: 'ISEeta', 'ISE', 'deviance'; default is "all"}
}
\value{
\item{h}{bootstrap bandwidth for the chosen method; if no method was specified, then it a list of three elements with entries corresponding to the estimated bandwidths on p-scale (\code{h$pscale}), on eta-scale (\code{h$etascale}) and for mean likelihood (\code{h$deviance})}
}

\examples{
data("01_Miranda")
h<- bandwidth_bootstrap( example01$r, example01$m, example01$x, c( 0.1, 10 ), 10 )
}

\keyword{nonparametric}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}