\name{bandwidth_cross_validation}
\alias{bandwidth_cross_validation}
\title{Cross-validation estimate of bandwidth}
\description{Finds the cross-validation bandwidth for the local polynomial estimator of the psychometric function (PF) with guessing and lapsing rates specified in lims.}

\usage{
h<-bandwidth_cross_validation( r, m, x, H, link = "logit", guessing = 0, lapsing = 0, K = 2, p = 1, ker = "dnorm", maxiter = 50, tol = 1e-6, method = "all");
}
\arguments{
\item{r}{number of successes in points \code{x}}
\item{m}{number of trials in points \code{x}}
\item{x}{stimulus levels}
\item{H}{minimum and maximum values of bandwidth to be considered}
\item{link}{name of the link function to be used; default is "logit"}
\item{guessing}{guessing rate; default is 0}
\item{lapsing}{lapsing rate; default is 0}
\item{K}{power parameter for Weibull and reverse Weibull link; default is 2}
\item{p}{order of the polynomial; default is 1}
\item{ker}{kernel function for weights; default "dnorm"}
\item{maxiter}{maximum number of iterations in Fisher scoring; default is 50}
\item{tol}{tolerance level at which to stop Fisher scoring; default is 1e-6}
\item{method}{loss function to be used in cross-validation: choose from: 'ISEeta', 'ISE', 'deviance'; default is "all"}}
\value{
\item{h}{cross-validation bandwidth for the chosen method; if no method was specified, then it a list of three elements with entries corresponding to the estimated bandwidths on p-scale (\code{h$pscale}), on eta-scale (\code{h$etascale}) and for mean likelihood (\code{h$deviance})}
}

\examples{
data("01_Miranda");
h<- bandwidth_cross_validation( example01$r, example01$m, example01$x, c( 0.1, 10 ) );
}

\keyword{nonparametric}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}