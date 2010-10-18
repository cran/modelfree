\name{threshold_slope}
\alias{threshold_slope}

\title{Threshold and slope of estimated psychometric function}

\description{Finds the approximate value of \code{x (=x_th)} for which the value of the psychometric function is equal thresh and the approximate value of slope in \code{x_th}}

\usage{
value<-threshold_slope( pfit, xfit, thresh = 0.5 );
}

\arguments{
\item{pfit}{estimated values of the psychometric function}
\item{xfit}{stimulus levels in which the function was estimated}
\item{thresh}{value for which to estimate threshold; default is 0.5}
}

\value{
\item{value}{Object with 2 elements:
x_th: estimated threshold
slope: estimated value of slope, i.e. derivative of pfit at x_th}
}

\examples{
data( "01_Miranda" );
xnew = 1.2 * (0:999)/999+0.1;
h = 0.2959;
fit <- locglmfit( xnew, example01$r, example01$m, example01$x, h );
value <- threshold_slope( fit$pfit , xnew );
}

\keyword{nonparametric}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}