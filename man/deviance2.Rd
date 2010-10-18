\name{deviance2}
\alias{deviance2}

\title{Deviance between data and fitted function}

\description{Calculates deviance for data Y and fitted values of the psychometric function f}

\usage{
D<-deviance2(r, m, pfit)
}

\arguments{
\item{r}{number of successes in points x}
\item{m}{number of trials in points x}
\item{pfit}{fitted values}
}

\value{
\item{D}{Deviance}
}

\examples{
data( "01_Miranda" );
h = 0.2959;
fit <- locglmfit( example01$x, example01$r, example01$m, example01$x, h );
Dev <- deviance2( example01$r, example01$m, fit$fitval );
}

\keyword{nonparametric}