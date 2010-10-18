\name{revweibull_link_private}
\alias{revweibull_link_private}

\title{Reverse Weibull link function with guessing and lapsing rates}

\description{THIS IS AN INTERNAL FUNCTION: USE \code{REVWEIBULL_LINK} FOR BEST RESULTS. Reverse Weibull link for use with GLM functions. The guessing rate and lapsing rate are fixed, and power parameter is set to be equal K, hence link is a function of only one variable}

\usage{
link<-revweibull_link_private( K, guessing = 0, lapsing = 0 );
}

\arguments{
\item{K}{power parameter for reverse Weibull link function}
\item{guessing}{guessing rate; default is 0}
\item{lapsing}{lapsing rate; default is 0}
}

\value{
\item{link}{Reverse Weibull link for use in all GLM functions}
}

\examples{
data( "01_Miranda" );
x <- example01$x;
r <- example01$r;
m <- example01$m;
glmdata <- data.frame( cbind( r/m ,m , x ) );
names( glmdata ) <- c( "resp", "m", "x" );
glmformula <- c( "resp ~ x" );
userlink<-revweibull_link_private( 20, 0, 0 );
fit <- glm( glmformula, data = glmdata, weights = m, family = binomial( userlink ) );
}

\keyword{nonparametric}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}