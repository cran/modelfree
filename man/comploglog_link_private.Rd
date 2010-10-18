\name{comploglog_link_private}
\alias{comploglog_link_private}

\title{Complementary log-log link function with guessing and lapsing rates}

\description{THIS IS AN INTERNAL FUNCTION: USE \code{COMPLOGLOG_LINK} FOR BEST RESULTS. Complementaty loglog link for use with GLM functions. The guessing and lapsing rate are fixed to values given in \code{lims}, hence \code{link} is a function of only one variable}

\usage{
comploglog_link_private( guessing, lapsing )
}

\arguments{
\item{guessing}{guessing rate}
\item{lapsing}{lapsing rate}
}

\value{
\item{link}{Complementary log-log link for use in all GLM functions}
}

\examples{
data( "01_Miranda" )
x <- example01$x
r <- example01$r
m <- example01$m
glmdata <- data.frame( cbind( r/m ,m , x ) )
names( glmdata ) <- c( "resp", "m", "x" )
glmformula <- c( "resp ~ x" )
userlink<-comploglog_link_private( 0.1, 0.1 )
fit <- glm( glmformula, data = glmdata, weights = m, family = binomial( userlink ) )
}

\keyword{nonparametric}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}