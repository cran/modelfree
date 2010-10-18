\name{probit_link}
\alias{probit_link}

\title{Probit link function with guessing and lapsing rates}

\description{Probit link for use with GLM functions. The guessing and lapsing rate are fixed to values given in \code{lims}, hence \code{link} is a function of only one variable}

\usage{
probit_link( guessing = 0, lapsing = 0 )
}

\arguments{
\item{guessing}{guessing rate; default is 0}
\item{lapsing}{lapsing rate; default is 0}
}

\value{
\item{link}{Probit link for use in all GLM functions}
}

\examples{
data( "01_Miranda" )
x <- example01$x
r <- example01$r
m <- example01$m
glmdata <- data.frame( cbind( r/m ,m , x ) )
names( glmdata ) <- c( "resp", "m", "x" )
glmformula <- c( "resp ~ x" )
userlink<-probit_link( 0.1, 0.1 )
fit <- glm( glmformula, data = glmdata, weights = m, family = binomial( userlink ) )
}

\keyword{nonparametric}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}
