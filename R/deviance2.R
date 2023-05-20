#' Deviance of a psychometric function
#'
#' This function calculates the deviance for the fitted values of the psychometric function pfit.
#' @usage deviance2( r, m, pfit )
#
#
# INPUT
#' @param r     number of successes
#' @param m     number of trials
#' @param pfit  fittd values
#
# OUTPUT
#' @returns  \verb{D  } deviance
#'
#' @examples
#' data("Carcagno")
#' x = Carcagno$x
#' r = Carcagno$r
#' m = Carcagno$m
#' plot( x, r / m, xlim = c( 1.95, 4.35 ), ylim = c( 0.24, 0.99 ), type = "p", pch="*" )
#' guess = 1/3; # guessing rate
#' laps = 0; # lapsing rate
#' val <- binomfit_lims( r, m, x, link = "probit", guessing = guess, lapsing = laps )
#' pfit<-predict( val$fit, data.frame( x = x ), type = "response" )
#' d2 = deviance2( r, m, pfit )
#'
#' @export
deviance2 <- function( r, m, pfit ) {
#
# The function calculates the deviance for the fitted values of the psychometric function pfit.
#
#
# INPUT
# r    - number of successes
# m    - number of trials
# pfit - fittd values
#
# OUTPUT
# D - deviance

# Both arguments are mandatory
    if( missing("pfit") || missing("r") || missing("m") ) {
        stop("Check input. First 3 arguments are mandatory");
    }

# adjustment to avoid degenerate values
    r[which(r >= m)]<-r[which(r >= m)] - .001;
    r[which(r <= 0)]<-.001;

    pfit[which(pfit >= 1)]<-1 - .001;
    pfit[which(pfit <= 0)]<-.001;

# deviance
    return(2 * sum( ( r * log( r / ( m * pfit) ) + ( m - r ) * log( ( m - r ) / ( m - m * pfit ) ) ) ));
}
