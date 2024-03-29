#' Log-log link with guessing and lapsing rates
#'
#' Creates a log-log link function; the guessing rate and lapsing rate are fixed, hence link is
#' a function of only one variable.
#
# OPTIONAL INPUT
#
#' @param guessing (optional) guessing rate; default is 0
#' @param lapsing (optional) lapsing rate; default is 0
#
# OUTPUT
#
#' @returns link log-log link for use in all GLM functions
loglog_link <-function( guessing = 0, lapsing = 0 ) {
#
# Log-log link for use with GLM functions
#
# The guessing rate and lapsing rate are fixed, hence link is
# a function of only one variable.
#
# OPTIONAL INPUT
#
# guessing - guessing rate; default is 0
# lapsing - lapsing rate; default is 0
#
# OUTPUT
#
# link - log-log link for use in all GLM functions

# LINK FUNCTIONS
# LOGLOG WITH LIMITS
# link
    loglogFL<-function( mu, g, l ) {
        mu <- pmax( pmin( l - .Machine$double.eps, mu) , g
              + .Machine$double.eps );
        return( -log( -log( ( mu - g ) / ( l - g ) ) ) );
    }

# derivative
    loglogFD<-function( eta, g, l ) {
        eta <- pmax( -log( -log( ( .Machine$double.eps ) / ( l - g ) ) ), eta );
        eta <- pmin( -log( -log( 1 - .Machine$double.eps / ( l - g ) ) ), eta );
        return( ( l - g ) * exp( -exp( -eta ) - eta  ) );
    }

# inverse link
    loglogFI<-function( eta, g, l ) {
        eta <- pmax( -log( -log( ( .Machine$double.eps ) / ( l - g ) ) ), eta );
        eta <- pmin( -log( -log( 1 - .Machine$double.eps / ( l - g ) ) ), eta );
        return( g + ( l - g ) * exp( -exp( -eta ) ) );
    }

# User-defined link
    linkuser <- function( guessing, lapsing ) {
        linkl <- "loglogFL";
        linkd <- "loglogFD";
        linki <- "loglogFI";
        linkfun <- function(mu)  eval( call( linkl, mu,  guessing, 1-lapsing ) );
        linkinv <- function(eta) eval( call( linki, eta, guessing, 1-lapsing ) );
        mu.eta  <- function(eta) eval( call( linkd, eta, guessing, 1-lapsing ) );
        link <- paste("loglog( ", "c(", guessing, ", ", 1-lapsing, ")", " )",
                sep = "" );
        structure(list(linkfun = linkfun, linkinv = linkinv,
                        mu.eta = mu.eta, name = link),
                  class = "link-glm" );
    }

# MAIN PROGRAM
# CHECK ROBUSTNESS OF INPUT PARAMETERS
    checkinput( 'guessingandlapsing', c( guessing, lapsing ) );

    return( linkuser( guessing, lapsing ) );
}
