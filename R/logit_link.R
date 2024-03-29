#' Logit link with guessing and lapsing rates
#
#' Creates a complementary logit link function; the guessing rate and lapsing rate are fixed, hence link is
#' a function of only one variable.
#
# OPTIONAL INPUT
#
#' @param guessing (optional) guessing rate; default is 0
#' @param lapsing (optional) lapsing rate; default is 0
#
# OUTPUT
#
#' @returns link logit link for use in all GLM functions
logit_link <- function( guessing = 0, lapsing = 0 ) {
#
# Logit link for use with GLM functions
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
# link - logit link for use in all GLM functions

# LINK FUNCTIONS
# LOGIT WITH LIMITS
# link
    logitFL<-function( mu, g, l ) {
        mu<-pmax( pmin( l - .Machine$double.eps, mu ), g
            + .Machine$double.eps );
        return( log( ( mu - g ) / ( l - mu ) ) );
    }

# derivative
    logitFD<-function( eta, g, l ) {
        eta <- pmax( log( .Machine$double.eps / ( l - g ) ), eta );
        eta <- pmin( log( ( l - g ) / .Machine$double.eps ), eta );
        return( (l - g) * exp( -eta ) / ( 1 + exp( -eta ) )^2 );
    }

# inverse link
    logitFI<-function( eta, g, l ) {
        eta <- pmax( log( .Machine$double.eps / ( l - g ) ), eta );
        eta <- pmin( log( ( l - g ) / .Machine$double.eps ), eta );
        return( g + ( l - g ) / ( 1 + exp( -eta ) ) );
    }

# User-defined link
    linkuser <- function( guessing, lapsing ) {
        linkl <- "logitFL";
        linkd <- "logitFD";
        linki <- "logitFI";
        linkfun <- function(mu)  eval( call( linkl, mu,  guessing, 1-lapsing ) );
        linkinv <- function(eta) eval( call( linki, eta, guessing, 1-lapsing ) );
        mu.eta  <- function(eta) eval( call( linkd, eta, guessing, 1-lapsing ) );
        link <- paste("logit( ", "c(", guessing, ", ", 1-lapsing, ")", " )",
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
