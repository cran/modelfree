#'Weibull model for the psychometric function
#'
#' This function finds the maximum likelihood estimates of the parameters
#' of the Weibull model for the psychometric function.
#'
#' @usage binom_weib( r, m, x, p = 1, initK = 2, guessing = 0, lapsing = 0 )
#
# INPUT
#
#' @param r number of successes at points x
#' @param m number of trials at points x
#' @param x stimulus levels
#
# OPTIONAL INPUT
#
#' @param p    (optional) degree of the polynomial; default is 1
#' @param initK (optional) initial value for K (power parameter in Weibull model); default is 2
#' @param guessing (optional) guessing rate; default is 0
#' @param lapsing  (optional) lapsing rate; default is 0
#
# OUTPUT
#
#' @returns \verb{b    } vector of estimated coefficients for the linear part
#' @returns \verb{K    } estiamte of the power parameter in the Weibull model
#' @returns \verb{fit  } glm object to be used in evaluation of fitted values
#'
#' @examples
#' data("Miranda_Henson")
#' x = Miranda_Henson$x
#' r = Miranda_Henson$r
#' m = Miranda_Henson$m
#' numxfit <- 199; # Number of new points to be generated minus 1
#' xfit <- (max(x)-min(x)) * (0:numxfit) / numxfit + min(x)
#' val <- binom_weib( r, m, x )
#' # Plot the fitted curve
#' plot( x, r / m, xlim = c( 0.1, 1.302 ), ylim = c( 0.0165, 0.965 ), type = "p", pch="*" )
#' pfit <- predict( val$fit, data.frame( x = xfit ), type = "response" )
#' lines(xfit, pfit, col = "red" )
#'
#' @export
binom_weib <- function( r, m, x, p = 1, initK = 2, guessing = 0, lapsing = 0 ) {
#
# This function finds the maximum likelihood estimates of the parameters
# of the Weibull model for the psychometric function.
#
# INPUT
#
# r    - number of successes at points x
# m    - number of trials at points x
# x    - stimulus levels
#
# OPTIONAL INPUT
#
# p     - degree of the polynomial; default is 1
# initK - initial value for K (power parameter in Weibull model); default is 2
# guessing - guessing rate; default is 0
# lapsing  - lapsing rate; default is 0
#
# OUTPUT
#
# Object with 3 components:
# b - vector of estimated coefficients for the linear part
# K - estiamte of the power parameter in the Weibull model
# fit - glm object to be used in evaluation of fitted values

# LIKELIHOOD
    likfun <- function( K ) {
        K = 0.05 + exp( K );

# fit
        fit <- glm( glmformula, data = glmdata, weights = m,
               family = binomial( weibull_link( K, guessing, lapsing ) ) );

# FITTED PROBABILITIES
        fitted <- as.numeric( predict( fit, type = "response" ) );
        fitted[which( fitted <= guessing )] <- guessing + .Machine$double.eps;
        fitted[which( fitted >= 1 - lapsing )] <- 1 - lapsing - .Machine$double.eps;

        return( c( -( t( r ) %*% log( fitted ) + t( m - r ) %*%
                log( 1 - fitted ) ) ) );
    }

# MAIN PROGRAM
# First 3 arguments are mandatory
    if( missing("x") || missing("r") || missing("m") ) {
        stop("Check input. First 3 arguments are mandatory");
    }

# CHECK ROBUSTNESS OF INPUT PARAMETERS
    checkdata<-list();
    checkdata[[1]] <- x;
    checkdata[[2]] <- r;
    checkdata[[3]] <- m;
    checkinput( "psychometricdata", checkdata );
    rm( checkdata )
    pn <- list()
	pn[[1]] <- p
	pn[[2]] <- x
    checkinput( "degreepolynomial", pn );
    checkinput( "guessingandlapsing", c( guessing, lapsing ) );
    checkinput( 'exponentk', initK );

# GLM settings
    glmdata <- data.frame( cbind( r/m ,m , x ) );
    names( glmdata ) <- c( "resp", "m", "x" );

# formula
    glmformula <- c( "resp ~ x" );
    if( p > 1 ) {
        for( pp in 2:p ) {
            glmformula <- paste( glmformula, " + I(x^", pp,")", sep = "");
        }
    }
    fit <- NULL;

	initK <- log( initK )
    suppressWarnings( K <- optim( initK, likfun )$par );
    K <- 0.05 + exp( K );

    fit <- glm( glmformula, data = glmdata, weights = m,
           family = binomial( weibull_link( K, guessing, lapsing ) ) );

    fit$df.residual <- length(x) - (p + 1) - 1

    b <- fit$coeff

        value <- NULL
    value$b <- b
    value$K <- K
    value$fit <- fit

    return( value );
}
