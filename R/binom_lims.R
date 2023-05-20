#' Maximum likelihood parameter estimates for a psychometric function with guessing and lapsing rates
#'
#' This function finds the maximum likelihood estimates of the parameters
#' of the psychometric function with guessing and lapsing rates, only
#' guessing rate, or only lapsing rate.
#'
#' @usage binom_lims( r, m, x, gl = "both", link = "logit", p = 1, K = 2, initval = NULL )
#
# INPUT
#
#' @param  r number of successes at points x
#' @param m number of trials at points x
#' @param x stimulus levels
#
# OPTIONAL INPUT
#
#' @param gl     (optional) indicator, calulate only guessing if "guessing", only lapsing if "lapsing" and both guessing and lapsing if "both"; default is "both"
#' @param link    (optional) name of the link function; default is "logit"
#' @param p       (optional) degree of the polynomial; default is 1
#' @param K       (optional) power parameter for Weibull and reverse Weibull link; default is 2
#' @param initval (optional) initial value for guessing and lapsing; default is c(.01 .01) if guessing and rates are estimated, and .01 if only guessing or only lapsing rate is estimated
#
# OUTPUT
#
#' @returns \verb{b         }estimated coefficients for the linear part
#' @returns \verb{guessing  } estimated guessing rate (if estimated)
#' @returns \verb{lapsing   } estimated lapsing rate (if estimated)
#' @returns \verb{fit       } glm object to be used in evaluation of fitted values
#' @examples
#' data("Baker_etal")
#' x = Baker_etal$x
#' r = Baker_etal$r
#' m = Baker_etal$m
#' plot( x, r / m, xlim = c( 0.16, 7.83 ), ylim = c( -0.01, 1.01 ), type = "p", pch="*" )
#' val <- binomfit_lims( r, m, x, link = "probit" )
#' numxfit <- 199; # Number of new points to be generated minus 1
#' xfit <- (max(x)-min(x)) * (0:numxfit) / numxfit + min(x)
#' # Plot the fitted curve
#' pfit<-predict( val$fit, data.frame( x = xfit ), type = "response" )
#' lines(xfit, pfit )
#' @export
binom_lims<-function( r, m, x, gl = "both", link = "logit", p = 1, K = 2,
                      initval = NULL ) {
#
# This function finds the maximum likelihood estimates of the parameters
# of the psychometric function with guessing and lapsing rates, only
# guessing rate, or only lapsing rate.
#
# INPUT
#
# r    - number of successes at points x
# m    - number of trials at points x
# x    - stimulus levels
#
# OPTIONAL INPUT
#
# gl      - indicator, calulate only guessing if "guessing", only lapsing if "lapsing"
#  and both guessing and lapsing if "both"; default is "both"
# link    - name of the link function; default is "logit"
# p       - degree of the polynomial; default is 1
# K       - power parameter for Weibull and reverse Weibull link; default is 2
# initval - initial value for guessing and lapsing; default is c(.01 .01) if guessing and
# lapsing rates are estimated, and .01 if only guessing or only lapsing rate is estimated
#
# OUTPUT
#
# Object with 3 or 4 components:
# b - estiamted coefficients for the linear part
# guessing - estimated guessing rate (if estimated)
# lapsing - estimated lapsing rate (if estimated)
# fit - glm object to be used in evaluation of fitted values

# MAIN PROGRAM
# First 3 arguments are mandatory
    if( missing("x") || missing("r") || missing("m") ) {
        stop("Check input. First 3 arguments are mandatory");
    }

    if( gl != "both"  && gl != "guessing" && gl != "lapsing") {
        stop( "Wrong value for guessing/lapsing indicator gl" );
    }
    if(is.null(initval)){

    if( gl == "guessing" || gl == "lapsing" ){
    	initval <- .01;
    	}else{
    		initval <- c(.01,.01);
    		}
    }

# CHECK ROBUSTNESS OF INPUT PARAMETERS
    checkdata<-list();
    checkdata[[1]] <- x;
    checkdata[[2]] <- r;
    checkdata[[3]] <- m;
    checkinput( "psychometricdata", checkdata );
    rm( checkdata )
    checkinput( "linkfunction", link );
    pn <- list()
	pn[[1]] <- p
	pn[[2]] <- x
    checkinput( "degreepolynomial", pn );
    if( link == "weibull"  || link == "revweibull") {
    	checkinput( 'exponentk', K );
    	}

    initval[initval==0] <- .Machine$double.eps


# Check initial values for guessing or guessing/lapsing and call internal
# functions binom_gl or binom_g depending on indicator gl
    if( gl == "both" ) {
    	if( length(initval) != 2){
    		stop( "initval must have two elements if both guessing and lapsing rates are being estimated")
    		}
        checkinput( "guessingandlapsing", initval );
        return( binom_gl( r, m, x, link, p, K, initval ) );
    }
    else {
    	if( gl == "guessing" ) {
# Check that initval is a positive scalar
        if( length( initval ) > 1 || initval <= 0 || initval >= 1 ) {
            stop( "Guessing rate must be a scalar between 0 and 1" );
        }
        return( binom_g( r, m, x, link, p, K, initval ) );
        }# end guessing
        else{
        	if( length( initval ) > 1 || initval <= 0 || initval >= 1 ) {
            stop( "Lapsing rate must be a scalar between 0 and 1" );
	        }
    	    return( binom_l( r, m, x, link, p, K, initval ) );
        	}# end lapsing

        }

}
