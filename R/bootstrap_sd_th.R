#' Bootstrap standard deviation for the slope of a psychometric function
#'
#' The function finds a bootstrap estimate of the standard deviation of the estimated
#' threshold for the local polynomial estimate of the psychometric function with
#' guessing and lapsing rates.
#'
#' @usage bootstrap_sd_th( TH, r, m, x, N, h0, X = (max(x)-min(x))*(0:999)/999+min(x),
#'                  link = "logit", guessing = 0, lapsing = 0, K = 2, p = 1,
#'                  ker = "dnorm", maxiter = 50, tol = 1e-6 )
#
# INPUT
#
#' @param  TH required threshold level
#' @param r number of successes at points x
#' @param m number of trials at points x
#' @param x stimulus levels
#' @param N number of bootstrap replications; N should be at least 200 for reliable results
#' @param h0 bandwidth
#
# OPTIONAL INPUT
#
#' @param X       (optional) set of values at which estimates of the psychometric function for the threshold estimation are to be obtained; if not given, 1000 equally spaced points from minimum to maximum of 'x' are used
#' @param link     (optional) name of the link function; default is "logit"
#' @param guessing (optional) guessing rate; default is 0
#' @param lapsing  (optional) lapsing rate; default is 0
#' @param K    (optional) power parameter for Weibull and reverse Weibull link; default is 2
#' @param p        (optional) degree of the polynomial; default is 1
#' @param ker      (optional) kernel function for weights; default is "dnorm"
#' @param maxiter  (optional) maximum number of iterations in Fisher scoring; default is 50
#' @param tol      (optional) tolerance level at which to stop Fisher scoring; default is 1e-6
#
# OUTPUT
#
#' @returns  \verb{sd   } bootstrap estimate of the standard deviation of the threshold estimator
#' @returns  \verb{th0  } threshold estimate
#'
#' @examples
#' \donttest{
#' data("Miranda_Henson")
#' x = Miranda_Henson$x
#' r = Miranda_Henson$r
#' m = Miranda_Henson$m
#' bwd_min <- min( diff( x ) )
#' bwd_max <- max( x ) - min( x )
#' bwd <- bandwidth_cross_validation( r, m, x, c( bwd_min, bwd_max ), method = "deviance" )
#' prob <- 0.5 # Required threshold level
##' # This might take a few minutes
#' niter <- 200 # Note number of bootstrap iterations should be at least 200
#' th_sl <- bootstrap_sd_sl( prob, r, m, x, niter, bwd ) # Be patient, slow process
#'}
#' @importFrom stats var
#' @export
bootstrap_sd_th <- function( TH, r, m, x, N, h0,
                             X = (max(x)-min(x))*(0:999)/999+min(x),
                             link = "logit", guessing = 0, lapsing = 0,
                             K = 2, p = 1, ker = "dnorm", maxiter = 50,
                             tol = 1e-6 ) {
#
# The function finds a bootstrap estimate of the standard deviation of the estimated
# threshold for the local polynomial estimate of the psychometric function with
# guessing and lapsing rates.
#
# INPUT
#
# TH   - required threshold level
# r    - number of successes at points x
# m    - number of trials at points x
# x    - stimulus levels
# N    - number of bootstrap replications; N should be at least 200 for reliable results
# h0   - bandwidth
#
# OPTIONAL INPUT
#
# X        - set of values at which estimates of the psychometric function for the threshold
#		estimation are to be obtained; if not given, 1000 equally spaced points from
#		minimum to maximum of 'x' are used
# link     - name of the link function; default is "logit"
# guessing - guessing rate; default is 0
# lapsing  - lapsing rate; default is 0
# K    - power parameter for Weibull and reverse Weibull link; default is 2
# p        - degree of the polynomial; default is 1
# ker      - kernel function for weights; default is "dnorm"
# maxiter  - maximum number of iterations in Fisher scoring; default is 50
# tol      - tolerance level at which to stop Fisher scoring; default is 1e-6
#
# OUTPUT
#
# Object with 2 components:
# sd  - bootstrap estimate of the standard deviation of the threshold
# estimator
# th0 - threshold estimate

# MAIN PROGRAM
# First 6 arguments are mandatory
    if( missing("TH") || missing("r") || missing("m") || missing("x") ||
        missing("N") || missing("h0") ) {
        stop("Check input. First 6 arguments are mandatory");
    }

# CHECK ROBUSTNESS OF INPUT PARAMETERS
    if( !is.double(TH) || length( TH ) > 1 ) {
        stop( "Threshold level must be scalar" );
    }

    checkdata<-list();
    checkdata[[1]] <- x;
    checkdata[[2]] <- r;
    checkdata[[3]] <- m;
    checkinput( "psychometricdata", checkdata );
    rm( checkdata )

    checkinput( "bootstrapreplications", N );

    if( !is.vector( X ) ) {
        stop("X (values where to estimate the PF) has to be a vector");
    }
    if( length( X ) < 2 ) {
        stop("At least 2 values needed for vector X");
    }
    if( min(x) > min(X) || max(x) < max(X) ) {
        stop("Supplied values of X are outside the range of stimulus levels x");
    }

    checkinput( "bandwidth", h0 );
    checkinput( "linkfunction", link );
    if( length( guessing ) > 1 ) {
        stop( "Guessing rate must be scalar" );
    }
    if( length( lapsing ) > 1 ) {
        stop( "Lapsing rate must be scalar" );
    }
    checkinput( "guessingandlapsing", c( guessing, lapsing ) );

    if( any( TH <= guessing )) {
        stop( "Threshold level should be greater than guessing rate" );
    }
	if( any( TH >= 1-lapsing ) ) {
        stop( "Threshold level should be smaller than 1 - lapsing rate" );
    }

    if (link == "weibull" || link == "revweibull"){
	    checkinput( "exponentk", K );
	    }
	pn <- list()
	pn[[1]] <- p
	pn[[2]] <- x
    checkinput( "degreepolynomial", pn );
    checkinput( "kernel", ker );
    checkinput( "maxiter", maxiter );
    checkinput( "tolerance", tol );

    if( N < 200 ) {
        warning( "number of bootstrap should be larger than 200\n otherwise results might be unreliable" );
    }

    n <- length( x );

# INITIAL ESTIMATE
# initial estimates with bandwidth h0
    f <- locglmfit( x, r, m, x, h0, FALSE, link, guessing, lapsing, K,
                    p, ker, maxiter, tol )$pfit;

# dense version for estimation of the threshold
    F <- locglmfit( X, r, m, x, h0, FALSE, link, guessing, lapsing, K,
                    p, ker, maxiter, tol )$pfit;

# THRESHOLD ESTIMATE
    sd <- NULL;
    th0 <- threshold_slope( F, X, TH )$x_th;

# BOOTSTRAP SAMPLING
# re-sampling
    samp <- matrix( 0, n, N );
	samp <- matrix(rbinom(N*n,m,f),n,N)

# exclude "degenerate samples" if min(M)>1
    if( min( m ) > 1 ) {
        ok <- NULL;
        for( i in 1:N ) {
            ok[i] = length( unique( samp[,i] ) ) > 3;
        }
        while( any( ok == FALSE ) ) {
            lis <- which( ok == 0 )
 			samp[,lis] <- matrix(rbinom(length(lis)*n,m,f),n,length(lis))
            for( i in 1:length( lis ) ) {
                ok[lis[i]] <- length( unique( samp[,lis[i]] ) ) > 3;
            }
        }
    }

# INITIATE VARIABLE IN WHICH DATA ARE STORED
    th_boot <- NULL;

# BOOTSTRAP ESTIMATES OF THE THRESHOLD

    for( i in 1:N ) {
        ftmp <- locglmfit_inter( X, samp[,i], m, x, h0, FALSE, link, guessing,
                           lapsing, K, p, ker, maxiter, tol )$pfit;
        th_boot[i] <- threshold_slope( ftmp, X, TH )$x_th;
     }


    sd <- sqrt( var( th_boot ) );

	value <- NULL
    value$sd <- sd
    value$th0 <- th0

    return( value );
}
