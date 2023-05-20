#' Bootstrap confidence interval for the slope of a psychometric function
#'
#' Finds a bootstrap estimate of a confidence interval at a significance level alpha
#' for the estimated slope for the local polynomial estimate of the psychometric function
#' with guessing and lapsing rates. The confidence interval is based on bootstrap percentiles.
#
#' See Efron & Tibshirani's "An introduction to the bootstrap", 1993
#'
#' @usage bootstrap_ci_sl( TH, r, m, x, N, h0, alpha = 0.05,
#'                  X = (max(x)-min(x))*(0:999)/999+min(x), link = "logit", guessing = 0,
#'                  lapsing = 0, K = 2, p = 1, ker = "dnorm", maxiter = 50, tol = 1e-6 )
#
# INPUT
#
#' @param  TH   required threshold level
#' @param  r    number of successes at points x
#' @param  m    number of trials at points x
#' @param  x    stimulus levels
#' @param  N    number of bootstrap replications; N should be at least 1000 for reliable results
#' @param  h0   bandwidth
#
# OPTIONAL INPUT
#
#' @param  alpha    (optional) significance level of the confidence interval; default is 0.05
#' @param  X        (optional) set of values at which estimates of the psychometric function for the slope estimation are to be obtained; if not given, 1000 equally spaced points from minimum to maximum of 'x' are used
#' @param  link     (optional) name of the link function; default is "logit"
#' @param  guessing (optional) guessing rate; default is 0
#' @param  lapsing  (optional) lapsing rate; default is 0
#' @param  K    (optional) power parameter for Weibull and reverse Weibull link; default is 2
#' @param  p        (optional) degree of the polynomial; default is 1
#' @param  ker      (optional) kernel function for weights; default is "dnorm"
#' @param  maxiter  (optional) maximum number of iterations in Fisher scoring; default is 50
#' @param  tol      (optional) tolerance level at which to stop Fisher scoring; default is 1e-6
#
# OUTPUT
#
#' @returns  \verb{ci   } confidence interval based on bootstrap percentiles
#' @returns  \verb{sl0  } slope estimate
#'
#'
#' @examples
#' \dontrun{
#' data("Miranda_Henson")
#' x = Miranda_Henson$x
#' r = Miranda_Henson$r
#' m = Miranda_Henson$m
#' bwd_min <- min( diff( x ) )
#' bwd_max <- max( x ) - min( x )
#' bwd <- bandwidth_cross_validation( r, m, x, c( bwd_min, bwd_max ), method = "deviance" )
#' prob <- 0.5 # Required threshold level
#' alpha <- 0.05 # Significance level for the confidence intervals
#' # This might take a few minutes
#' niter <- 1000 # Note number of bootstrap iterations should be at least 1000
#' ci_sl <- bootstrap_ci_sl( prob, r, m, x, niter, bwd, alpha ) # Be patient, slow process
#'}
#'
#' @importFrom stats quantile
#' @export
bootstrap_ci_sl <- function( TH, r, m, x, N, h0, alpha = 0.05,
                             X = (max(x)-min(x))*(0:999)/999+min(x),
                             link = "logit", guessing = 0, lapsing = 0,
                             K = 2, p = 1, ker = "dnorm", maxiter = 50,
                             tol = 1e-6 ) {
#
# Finds a bootstrap estimate of a confidence interval at a significance level alpha
# for the estimated slope for the local polynomial estimate of the psychometric function
# with guessing and lapsing rates. The confidence interval is based on bootstrap percentiles.
#
# See Efron & Tibshirani's "An introduction to the bootstrap", 1993
#
# INPUT
#
# TH   - required threshold level
# r    - number of successes at points x
# m    - number of trials at points x
# x    - stimulus levels
# N    - number of bootstrap replications; N should be at least 1000 for reliable results
# h0   - bandwidth
#
# OPTIONAL INPUT
#
# alpha    - significance level of the confidence interval; default is 0.05
# X        - set of values at which estimates of the psychometric function
# 		for the slope estimation are to be obtained; if not given, 1000
# 		equally spaced points from minimum to maximum of 'x' are used
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
# ci  - confidence interval based on bootstrap percentiles
# sl0 - slope estimate

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
    if( length(alpha)>1 ){
    	  stop('Significance level must be scalar');
    	}

    if( alpha <= 0 || alpha > 0.5 ) {
        stop('Significance level must be between 0 and 0.5');
    }
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

    if( N < 1000 ) {
        warning( "number of bootstrap should be larger than 1000\n otherwise results might be unreliable" );
    }


    n <- length( x );

# INITIAL ESTIMATE
# initial estimates with bandiwdth h0
options(warn=-1)
    f <- locglmfit( x, r, m, x, h0, FALSE, link, guessing, lapsing, K,
                    p, ker, maxiter, tol )$pfit;

# dense version for estimation of the threshold
    F <- locglmfit( X, r, m, x, h0, FALSE, link, guessing, lapsing, K,
                    p, ker, maxiter, tol )$pfit;

# THRESHOLD ESTIMATE
		ci <- NULL
        sl0 <- threshold_slope( F, X, TH )$slope;

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
    sl_boot <- NULL;

# BOOTSTRAP ESTIMATES OF THE THRESHOLD

    for( i in 1:N ) {
    	ftmp <- locglmfit( X, samp[,i], m, x, h0, FALSE, link, guessing,
                           lapsing, K, p, ker, maxiter, tol )$pfit;
        sl_boot[i] <- threshold_slope( ftmp, X, TH )$slope;
     }

options(warn=0)

    ci[1] <- quantile( sl_boot, probs = alpha / 2 );
    ci[2] <- quantile( sl_boot, probs = 1 - alpha / 2 );

	value <- NULL
    value$ci <- ci
    value$sl0 <- sl0

    return( value );
}
