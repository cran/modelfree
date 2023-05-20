#' Local polynomial estimator of a psychometric function
#'
#' Local polynomial estimator for the psychometric function and eta function (psychometric function
#' transformed by link) for binomial data; also returns the hat matrix H.
#'
#' @usage locglmfit( xfit, r, m, x, h, returnH = FALSE, link = "logit",
#'                guessing = 0, lapsing = 0, K = 2, p = 1,
#'                ker = "dnorm", maxiter = 50, tol = 1e-6 )
#'
#' @param xfit points at which to calculate the estimate pfit
#' @param r    number of successes at points x
#' @param m    number of trials at points x
#' @param x    stimulus levels
#' @param h    bandwidth(s)
#' @param returnH  (optional) logical, if TRUE then hat matrix is calculated; default is FALSE
#' @param link     (optional) name of the link function; default is 'logit'
#' @param guessing (optional) guessing rate; default is 0
#' @param lapsing  (optional) lapsing rate; default is 0
#' @param K    (optional) power parameter for Weibull and reverse Weibull link; default is 2
#' @param p        (optional) degree of the polynomial; default is 1
#' @param ker      (optional) kernel function for weights; default is 'dnorm'
#' @param maxiter  (optional) maximum number of iterations in Fisher scoring; default is 50
#' @param tol      (optional) tolerance level at which to stop Fisher scoring; default is 1e-6
#'
#' @returns \verb{pfit    }    value of the local polynomial estimate at points xfit
#' @returns \verb{etafit  }  estimate of eta (link of pfit)
#' @returns \verb{H       }       hat matrix (OPTIONAL)
#'
#' @examples
#' data("Miranda_Henson")
#' x = Miranda_Henson$x
#' r = Miranda_Henson$r
#' m = Miranda_Henson$m
#' numxfit <- 199; # Number of new points to be generated minus 1
#' xfit <- (max(x)-min(x)) * (0:numxfit) / numxfit + min(x)
#' # Find a plug-in bandwidth
#' bwd <- bandwidth_plugin( r, m, x)
#' pfit <- locglmfit( xfit, r, m, x, bwd )$pfit
#' # Plot the fitted curve
#' plot( x, r / m, xlim = c( 0.1, 1.302 ), ylim = c( 0.0165, 0.965 ), type = "p", pch="*" )
#' lines(xfit, pfit )
#'
#' @export
locglmfit<-function( xfit, r, m, x, h, returnH = FALSE, link = "logit",
                     guessing = 0, lapsing = 0, K = 2, p = 1,
                     ker = "dnorm", maxiter = 50, tol = 1e-6 ) {
#
# Local polynomial estimator for the psychometric function and eta function (psychometric function
# transformed by link) for binomial data; also returns the hat matrix H. Actual calculations are
# done in LOCGLMFIT_PRIVATE or LOCGLMFIT_SPARSE_PRIVATE depending on the size of the data set.
# Here, the data are split into several parts to speed up the calculations.
#
#
#INPUT
#
# xfit - points at which to calculate the estimate pfit
# r    - number of successes at points x
# m    - number of trials at points x
# x    - stimulus levels
# h    - bandwidth(s)
#
# OPTIONAL INPUT
#
# returnH  - logical, if TRUE then hat matrix is calculated; default is FALSE
# link     - name of the link function; default is 'logit'
# guessing - guessing rate; default is 0
# lapsing  - lapsing rate; default is 0
# K    - power parameter for Weibull and reverse Weibull link; default is 2
# p        - degree of the polynomial; default is 1
# ker      - kernel function for weights; default is 'dnorm'
# maxiter  - maximum number of iterations in Fisher scoring; default is 50
# tol      - tolerance level at which to stop Fisher scoring; default is 1e-6
#
# OUTPUT
#
# pfit    - value of the local polynomial estimate at points xfit
# etafit  - estimate of eta (link of pfit)
# H       - hat matrix (OPTIONAL)

# First 5 arguments are mandatory
    if( missing("xfit") || missing("r") || missing("m") || missing("x") || missing("h") ) {
        stop("Check input. First 5 arguments are mandatory");
    }

# CHECK ROBUSTNESS OF INPUT PARAMETERS
    if ( !is.vector( xfit ) ) {
        stop("Vector xfit should be a column vector");
    }
    checkdata<-list();
    checkdata[[1]] <- x;
    checkdata[[2]] <- r;
    checkdata[[3]] <- m;
    checkinput("checkdesignpoints",x)
    checkinput( "psychometricdata", checkdata );
    rm( checkdata )
    if ( !is.vector( h ) ) {
        stop("Bandwidths h should be a vector");
    }
    if( !( length( h ) == 1 || length( h ) == length( xfit ) ) ) {
        stop( "Bandwidth h must be either a scalar or a vector with the same number of elements as xfit" );
    }
    checkinput( "linkfunction", link );
    if( length( guessing ) > 1 ) {
        stop( "Guessing rate must be a scalar" );
    }
    if( length( lapsing ) > 1 ) {
        stop( "Lapsing rate must be a scalar" );
    }

    checkinput( 'guessingandlapsing', c(guessing, lapsing) );
   if (link == "weibull" || link == "revweibull"){
	    checkinput( "exponentk", K );
	    }
    pn <- list()
	pn[[1]] <- p
	pn[[2]] <- x
    checkinput( "degreepolynomial", pn );
    checkinput( 'kernel', ker );
    checkinput( 'maxiter', maxiter );
    checkinput( 'tolerance', tol );

# INITIAL VALUES
    split <- 20;

    Lxfit <- length(xfit);
    Lx <- length(x);

	 value <- NULL
    pfit <- NULL;
    etafit    <- NULL;
    if( returnH  ) H <- NULL;

	if( link == "logit"      ||
        link == "probit"     ||
        link == "loglog"     ||
        link == "comploglog" ||
        link == "weibull"    ||
        link == "revweibull" ) {

    			link <- paste( link, "_link_private", sep = "" );
	}

    if( Lx > 15 ) {
# big data
# First try to load package SparseM
      if (requireNamespace("SparseM", quietly = TRUE)) {
        fun_estim <- locglmfit_sparse_private;
      } else {
        fun_estim <- locglmfit_private;
        message("Package SparseM not installed. No sparse matrices used");
        message("SparseM can be found at CRAN web site http://cran.r-project.org/");
      }
 #       options(warn = -1);
  #      existSparseM <- library( SparseM, logical.return = TRUE );
   #     options(warn = 0);
    #    if( existSparseM ) {
     #       fun_estim <- locglmfit_sparse_private;
      #  }
       # else {
        #    fun_estim <- locglmfit_private;
         #   message("Package SparseM not installed. No sparse matrices used");
          #  message("SparseM can be found at CRAN web site http://cran.r-project.org/");
        #}
    }
    else {
# small data
        fun_estim <- locglmfit_private;
    }

# SPLIT AND EVALUATION
################################################################ SCALAR h
    if( length( h ) == 1 ) {
# with Hat matrix
        if( returnH ) {
            if( Lxfit <= split ) {
# small x

                value <- fun_estim( xfit, r, m, x, h, returnH, link, guessing,
                       lapsing, K, p, ker, maxiter, tol);
            }
            else {
# large x
# number of parts into which the fitting is divided
                fLx = floor( Lxfit / split );
# initialise output
                    for( i in 0:(fLx-1) ) {
# part of the fit
                        value1 <- fun_estim( xfit[i*split+c(1:split)], r, m, x, h,
                                returnH, link, guessing, lapsing, K, p, ker,
                                maxiter, tol );
# put the fits together
                        pfit <- c( pfit, value1$pfit );
                        etafit    <- c( etafit,    value1$etafit );
                        H      <- rbind(H,   value1$H );
                    }
# final part of the fit
                    if( ( split * fLx ) < Lxfit ) {
                        value1 <- fun_estim( xfit[c(1+split*fLx):Lxfit], r, m, x, h,
                                returnH, link, guessing, lapsing, K, p, ker,
                                maxiter, tol );

# put the fits together
                        pfit <- c( pfit, value1$pfit );
                        etafit    <- c( etafit,    value1$etafit );
                        H      <- rbind(H,   value1$H );
                    }
#   values to return

value$pfit <- pfit
value$etafit <- etafit
value$H <- H

				}
        }
        else { # no Hat matrix

            if( Lxfit <= split ) {
# small x

                value <- fun_estim( xfit, r, m, x, h,returnH, link, guessing,
                       lapsing, K, p, ker, maxiter, tol );
            }
            else {
# large x
# number of parts into which the fitting is divided
                fLx = floor( Lxfit / split );

# initialise output
                for( i in 0:(fLx-1) ) {
# part of the fit

                    value1 <- fun_estim(  xfit[i*split + c(1:split)], r, m, x, h,
                                returnH, link, guessing, lapsing, K, p, ker,
                                maxiter, tol );
# put the fits together
                    pfit <- c( pfit, value1$pfit );
                    etafit    <- c( etafit,    value1$etafit );
                }
# final part of the fit
                if( ( split * fLx ) < Lxfit ) {
                    value1 <- fun_estim( xfit[c(1+split*fLx):Lxfit], r, m, x, h,
                            returnH, link, guessing, lapsing, K, p, ker,
                            maxiter, tol );
# put the fits together
                    pfit <- c( pfit, value1$pfit );
                    etafit    <- c( etafit,    value1$etafit );
                }

#   values to return

value$pfit <- pfit
value$etafit <- etafit

			}
        } # if( returnH )
    }

 ################################################################ VECTOR h
    else { # if( length( h ) == 1 )

# with Hat matrix
        if( returnH ) {
            if ( Lxfit <= split ) {
# small x
                value <- fun_estim( xfit, r, m, x, h, returnH, link, guessing,
                       lapsing, K, p, ker, maxiter, tol );

            }
            else {
# large x
# number of parts into which the fitting is divided
                fLx = floor( Lxfit / split );
                    for( i in 0:(fLx-1) ) {
# part of the fit
                        value1 <- fun_estim( xfit[i*split + c(1:split)], r, m, x,
                                h[i*split + c(1:split)], returnH, link,
                                guessing, lapsing, K, p, ker, maxiter,
                                tol );
# put the fits together
                        pfit <- c( pfit, value1$pfit );
                        etafit    <- c( etafit,    value1$etafit );
                        H      <- rbind(H,   value1$H );
                    }
# final part of the fit
                    if( ( split * fLx ) < Lxfit ) {
                        value1 <- fun_estim( xfit[c(1+split*fLx):Lxfit], r, m, x,
                                h[c(1+split*fLx):Lxfit], returnH, link,
                                guessing, lapsing, K, p, ker, maxiter,
                                tol );
# put the fits together
                        pfit <- c( pfit, value1$pfit );
                        etafit    <- c( etafit,    value1$etafit );
                        H      <- rbind(H,   value1$H );
                    }
#   values to return

value$pfit <- pfit
value$etafit <- etafit
value$H <- H

				}
        }# end if(retunH)
        else { # no Hat matrix
            if( Lxfit <= split ) {
# small x
                value <- fun_estim( xfit, r, m, x, h, returnH, link, guessing,
                       lapsing, K, p, ker, maxiter, tol );
            }
            else {
# large x
# number of parts into which the fitting is divided
            fLx = floor( Lxfit / split );

# initialise output
                for( i in 0:(fLx-1) ) {
# part of the fit
                    value1 <- fun_estim( xfit[i*split + c(1:split)], r, m, x,
                            h[i*split + c(1:split)], returnH, link,
                            guessing, lapsing, K, p, ker, maxiter, tol );
# put the fits together
                    pfit <- c( pfit, value1$pfit );
                    etafit    <- c( etafit,    value1$etafit );
                }
# final part of the fit
                if( ( split * fLx ) < Lxfit ) {
                    value1 <- fun_estim( xfit[c(1+split*fLx):Lxfit], r, m, x,
                            h[c(1+split*fLx):Lxfit], returnH, link,
                            guessing, lapsing, K, p, ker, maxiter, tol );
# put the fits together
                    pfit <- c( pfit, value1$pfit );
                    etafit    <- c( etafit,    value1$etafit );
                }

#   values to return

value$pfit <- pfit
value$etafit <- etafit

			}
        } # if( returnH )

    } # if( length( h ) == 1 )

   return( value )
}
