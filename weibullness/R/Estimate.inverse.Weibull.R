###################################################################
# Inverse Weibull
###################################################################
invweibull.mle <- function(x, interval, tol = .Machine$double.eps^0.25, maxiter = 1000, trace = 0 ) {
   y = 1/x
   if (missing(interval)) {
     meanlog = mean(log(y))
     lower = 1 / ( log(max(y)) - meanlog )
     upper.rev = sum( (y^lower)*log(y) ) / sum( y^lower ) - meanlog
     upper = 1 / upper.rev
     if (is.nan(upper)) upper=.Machine$double.xmax^0.2
     interval = c(lower,upper)
  }
  EEweibull = function(beta,y) {
     TINY = .Machine$double.neg.eps 
     ybeta = y^beta
     if ( sum(ybeta) < TINY ) return( mean(log(y)) - 1/beta - mean(log(y)) )
     sum(log(y)*(ybeta)) / sum(ybeta) - 1/beta - mean(log(y))
  }
  tmp = uniroot(EEweibull, interval=interval, y=y,tol=tol, extendInt="upX", maxiter=maxiter,trace=trace)
  beta = tmp$root
  theta  = 1 / mean(y^beta)^(1/beta)
  structure(list(shape = beta, scale = theta, threshold = 0), class = "weibull.estimate")
}

