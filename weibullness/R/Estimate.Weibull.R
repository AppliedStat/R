# Based on Farnum/Booth:1997
weibull.mle <- 
function(x, threshold, interval, interval.threshold, extendInt="downX",
         a, tol=.Machine$double.eps^0.25, maxiter=1000, trace=0) 
{
  if (missing(threshold)) { 
     threshold = weibull.threshold(x, a, interval.threshold, extendInt)
  }
  x = x - threshold

  ## TINY = .Machine$double.neg.eps
  ## if ( any (x < TINY) ) stop("The data should be positive")
  if (missing(interval)) {
     meanlog = mean(log(x))
     lower = 1 / ( log(max(x)) - meanlog )
     upper.rev = sum( (x^lower)*log(x) ) / sum( x^lower ) - meanlog
     upper = 1 / upper.rev
     if (is.nan(upper)) upper=.Machine$double.xmax^0.2
     interval = c(lower,upper)
  }
  EEweibull = function(alpha,x) {
     TINY = .Machine$double.neg.eps
     xalpha = x^alpha
     if ( sum(xalpha) < TINY ) return( mean(log(x)) - 1/alpha - mean(log(x)) )
     sum(log(x)*(xalpha)) / sum(xalpha) - 1/alpha - mean(log(x))
  }
  tmp = uniroot(EEweibull, interval=interval, x=x, extendInt="upX", tol=tol,maxiter=maxiter,trace=trace)
  alpha = tmp$root   
  beta  = mean(x^alpha)^(1/alpha)
  structure(list(shape=alpha,scale=beta,threshold=threshold),class="weibull.estimate")
}
#------------------------------
# We need n when the data are right-censored at the max. obs
weibull.wp <- function(x,n, a){   ## It was a=0.5
  x = sort(x)
  r = length(x)
  if ( missing(n) ) n=r
  if (missing(a)) { a = ifelse(n <= 10, 3/8, 1/2) }
  if ( n < r ) stop("n cannot be smaller than the number of observations")
  R = 1-ppoints(n,a=a)
  LM=lm(log(-log( R[1:r] ))~log(x))
  B = as.numeric( coef (LM) )
  alpha = B[2]
  beta  = exp(-B[1]/B[2])
  structure(list(shape=alpha,scale=beta), class="weibull.estimate")
}
#------------------------------
weibull.rm <- function(x, a){   ## Repeated median (added on 2022-07-10)
  x = sort( log(unique(x)) )
  n = length(x)
  if (missing(a)) { a = ifelse(n <= 10, 3/8, 1/2) }
  y = log(-log(1-ppoints(n, a=a)))
  DX = outer(x,x,"-"); DY = outer(y,y,"-")
  diag(DX) = NA
  alpha.hat = median( apply(DY/DX,2, median, na.rm=TRUE) )  # slope 
  intercept = median(y-alpha.hat*x) # This intercept is better than the below. 
  ## intercept = median( apply((outer(x,y,"*")-outer(y,x,"*"))/DX, 2, median, na.rm=TRUE) )
  beta.hat = exp(-intercept/alpha.hat)
  structure( list(shape=alpha.hat, scale=beta.hat),  class = "weibull.estimate")
}
#------------------------------

#------------------------------
weibull.threshold <- 
function(x, a, interval.threshold, extendInt="downX") {
   n = length(x)
   if (missing(a)) { a = ifelse(n <= 10, 3/8, 1/2) }
   EE.weibull.threshold = function(threshold) {
      x = sort(x)
      v = log(-log(1-ppoints(n,a=a)))
      u1= log(x-threshold)
      w = -1/(x-threshold)

      u0 = u1 - mean(u1)
      v0 = v - mean(v)
      w0 = w - mean(w)
      return( sum(w0*v0)/sum(u0*v0) - sum(u0*w0)/sum(u0*u0) )
   }
   minx = min(x)
   TINY = .Machine$double.neg.eps^0.5 
   SD = mad(x)
   if  (missing(interval.threshold)) {
      LOWER = minx - TINY - SD
      UPPER = minx - TINY 
      interval.threshold = c(LOWER,UPPER) 
   }
   tmp = uniroot(EE.weibull.threshold, 
                 interval=interval.threshold, extendInt=extendInt)
   ans = tmp$root
   names(ans) = "threshold"
   return(ans)
}
#------------------------------

#------------------------------
print.weibull.estimate <-
function(x, digits = getOption("digits"), ...)
{ 
   ans = format(x, digits=digits) 
   dn  = dimnames(ans)
   print(ans, quote=FALSE)
   invisible(x)
}

###################################################################
# Interval Censored 
###################################################################

# =================================================================
# Incomplete Upper Gamma Function 
#------------------------------------------------------------------
lGAMMA.incomplete = function(a,b) {
    if( a <= 0 ) a = .Machine$double.neg.eps
    return( pgamma(b,a,lower.tail=FALSE,log.p=TRUE) + lgamma(a) )
}
GAMMA.incomplete = function(a,b) { exp(lGAMMA.incomplete(a,b)) }
# =================================================================

# ===================================================================
# Version   : 1.0,   Dec. 25, 2016
#             1.1,   Feb. 14, 2020
#             1.2,   Jan.  4, 2023
# NOTE: Euler-Mascheroni constant gam = -digamma(1)
# ===================================================================
U.i.s = function(kappa.s, theta.s, ai, bi) { 
   TINY = .Machine$double.neg.eps

   if ( ai > bi ) return(0.0)
   if ( abs(bi-ai) < TINY ) return(log(ai))

   t.ai = (ai/theta.s)^kappa.s
   exp.t.ai = exp(-t.ai)

   t.bi = (bi/theta.s)^kappa.s
   exp.t.bi = exp(-t.bi)

   D.i.s = exp.t.ai - exp.t.bi

   if (ai>0) {
      tmpai = log(ai)*exp.t.ai + GAMMA.incomplete(0,t.ai)/kappa.s
   } else {
      tmpai = log(theta.s) + digamma(1)/kappa.s
   }

   BIG = .Machine$double.xmax^0.2
   if (bi >= BIG) {   ## .Machine$double.xmax
      tmpbi = 0
   } else {
      tmpbi = log(bi)*exp.t.bi + GAMMA.incomplete(0,t.bi)/kappa.s
   }

   return( (tmpai-tmpbi)/D.i.s )
}
#
U.i.s.0.inf = function(kappa.s, theta.s) { log(theta.s) + digamma(1)/kappa.s }
#
# U.i.s (2, 3, 0, Inf);  U.i.s.0.inf(2, 3)
# U.i.s (2, 3, 3, 3);    U.i.s (2, 3, 3, 3.000)
U.i.s.0.bi = function(kappa.s, theta.s, bi) { 
   A = log(theta.s) + digamma(1)/kappa.s
   t.bi = (bi/theta.s)^kappa.s
   tmp = log(bi)*exp(-t.bi) + GAMMA.incomplete(0, t.bi)/kappa.s
   (A - tmp) / (1- exp(-t.bi))
   ## A  /(1- exp(-t.bi)) - tmp/exp(-t.bi) 
}
# U.i.s (2,3, 0.0, 3) ;    U.i.s.0.bi(2,3, 3)
#
U.i.s.ai.inf = function(kappa.s, theta.s, ai) { 
   t.ai = (ai/theta.s)^kappa.s
   tmp = log(ai)*exp(-t.ai) + GAMMA.incomplete(0, t.ai)/kappa.s
   tmp / exp(-t.ai)
}
# U.i.s (2, 3, 1, Inf) ;    U.i.s.ai.inf(2,3,1)


# ===================================================================
# Version   : 1.0,   Dec. 25, 2016
#             1.1,   Feb. 14, 2020
#             1.2,   Jan.  4, 2023
# NOTE: Euler-Mascheroni constant gam = -digamma(1)
# ===================================================================
V.i.s = function(kappa, kappa.s, theta.s, ai, bi) { 
   TINY = .Machine$double.neg.eps
   if ( ai > bi ) return(0.0)
   if ( abs(bi-ai) < TINY ) return(ai^kappa)

   t.ai = (ai/theta.s)^kappa.s
   t.bi = (bi/theta.s)^kappa.s

   D.i.s = exp(-t.ai) - exp(-t.bi)

   tmpai = GAMMA.incomplete( (kappa+kappa.s)/kappa.s, (ai/theta.s)^kappa.s )
   tmpbi = GAMMA.incomplete( (kappa+kappa.s)/kappa.s, (bi/theta.s)^kappa.s )

   theta.s^kappa * (tmpai-tmpbi)/D.i.s
}
# V.i.s (1, 2, 3, 1, Inf);  V.i.s (1, 2, 3, 0, Inf);  V.i.s (1, 2, 3, 0, 4)
# V.i.s (1, 2, 3, 9, 9);   V.i.s (1, 2, 3, 8.9999, 9)
# V.i.s (1, 2, 3, 1, 1);   V.i.s (1, 2, 3, 0.9999, 1)

##################################################################
## V.i.s = function(kappa, kappa.s, theta.s, ai, bi) 
## U.i.s = function(kappa.s, theta.s, ai, bi) 
weibull.ic = function(X, start=c(1,1), maxits=10000, eps=1E-5){
   kappa = start[1]
   theta = start[2]
   ij = dim(X)
   n  = ij[1]
   if ( ij[2] > 2 ) stop(" Warning: The data X should be n x 2 matrix.")
   if ( any(X[,1] > X[,2]) )  stop(" Warning: The data should satisfy a <= b.");
   iter = 0
   converged = FALSE

   # Start the EM
   colnames(X) = NULL; rownames(X) = NULL

   ai = X[,1];  bi = X[,2]

   fn = function(newkappa) {
        sumU.i.s = 0
        sumV.i.s = 0
        for ( i in seq_len(n) ) sumU.i.s = sumU.i.s + U.i.s(kappa, theta, ai[i], bi[i])
        for ( i in seq_len(n) ) sumV.i.s = sumV.i.s + V.i.s(newkappa, kappa, theta, ai[i], bi[i])
        -1 * ( n*log(newkappa) + newkappa*sumU.i.s - n*log(sumV.i.s) )
   }

   while( (iter<maxits)&(!converged) ) {
      OPT = nlm(fn, kappa); newkappa = OPT$estimate
      ## OPT = nlminb(kappa,fn); newkappa = OPT$par

      meanV.i.s = 0
      for ( i in seq_len(n) ) meanV.i.s = meanV.i.s + V.i.s(newkappa, kappa, theta, ai[i], bi[i])/n
      newtheta = meanV.i.s^(1/newkappa)

      # assess convergence
      converged = (abs(newkappa-kappa)<eps*abs(newkappa)) & (abs(newtheta -theta)<eps*abs(newtheta))
      iter = iter+1
      kappa   = newkappa
      theta   = newtheta
   }
   list( shape=newkappa, scale=newtheta, iter=iter, conv=converged )
}
###################################################################
