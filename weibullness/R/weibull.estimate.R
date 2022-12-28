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
     upper = sum( (x^lower)*log(x) ) / sum( x^lower ) - meanlog
     interval = c(lower,1/upper)
  }
  EEweibull = function(alpha,x) {
     xalpha = x^alpha
     sum(log(x)*(xalpha)) / sum(xalpha) - 1/alpha - mean(log(x))
  }
  tmp = uniroot(EEweibull, interval=interval, x=x,tol=tol,maxiter=maxiter,trace=trace)
  alpha = tmp$root   
  beta  = mean(x^alpha)^(1/alpha)
  structure(list(shape=alpha,scale=beta,threshold=threshold),class="weibull.estimate"  )
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
   if  (missing(interval.threshold)) {
      SD = sd(x)
      TINY = SD * .Machine$double.neg.eps^0.25 
      LOWER = min(x) - SD
      UPPER = min(x) - TINY 
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
#######################
