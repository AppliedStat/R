#=======================================================================
#  Function   : Laplace (Double Exponential) Distribution
#-----------------------------------------------------------------------
#
#  Programmer : Park, Chanseok
#  Date       : Mar. 15, 1999
#  
#  Usage      : dlaplace(x, location=0, rate=1) 
#               plaplace(q, location=0, rate=1)
#               qlaplace(p, location=0, rate=1)
#               rlaplace(n, location=0, rate=1)
#   CF :  help.findsum (".Data")
#=======================================================================
dlaplace <-
function(x, location=0, scale=1, rate=1/scale, log=FALSE) {
   k <- max(lx<-length(x),lloc<-length(location),lrate<-length(rate))
   if (lx   < k)        x = rep(x, length=k) 
   if (lloc < k) location = rep(location, length=k) 
   if (lrate< k)     rate = rep(rate, length=k) 
   if (log == FALSE) {
       y = 0.5*rate*exp(-abs(x-location)*rate);
   } else {
       y = log(0.5) + log(rate) - abs(x-location)*rate 
   }
   if(!is.null(Names <- names(x))) names(y) = rep(Names,length=length(y))
   y
}
##----------
plaplace <-
function(q, location=0, scale=1, rate=1/scale, lower.tail = TRUE, log.p=FALSE) {
   k <- max(lx<-length(q),lloc<-length(location),lrate<-length(rate))
   if (lx   < k)        q = rep(q, length=k) 
   if (lloc < k) location = rep(location, length=k) 
   if (lrate< k)     rate = rep(rate, length=k) 
   y = 0.5*exp(-abs(q-location)*rate);
   y = 0.5 - sign(q-location)*(y-0.5);

   if(!is.null(Names<-names(q))) names(y) = rep(Names,length=length(y))

   if (lower.tail == FALSE) y = 1-y
   if (log.p == FALSE) {
      y
   } else {
     log(y) 
   }
}
##----------
qlaplace <- 
function(p, location=0, scale=1, rate=1/scale, lower.tail = TRUE, log.p = FALSE) {
   k <- max(lp<-length(p),lloc<-length(location),lrate<-length(rate))
   if (lp   < k)        p = rep(p, length=k)
   if (lloc < k) location = rep(location, length=k)
   if (lrate< k)     rate = rep(rate, length=k)

   if (log.p == TRUE) p = exp(p)
   if (lower.tail == FALSE) p = 1-p

   index <- sign(0.5-p);
   q <- index*log(1+index*(2*p-1)) / rate + location;
   if(!is.null(Names<-names(p))) names(q) = rep(Names,length=length(q))
   q
}
##----------
rlaplace <- 
function(n, location=0, scale=1, rate=1/scale) {
   lloc<-length(location); lrate<-length(rate)
   if (lloc < n) location = rep(location, length=n)
   if (lrate< n)     rate = rep(rate, length=n)
   y = 2*runif(n)-1
   index = sign(y);
   location - index*log(1-index*(y))/rate;
}
#=======================================================================

