globalVariables( c("n.times.eBias.of.mad","n.times.eBias.of.shamos",
                   "n.times.eVar.of.median", 
                   "n.times.eVar.of.HL1",
                   "n.times.eVar.of.HL2",
                   "n.times.eVar.of.HL3",
                   "n.times.eVar.of.mad", "n.times.eVar.of.shamos") )
## load("EmpiricalResults.rda")
#============================================
# Finite Breakdown 
#============================================
finite.breakdown <-
function(n, method=c("mean","median","HL1","HL2","HL3","sd","range","mad","shamos")) { 
  n = as.integer(n)
  if ( n <= 1L ) return(0)

  method = match.arg(method)
  if ( method %in% c("mean","sd", "range") ) return(0.0)


  switch (method,
     HL1 = { tmp = floor( (n*(n-1)-2)/4 )
               m = floor( ((2*n-1) - sqrt((2*n-1)^2-8*tmp))/2 ) },
     HL2 = { tmp = floor( (n*(n+1)-2)/4 )
               m = floor( ((2*n+1) - sqrt((2*n+1)^2-8*tmp))/2 ) },
     HL3 =  {tmp = floor(0.5*(n^2-1))
               m = floor( n - sqrt(n^2-tmp)) },
     median={  m = floor(0.5*(n-1)) },
     mad   ={  m = floor(0.5*(n-1)) },
     shamos={tmp = floor( (n*(n-1)-2)/4 )
               m = floor( ((2*n-1) - sqrt((2*n-1)^2-8*tmp))/2 ) } )
  return( m/n )
}
#============================================

#============================================
# Hodges-Lehmann estimators
# HL1: median over i < j
# HL2: median over i <= j
# HL3: median over all (i,j)
#============================================
HL <- function(x, method=c("HL1","HL2","HL3"),na.rm=FALSE) {
   method = match.arg(method)
   if (na.rm) x <- x[!is.na(x)]
   xx = outer(x,x, "+")

   HL.estimation =  switch (method,
      HL1 = 0.5*median(xx[lower.tri(xx,diag=FALSE)]),
      HL2 = 0.5*median(xx[lower.tri(xx,diag=TRUE)]),
      HL3 = 0.5*median(xx) 
   )

   return( HL.estimation )
}
#============================================

#============================================
# Shamos (i < j) 
#============================================
shamos <- function(x,constant=1.048358, na.rm=FALSE, IncludeEqual=FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  w1 = outer(x,x, "-")
  w2 = abs( w1[lower.tri(w1,diag=IncludeEqual)]  )
  constant * median( w2 )
}
#============================================

#============================================
# Unbiasing factors for sigma 
#           based on the empirical/estimated biases  
#============================================
c4.factor = function(n, method=c("sd","range","mad","shamos")) {
  n = as.integer(n)
  if ( n <= 1L ) return(NA)
  
  method = match.arg(method)

  if ( n <= 100L ) {
      switch(method,
             sd = return(sqrt(2/(n-1))*exp(lgamma(n/2)-lgamma((n-1)/2))),
          range = { tmp = integrate ( 
                          function(x){Phi=pnorm(x);(1-(1-Phi)^n-Phi^n)}, 
                          lower=0,upper=Inf); return(2*tmp$value)}, 
            mad = return(1+n.times.eBias.of.mad[n]/n),
         shamos = return(1+n.times.eBias.of.shamos[n]/n) )
  } else {
      switch(method,
             sd = return(sqrt(2/(n-1))*exp(lgamma(n/2)-lgamma((n-1)/2))),
          range = { tmp = integrate ( 
                          function(x){Phi=pnorm(x);(1-(1-Phi)^n-Phi^n)}, 
                          lower=0,upper=Inf); return(2*tmp$value)}, 
            mad = return(1 - 0.76213/n - 0.86413/n^2),
         shamos = return(1 + 0.414253297/n + 0.442396799/n^2) )
  }
}
#-----
# Unbiased mad
mad.unbiased = function(x,center=median(x),constant=1.4826,na.rm=FALSE){
  if (na.rm) x = x[!is.na(x)]
  n = length(x)
  if (n == 1) return(0)
  constant*median(abs(x-center))/c4.factor(n,"mad") 
}
#-----
# Unbiased Shamos 
shamos.unbiased = function(x,constant=1.048358,na.rm=FALSE,IncludeEqual=FALSE){
  if (na.rm) x = x[!is.na(x)]
  n = length(x)
  if (n == 1) return(0)
  w1 = outer(x,x, "-")
  w2 = abs( w1[lower.tri(w1,diag=IncludeEqual)]  )
  constant * median( w2 ) / c4.factor(n,"shamos") 
}
#============================================

#============================================
#
# RE under N(0,1)
#
RE <-
function(n, method=c("mean", "median","HL1","HL2","HL3", 
        "sd", "range", "mad","shamos"), correction=FALSE) {
  n = as.integer(n)
  if ( n <= 0L ) return(NA)

  method = match.arg(method)
  if (method %in% c("mean", "sd")) return(1.0) 

  if ( n == 1L ) {
     switch (method, 
          HL1 = return(NA),
        range = return(NA),
          mad = return(NA),
       shamos = return(NA) )
  }

     mu = c4.factor(n, "range")
#-----------------------------------
  d3 = function(n) {
       f= function(x,y) {
          (y^2) * exp(-(x^2+(x+y)^2)/2) * (pnorm(x+y)-pnorm(x))^(n-2)
       }
       tmp = integrate(function(y) {
             sapply(y,function(y) {
                  integrate(function(x){f(x,y)},-Inf,Inf)$value})}, 0, Inf)
       sqrt( n*(n-1)/2/pi*tmp$value - c4.factor(n, "range")^2 )
  }
  #------------------
  if ( n <= 100L ) {
    re = switch(method,
           median = (1/n.times.eVar.of.median[n]), 
           HL1    = (1/n.times.eVar.of.HL1[n]), 
           HL2    = (1/n.times.eVar.of.HL2[n]), 
           HL3    = (1/n.times.eVar.of.HL3[n]), 
           range  = ((1-c4.factor(n,"sd")^2) / d3(n)^2 ), 
           mad    = ((1-c4.factor(n,"sd")^2)/(n.times.eVar.of.mad[n]/n)),
           shamos = ((1-c4.factor(n,"sd")^2)/(n.times.eVar.of.shamos[n]/n)) )
  } else {
    re = switch(method,
           median ={tmp=ifelse(n%%2L==1L,-0.6589-0.9430/n,-2.195+1.929/n);
                   (1/(1.57+tmp/n))},
           HL1    = (1/(1.0472 + 0.1127/n + 0.8365/n^2)), 
           HL2    = (1/(1.0472 + 0.2923/n + 0.2258/n^2)), 
           HL3    = (1/(1.0472 + 0.2022/n + 0.4343/n^2)), 
           range  = ((1-c4.factor(n,"sd")^2) / d3(n)^2 ), 
           mad    = {tmp=ifelse(n%%2L==1L,0.2996-149.357/n,-2.417-153.01/n);
                     (1/(2.7027+tmp/n))},
           shamos = (1/(1.15875+2.822/n+12.238/n^2)) )
  }
  if (method %in% c("median","HL1", "HL2", "HL3")) return(re)
  return(ifelse(correction, re*(c4.factor(n,method)/c4.factor(n,"sd"))^2,re))
}
#============================================

#
#============================================
# eVar under N(0,1)
#============================================
# 
evar = function(n, method=c("median","HL1","HL2","HL3","mad","shamos")) {
  n = as.integer(n)
  if ( n <= 0L ) return(NA)
  if ( n == 1L && method == "HL1" ) return(NA)

  method = match.arg(method)

  if ( n <= 100L ) {
      switch(method,
             median = return(n.times.eVar.of.median[n]/n), 
             HL1    = return(n.times.eVar.of.HL1[n]/n),
             HL2    = return(n.times.eVar.of.HL2[n]/n),
             HL3    = return(n.times.eVar.of.HL3[n]/n),
             mad    = return(n.times.eVar.of.mad[n]/n),
             shamos = return(n.times.eVar.of.shamos[n]/n) )
  } else {
      switch(method,
         median = return((1/RE(n,method="median"))/n),
         HL1    = return((1/RE(n,method="HL1"))/n),
         HL2    = return((1/RE(n,method="HL2"))/n),
         HL3    = return((1/RE(n,method="HL3"))/n),
         mad    = return((1/RE(n,method="mad"))*(1-c4.factor(n,"sd")^2)),
         shamos = return((1/RE(n,method="shamos"))*(1-c4.factor(n,"sd")^2))    )
  } 
}

#============================================
# Unbiasing factors for sigma^2 
#           based on the empirical/estimated biases and variances
#============================================
w4.factor =  function(n, method=c("mad2","shamos2")) {
   method = match.arg(method)
   switch(method,
          mad2 = return(evar(n, method="mad") + c4.factor(n,"mad")^2),
       shamos2 = return(evar(n, method="shamos") + c4.factor(n,"shamos")^2)
   )
}
#-----
# Unbiased mad^2
mad2.unbiased = function(x,center=median(x),constant=1.4826,na.rm=FALSE){
  if (na.rm) x = x[!is.na(x)]
  n = length(x)
  if (n == 1) return(0)
  mad(x, center, constant, na.rm)^2 / w4.factor(n,"mad2") 
}
#-----
# Unbiased Shamos^2
shamos2.unbiased = function(x,constant=1.048358,na.rm=FALSE,IncludeEqual=FALSE){
  if (na.rm) x = x[!is.na(x)]
  n = length(x)
  if (n == 1) return(0)

  shamos(x, constant, na.rm, IncludeEqual)^2 / w4.factor(n,"shamos2")
}
#============================================
# The American Standard uses "3*sigma"
# The British Standard  uses "3.09*sigma"
# Control Limits Factor for the sample standard deviation (s). 
A3.factor = function(n, method=c("sd","range","mad","shamos"), sigma.factor=3){
  n = as.integer(n)
  if ( n <= 1L ) return(NA)
  method = match.arg(method)
  switch(method,
         sd = return(sigma.factor/c4.factor(n,"sd")/sqrt(n)),
      range = return(sigma.factor/c4.factor(n,"range")/sqrt(n)),
        mad = return(sigma.factor/c4.factor(n,"mad")/sqrt(n)),
     shamos = return(sigma.factor/c4.factor(n,"shamos")/sqrt(n)) )
}
#============================================

