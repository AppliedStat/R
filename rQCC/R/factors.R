globalVariables( c("n.times.eBias.of.mad","n.times.eBias.of.shamos",
                   "n.times.eVar.of.median",
                   "n.times.eVar.of.HL1",
                   "n.times.eVar.of.HL2",
                   "n.times.eVar.of.HL3",
                   "n.times.eVar.of.mad", "n.times.eVar.of.shamos") )
## load("EmpiricalResults.rda")

#============================================
# Unbiasing factors for sigma 
#           based on the empirical/estimated biases  
#============================================
c4.factor = function(n, estimator=c("sd","range","mad","shamos")) 
{
  n = as.integer(n)
  if ( n <= 1L ) return(NA)

  estimator = match.arg(estimator)

  if ( n <= 100L ) {
      switch(estimator,
             sd = return(sqrt(2/(n-1))*exp(lgamma(n/2)-lgamma((n-1)/2))),
          range = { tmp = integrate (
                          function(x){Phi=pnorm(x);(1-(1-Phi)^n-Phi^n)},
                          lower=0,upper=Inf,subdivisions=500L); 
                          return(2*tmp$value)},
            mad = return(1+n.times.eBias.of.mad[n]/n),
         shamos = return(1+n.times.eBias.of.shamos[n]/n) )
  } else {
      switch(estimator,
             sd = return(sqrt(2/(n-1))*exp(lgamma(n/2)-lgamma((n-1)/2))),
          range = { tmp = integrate (
                          function(x){Phi=pnorm(x);(1-(1-Phi)^n-Phi^n)},
                          lower=0,upper=Inf,subdivisions=500L); 
                          return(2*tmp$value)},
            mad = return(1 - 0.76213/n - 0.86413/n^2),
         shamos = return(1 + 0.414253297/n + 0.442396799/n^2) )
  }
}
#--------------------------------------------

#============================================
# Unbiasing factors for sigma^2 
#           based on the empirical/estimated biases and variances
#============================================
w4.factor =  function(n, estimator=c("mad2","shamos2")) 
{
   estimator = match.arg(estimator)
   switch(estimator,
          mad2 = return(evar(n, estimator="mad") + c4.factor(n,"mad")^2),
       shamos2 = return(evar(n, estimator="shamos") + c4.factor(n,"shamos")^2)
   )
}
#--------------------------------------------


#============================================
# Factors for constructing control charts 
factors.cc = function(n,
  factor=c("A", "A1","A2","A3","B1","B2","B3","B4","B5","B6",
           "c2","c4","d2","d3","D1","D2", "D3","D4", "E1","E2","E3"), sigmaFactor=3)
{
  n = as.integer(n)
  if ( n <= 1L ) return(NA)
  factor = match.arg(factor)

  c2fn = function(n) {sqrt(2/(n  )) * exp(lgamma(n/2) - lgamma((n-1)/2))}
  c4fn = function(n) {sqrt(2/(n-1)) * exp(lgamma(n/2) - lgamma((n-1)/2))}
  d3fn = function(n) {
     if( n <= 100L ) {
        joint2 = function(x,y) {y^2*exp(-(x^2+(x+y)^2)/2)*(pnorm(x+y)-pnorm(x))^(n-2)}
        tmp = integrate(function(y) {
        sapply(y,function(y) {
        integrate(function(x) joint2(x,y),-Inf,Inf)$value})},     0,Inf)
        sqrt( n*(n-1)/2/pi*tmp$value - c4.factor(n,estimator="range")^2 )
     } else {
       exp( 0.73784298 + 0.06390565*log(n) - 0.71491753*sqrt(log(n)) )
     }
  }

  res =  switch (factor,
      A  = {sigmaFactor/sqrt(n)},
      A1 = {sigmaFactor/c2fn(n)/sqrt(n)},
      A2 = {d2n=c4.factor(n,estimator="range"); sigmaFactor/d2n/sqrt(n)},
      A3 = {sigmaFactor/c4fn(n)/sqrt(n)},
      c2 = {c2fn(n)},
      c4 = {c4fn(n)},
      B1 = {c2n=c2fn(n); max(0,c2n-sigmaFactor*sqrt((n-1)/n-c2n^2))},
      B2 = {c2n=c2fn(n); c2n+sigmaFactor*sqrt((n-1)/n-c2n^2)},
      B3 = {c4n=c4fn(n); max(0,1-sigmaFactor/c4n*sqrt(1-c4n^2))},
      B4 = {c4n=c4fn(n);       1+sigmaFactor/c4n*sqrt(1-c4n^2)},
      B5 = {c4n=c4fn(n); max(0,c4n-sigmaFactor*sqrt(1-c4n^2))},
      B6 = {c4n=c4fn(n);       c4n+sigmaFactor*sqrt(1-c4n^2)}, 
      d2 = c4.factor(n,estimator="range"),
      d3 = d3fn(n),
      D1 = {d2n=c4.factor(n,estimator="range"); max(0,d2n-sigmaFactor*d3fn(n))},     
      D2 = {d2n=c4.factor(n,estimator="range"); d2n+sigmaFactor*d3fn(n)},
      D3 = {d2n=c4.factor(n,estimator="range"); max(0,1-sigmaFactor*d3fn(n)/d2n)},   
      D4 = {d2n=c4.factor(n,estimator="range"); 1+sigmaFactor*d3fn(n)/d2n},
      E1 = {sigmaFactor/c2fn(n)},
      E2 = {d2n=c4.factor(n,estimator="range"); sigmaFactor/d2n},
      E3 = {sigmaFactor/c4fn(n)}
  )   
  names(res) = factor
  return(res)
} 
#============================================


#============================================
# RE of the specific estimator with respect to the baseline estimator.
# Sample(s) are from N(0,1).
#
RE <-
function(n, 
  estimator=c("mean","median","HL1","HL2","HL3","sd","range","mad", "shamos"),
  poolType=c("A","B","C"), baseEstimator, basePoolType, 
  correction = TRUE, correctionBase) 
{
   estimator = match.arg(estimator)
   poolType = match.arg(poolType)

   if(missing(baseEstimator) && (estimator %in% c("mean","median","HL1","HL2","HL3"))){
       baseEstimator="mean"
   }
   if(missing(baseEstimator) && (estimator %in% c("sd","range","mad","shamos")) ) {
       baseEstimator="sd"
   }
   if(missing(basePoolType)) {basePoolType=poolType}
   if(missing(correctionBase)) {correctionBase=correction}
   NUME=evar(n,estimator=baseEstimator,poolType=basePoolType,correction=correctionBase)
   DENO=evar(n,estimator=estimator,poolType=poolType,correction=correction)
   re = NUME / DENO
   names(re) = "RE"
   if ( length(n) == 1L ) {
      message(paste0("RE of ", estimator, " with respect to ", baseEstimator))
   } else {
      message(paste0("RE of ", estimator,   " (Type=",poolType,")",
                     " with respect to ", baseEstimator, " (Type=",basePoolType,")") )
   }
   return(re)
}
#============================================
 
#============================================
# eVar under N(0,1)
#
# m >= 2 : only unbiased case is considered
# When m=1, it is compatible with evar ***
# poolType = Type of pooling
# ** For range, the calculation is not accurate when n > 100.
#--------------------------------------------
# source("/home/cp/MyFiles/research/R-package/rQCC-mw22-mw4/rQCC/R/factors.R")
evar <-
function(n, estimator=c("mean","median","HL1","HL2","HL3", "sd", "range", "mad","shamos"),
         poolType=c("A","B","C"), correction=TRUE)
{
   n = as.integer(n)
   N = sum(n)

   estimator = match.arg(estimator)
   poolType = match.arg(poolType)

   ## if ( any(n <= 0L) ) return(NA)
   if ( any(n <= 0L) )  {
      message("WARNING: sample size cannot be smaller than one.")
      return(NA)
   }

   ## if ( any(n <= 1L) && (estimator %in% c("HL1","sd", "range", "shamos")) )   return(NA)
   if ( any(n <= 1L) && (estimator %in% c("HL1","sd", "range", "shamos")) ) {
       message("WARNING: for HL1, sd, range and shamos, all the sample size should be at least two.")
      return(NA)
   }

   #-----------------------------------
   m = length(n)
   V = c4.factor.estimator = numeric(m)
   u = ifelse(m >= 2L, 1, as.numeric(correction))
   #-----------------------------------
   for( i in seq_len(m) ) {
      ni = n[i]
      c4sd = c4.factor(ni,"sd")
      one.m.c4sdsq = 1-c4sd^2
      d3 = factors.cc(ni,factor="d3")
      d3sq = d3^2
      if ( ni == 1L ) {
        V[i] = switch(estimator,
               mean   = 1/ni,
               median = (n.times.eVar.of.median[ni]/ni),
               HL1    = (n.times.eVar.of.HL1[ni]/ni),
               HL2    = (n.times.eVar.of.HL2[ni]/ni),
               HL3    = (n.times.eVar.of.HL3[ni]/ni),
               range  = 0,
               mad    = 0,
               shamos = 0
             )
      } else if ( ni <= 100L ) {
        V[i] = switch(estimator,
               mean   = 1/ni,
               median = (n.times.eVar.of.median[ni]/ni),
               HL1    = (n.times.eVar.of.HL1[ni]/ni),
               HL2    = (n.times.eVar.of.HL2[ni]/ni),
               HL3    = (n.times.eVar.of.HL3[ni]/ni),
               sd     = one.m.c4sdsq / (1-u+u*c4sd^2),
               range  = d3sq / (1-u+u*(c4.factor(ni,"range"))^2),
               mad    = (n.times.eVar.of.mad[ni]/ni)    /
                        (1-u+u*(c4.factor(ni,"mad"))^2),
               shamos = (n.times.eVar.of.shamos[ni]/ni) /
                        (1-u+u*(c4.factor(ni,"shamos"))^2)
             )
      } else {
        V[i] = switch(estimator,
               mean   = 1/ni,
               median ={tmp=ifelse(ni%%2L==1L,-0.6589-0.9430/ni,
                            -2.195+1.929/ni); (1.57+tmp/ni)/ni},
               HL1    = ((1.0472 + 0.1127/ni + 0.8365/ni^2)/ni),
               HL2    = ((1.0472 + 0.2923/ni + 0.2258/ni^2)/ni),
               HL3    = ((1.0472 + 0.2022/ni + 0.4343/ni^2)/ni),
               sd     = one.m.c4sdsq / (1-u+u*c4sd^2),
               range  = d3sq / (1-u+u*(c4.factor(ni,"range"))^2),
               mad    = {tmp=ifelse(ni%%2L==1L,0.2996-149.357/ni,
                             -2.417-153.01/ni);
                        (one.m.c4sdsq*(2.7027+tmp/ni)) /
                        (1-u+u*(c4.factor(ni,"mad"))^2) },
               shamos = (one.m.c4sdsq*(1.15875+2.822/ni+12.238/ni^2)) /
                        (1-u+u*(c4.factor(ni,"shamos"))^2)
             )
      }
   }
   if ( m == 1L ) return(V)
   if (estimator %in% c("mean","median","HL1","HL2","HL3") )  {
       # location estimators are all unbiased.
       Vpool = switch(poolType,
              A = ( sum(V)/(m^2) ),
              B = ( sum(n^2*V)/N^2 ),
              C = ( 1 / sum(1/V) )
      )
      return(Vpool)
   }
   if (estimator %in% c("sd","range","mad","shamos") )  {
       for( i in seq_len(m) ) c4.factor.estimator[i] = c4.factor(n[i],estimator)
       # NOTE: V (m>1) is unibased here. V_i/C_i^2(paper) -> V_i here.
       Vpool = switch(poolType,
             A = ( sum(V)/(m^2) ),
             B = ( sum(V*c4.factor.estimator^2) / sum(c4.factor.estimator)^2 ),
             C = ( 1 / sum(1/V) )
      )
      return(Vpool)
   }
}
#============================================

