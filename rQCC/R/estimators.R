#============================================
# Finite Breakdown 
#============================================
finite.breakdown <-
function(n, estimator=c("mean","median","HL1","HL2","HL3","sd","range","mad","shamos")) 
{ 
  n = as.integer(n)
  k = length(n)
  estimator = match.arg(estimator)
  breakdown = numeric(k)

  for ( i in 1L:k ) { 
      if ( estimator %in% c("mean","sd", "range") ) {
         breakdown[i] = 0
      } else if ( n[i] <= 1L ) {
         breakdown[i] = ifelse( n[i] <= 0L, NA, 0) 
      } else {
          switch (estimator,
           HL1 = { tmp = floor( (n[i]*(n[i]-1)-2)/4 )
                     m = floor( ((2*n[i]-1) - sqrt((2*n[i]-1)^2-8*tmp))/2 )},
           HL2 = { tmp = floor( (n[i]*(n[i]+1)-2)/4 )
                     m = floor( ((2*n[i]+1) - sqrt((2*n[i]+1)^2-8*tmp))/2 )},
           HL3 =  {tmp = floor(0.5*(n[i]^2-1))
                     m = floor( n[i] - sqrt(n[i]^2-tmp)) },
           median={  m = floor(0.5*(n[i]-1)) },
           mad   ={  m = floor(0.5*(n[i]-1)) },
           shamos={tmp = floor( (n[i]*(n[i]-1)-2)/4 )
                     m = floor( ((2*n[i]-1)-sqrt((2*n[i]-1)^2-8*tmp))/2) }
          )
          breakdown[i] = m/n[i] 
      }
  }
#message(paste0("Finite-sample breakdown point of ",estimator, " estimator"))
  names(breakdown) = paste0("n=", n)
  return(breakdown)
}
#============================================

#============================================
# Hodges-Lehmann estimators
# HL1: median over i < j
# HL2: median over i <= j
# HL3: median over all (i,j)
#============================================
HL <- function(x, estimator=c("HL1","HL2","HL3"),na.rm=FALSE) 
{
   estimator = match.arg(estimator)
   if (na.rm) x <- x[!is.na(x)]
   xx = outer(x,x, "+")

   HL.estimation =  switch (estimator,
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
shamos <- function(x,constant=1.048358, na.rm=FALSE, IncludeEqual=FALSE) 
{
  if (na.rm) x <- x[!is.na(x)]
  w1 = outer(x,x, "-")
  w2 = abs( w1[lower.tri(w1,diag=IncludeEqual)]  )
  constant * median( w2 )
}
#============================================

#============================================
# Unbiased scale estimators
#============================================
# Unbiased sd
sd.unbiased = function(x, na.rm = FALSE) 
{
    sqrt(var(if (is.vector(x)||is.factor(x)) x else as.double(x), 
         na.rm = na.rm)) / c4.factor(sum(!is.na(x)), estimator="sd")
}
#-----
# Unbiased mad
mad.unbiased = function(x,center=median(x),constant=1.4826,na.rm=FALSE)
{
  if (na.rm) x = x[!is.na(x)]
  n = length(x)
  if (n == 1L) return(0)
  constant*median(abs(x-center))/c4.factor(n,"mad") 
}
#-----
# Unbiased Shamos 
shamos.unbiased = function(x,constant=1.048358,na.rm=FALSE,IncludeEqual=FALSE)
{
  if (na.rm) x = x[!is.na(x)]
  n = length(x)
  if (n == 1L) return(0)
  w1 = outer(x,x, "-")
  w2 = abs( w1[lower.tri(w1,diag=IncludeEqual)]  )
  constant * median( w2 ) / c4.factor(n,"shamos") 
}
#-----
# Unbiased mad^2
mad2.unbiased = function(x,center=median(x),constant=1.4826,na.rm=FALSE)
{
  if (na.rm) x = x[!is.na(x)]
  n = length(x)
  if (n == 1L) return(0)
  mad(x, center, constant, na.rm)^2 / w4.factor(n,"mad2") 
}
#-----
# Unbiased Shamos^2
shamos2.unbiased = function(x,constant=1.048358,na.rm=FALSE,IncludeEqual=FALSE)
{
  if (na.rm) x = x[!is.na(x)]
  n = length(x)
  if (n == 1L) return(0)

  shamos(x, constant, na.rm, IncludeEqual)^2 / w4.factor(n,"shamos2")
}
#============================================



#============================================
# Calculate pooled estimators (pooling unbiased estimators)
#============================================
pooledEstimator <-
function(x, estimator=c("mean","median","HL1","HL2","HL3","sd","range","mad","shamos"), 
            poolType=c("A","B","C") )
{
   estimator = match.arg(estimator)
   poolType = match.arg(poolType)
   
   ## base:::mean.default # Don't use below b/c numeric in list is regarded as non-numeric.
   ## if(!is.numeric(x)) {
   ##     warning("argument is not numeric: returning NA")
   ##     return(NA_real_)
   ## }
   if ( !is.list(x) )  {  # Don't use ifelse
      if (is.matrix(x)) x = unname(split(x,row(x)))
      if (is.atomic(x)) x = unname(split(x,1))
   } 
         
   m = length(x)
   n = sapply(x,length)
   N = sum(n)
   W = numeric(m)

   # Estimates 
   EST = switch(estimator,
          mean   = sapply(x, mean),
          median = sapply(x, median),
          HL1    = mapply(HL, x, "HL1"),
          HL2    = mapply(HL, x, "HL2"),
          HL3    = mapply(HL, x, "HL3"),
          sd     = sapply(x, sd.unbiased),
          range  ={t1=apply(sapply(x,range),2,diff);
                   for(i in 1L:m)W[i]=c4.factor(n[i],"range");
                   t1/W},
          mad    = sapply(x, mad.unbiased),
          shamos = sapply(x, shamos.unbiased)
   )

   if (estimator %in% c("mean","median","HL1","HL2","HL3") ) {
      # Location
      WEIGHT <-
        switch(poolType,
               A = rep(1/m,m),
               B = n / N ,
               C = {for(i in 1L:m) W[i]=1/evar(n[i],estimator); W} )
   } else {
      # Scale 
      WEIGHT <-
        switch(poolType,
               A = rep(1/m,m),
               B ={for(i in seq_len(m)) W[i]=c4.factor(n[i],estimator);W},
               C ={for(i in seq_len(m)) W[i]=evar(n[i],estimator,poolType="C",correction=TRUE); 1/W} 
        )
   }

   return ( weighted.mean(EST,WEIGHT) )
}
#============================================

