
globalVariables( c("Quantiles.TA","Quantiles.TB") )

#=============================================================
# Robustified t-test using the TA and TB,
#     where TA (median and MAD) and TB (HL and Shamos).
#=============================================================
# 
#=============================================================
#- Hodges-Lehmann location estimator (i<j)
HL.estimate <- function(x, na.rm=FALSE, IncludeEqual=FALSE) {
  if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
      warning("argument is not numeric or logical: returning NA")
      return(NA_real_)
  }
  if (na.rm) x <- x[!is.na(x)]
  w = outer(x,x, "+")
  0.5 * median(w[lower.tri(w,diag=IncludeEqual)]) 
}
#=============================================================
# Only one sample. 
# Two-sample case (future work). 
rt.test <- 
function (x, alternative = c("two.sided", "less", "greater"), 
          mu=0, test.stat=c("TA", "TB"), conf.level=0.95 ) 
{
   DNAME = deparse(substitute(x))
   stopifnot(is.numeric(x)) 
   x = sort(x[complete.cases(x)])
   n =  length(x)

   if ((n < 4L || n > 100L)) stop("sample size must be between 4 and 100")
   if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
        conf.level<=0 || conf.level > 0.998 )) 
        stop("'conf.level' must be a single number between 0 and 0.998")
   alternative <- match.arg(alternative)
   test.stat   <- match.arg(test.stat) 
   
   if (test.stat == "TA") {   
       alphas = as.numeric( colnames(Quantiles.TA) )
       quantiles = Quantiles.TA[n-3,]

       med = median(x); const=0.7978845608028654057264 
       std.err =  mad(x, center=med) / (sqrt(n)*const)
       tstat = (med-mu) / std.err 
       estimate = median(x)
       names(estimate) = " median of x" 
       names(tstat) = "TA"
       method = "Robustified t-test using the median and MAD (TA)"
   } else if (test.stat == "TB") { 
       alphas = as.numeric( colnames(Quantiles.TB) )
       quantiles = Quantiles.TB[n-3,]

       const=0.9321290502845152   # const = sqrt( 6/pi) *qnorm(3/4) 
       w1=outer(x,x,"-"); w2=abs(w1[lower.tri(w1,diag=FALSE)]) # Shamos part
       std.err = median(w2) / (sqrt(n)*const)
       tstat = (HL.estimate(x)-mu) / std.err 

       estimate = HL.estimate(x)
       names(estimate) = " Hodges-Lehmann estimate of x" 
       names(tstat) = "TB"
       method = "Robustified t-test using the Hodges-Lehmann and Shamos (TB)"
   } else {
       stop("Choose the test statistic: \"TA\" or \"TB\".")
   }
 
   cdf.approx = approxfun(x=quantiles, y=alphas, yleft=0.0, yright=1.0) 
   cdf.function = function(q) ifelse( q>=0, cdf.approx(q), 1-cdf.approx(abs(q)) )
   quantile.function = function(conf.level) approxfun(x=conf.level, y=quantiles)  ###hhh

   if (alternative == "less") {
      pval = cdf.function(tstat) 
      cint = c(-Inf, tstat + q.robustified.t(conf.level,n,test.stat) )  
   }
   else if (alternative == "greater") {
      pval = 1 - cdf.function(tstat) 
      cint = c(tstat - q.robustified.t(conf.level,n,test.stat), Inf )   
   } 
   else {
      pval = 2*(1-cdf.function(abs(tstat)))
      cint = q.robustified.t((1+conf.level)/2, n,test.stat)
      cint = tstat + c(-cint,cint) 
   }
   cint = mu + cint*std.err 

   names(n)  = "sample size"
   names(mu) = "mean" 
   attr(cint, "conf.level") = conf.level 

   RVAL = list(statistic=tstat, p.value = pval, parameters=n, 
               conf.int=cint, sample.size=n,  estimate=estimate, null.value=mu,
               alternative=alternative, method=method, data.name=DNAME)
   class(RVAL) = "htest"
   return(RVAL)
}

## quantiles 
q.robustified.t = function (p, n, test.stat=c("TA","TB"), lower.tail=TRUE) { 
   if ((n < 4L || n > 100L)) stop("sample size must be between 4 and 100")
   test.stat  <- match.arg(test.stat) 

   if (test.stat == "TA") {
       alphas = as.numeric( colnames(Quantiles.TA) )
       quantiles = Quantiles.TA[n-3,]
   } else if (test.stat == "TB") {
       alphas = as.numeric( colnames(Quantiles.TB) )
       quantiles = Quantiles.TB[n-3,]
   } else {
       stop("Choose the test statistic: \"TA\" or \"TB\".")
   }

   quantile.approx = approxfun(x=alphas, y=quantiles) 
   index1 = (p >= 0.5) & (p <= 0.999)
   index2 = (p < 0.5)  & (p >= 0.001)
   index3 = p < 0.001 
   index4 = p > 0.999 
   q1 =  quantile.approx( p[index1] )
   q2 = -quantile.approx(1- p[index2] )
   q = numeric( length(p) )
   q[index1] = q1 
   q[index2] = q2 
   q[index3] = -Inf
   q[index4] =  Inf
   if(lower.tail) return(q) else return(-q) 
}

