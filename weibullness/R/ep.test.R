globalVariables("Exponential.ANOVA.Quantiles")
# Exponential GOF Test from the Exponential ANOVA
ep.test <- function (x,a) {
   DNAME = deparse(substitute(x))
   stopifnot(is.numeric(x)) 
   x = sort(x[complete.cases(x)])
   n =  length(x)
   if ((n < 3L || n > 1000L)) stop("sample size must be between 3 and 1000")
   if (missing(a))  {a= ifelse(n <= 10, 3/8, 1/2)}

   w = n /(n-1) * (mean(x)-min(x))^2 / ((n-1)*var(x))

   alphas = as.numeric( colnames(Exponential.ANOVA.Quantiles) )
   quantiles = Exponential.ANOVA.Quantiles[n-2,]

   pvalue.function = approxfun(x=quantiles, y=alphas, yleft=0.0, yright=1.0) 
   pval = (1-pvalue.function(w))  ## <- 

   RVAL = list(statistic = c(W=w), p.value = pval, sample.size=n, 
           method = "Exponential GOF test from ANOVA for Exponential", data.name = DNAME)
    class(RVAL) = "htest"
    return(RVAL)
}
## ep.test(1:3)

# p-value with w and n
ep.test.pvalue <- function (w,n) { 
   if (n < 3L || n > 1000L) stop("sample size must be between 3 and 1000")

   alphas = as.numeric( colnames(Exponential.ANOVA.Quantiles) )
   quantiles = Exponential.ANOVA.Quantiles[n-2,]

   pvalue.function = approxfun(x=quantiles, y=alphas, yleft=0.0, yright=1.0) 
   pval = (1-pvalue.function(w))  ## <- 

   RVAL = list(statistic = c(W = w), p.value = pval,        
           method = "p-value value for the Exponential GOF test",
           data.name = NULL)
    class(RVAL) = "htest"
    return(RVAL)
}

# critical value with alpha and n
ep.test.critical <- function (alpha, n) { 
   if (n < 3L || n > 1000L) stop("sample size must be between 3 and 1000")
   if (alpha <= 0 || alpha >= 1) stop("Significance level must be between 0 and 1")
   alpha = (1-alpha)      ###
   AL = as.integer( round(alpha*1000) ) 

   alphas = as.numeric( colnames(Exponential.ANOVA.Quantiles) )
   quantiles = Exponential.ANOVA.Quantiles[n-2,AL+1]

   RVAL = list(sample.size=n, alpha= (1-AL/1000), critical.value=quantiles, 
               title="Critical value for the Exponential GOF test", 
               data.name =Exponential.ANOVA.Quantiles)
   class(RVAL) = "ep.test.critical"
   return(RVAL)
}

## critical value for ep.test 
print.ep.test.critical <- function (x,...) { 
   cat("\n       ", x$title, "\n\n")
   cat("significance level = ",x$alpha, ", ",
       "sample size = ", x$sample.size,"\n\n",sep="" )
   cat("critical value =",   x$critical.value, "\n\n")
}

