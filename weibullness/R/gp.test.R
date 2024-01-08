globalVariables("Gumbel.Plot.Quantiles")
# Gumbel Probability Plot
gp.test <- function (x,a) {
   DNAME = deparse(substitute(x))
   stopifnot(is.numeric(x)) 
   x = sort(x[complete.cases(x)])
   n =  length(x)
   if ((n < 3L || n > 1000L)) stop("sample size must be between 3 and 1000")
   if (missing(a))  {a= ifelse(n <= 10, 3/8, 1/2)}

   u = x
   v =-log(-log(ppoints(n,a=a)))
   r = cor(u,v) 

   alphas = as.numeric( colnames(Gumbel.Plot.Quantiles) )
   quantiles = Gumbel.Plot.Quantiles[n-2,]

   pvalue.function = approxfun(x=quantiles, y=alphas, yleft=0.0, yright=1.0) 
   pval = pvalue.function(r)

   RVAL = list(statistic = c(correlation=r), p.value = pval, sample.size=n, 
           method = "Gumbel goodness-of-fit test from the Gumbel probability plot", data.name = DNAME)
    class(RVAL) = "htest"
    return(RVAL)
}
## gp.test(1:3)

# p-value with r and n
gp.test.pvalue <- function (r,n) { 
   if (n < 3L || n > 1000L) stop("sample size must be between 3 and 1000")

   alphas = as.numeric( colnames(Gumbel.Plot.Quantiles) )
   quantiles = Gumbel.Plot.Quantiles[n-2,]

   pvalue.function = approxfun(x=quantiles, y=alphas, yleft=0.0, yright=1.0) 
   pval = pvalue.function(r)

   RVAL = list(statistic = c(correlation = r), p.value = pval,        
           method = "p-value value for the Gumbel goodness-of-fit test",
           data.name = NULL)
    class(RVAL) = "htest"
    return(RVAL)
}

## critical value for gp.test 
print.gp.test.critical <- function (x,...) { 
   cat("\n       ", x$title, "\n\n") 
   cat("significance level = ",x$alpha, ", ", 
       "sample size = ", x$sample.size,"\n\n",sep="" )
   cat("critical value =",   x$critical.value, "\n\n")
}

# critical value with alpha and n
gp.test.critical <- function (alpha, n) { 
   if (n < 3L || n > 1000L) stop("sample size must be between 3 and 1000")
   if (alpha <= 0 || alpha >= 1) stop("Significance level must be between 0 and 1")
   AL = as.integer( round(alpha*1000) ) 

   alphas = as.numeric( colnames(Gumbel.Plot.Quantiles) )
   quantiles = Gumbel.Plot.Quantiles[n-2,AL+1]

   RVAL = list(sample.size=n, alpha=AL/1000, critical.value=quantiles, 
               title="Critical value for the Gumbel goodness-of-fit test", 
               data.name =Gumbel.Plot.Quantiles)
   class(RVAL) = "gp.test.critical"
   return(RVAL)
}
