globalVariables("bs.plot.quantiles")
# Goodness of fit Test for Birnbaum-Saunders Distribution from the Probability Plot
bs.test <- function (x,a) {
   DNAME = deparse(substitute(x))
   stopifnot(is.numeric(x)) 
   x = sort(x[complete.cases(x)])
   n =  length(x)
   if ((n < 3L || n > 1000L)) stop("sample size must be between 3 and 1000")
   if (missing(a))  {a= ifelse(n <= 10, 3/8, 1/2)}

   # u = x
   # v = sqrt(x) * qnorm(ppoints(n,a=a))
   # r = cor(u,v) 
   r = cor( x, sqrt(x)*qnorm(ppoints(n,a=a)) )  

   alphas = as.numeric( colnames(bs.plot.quantiles) )
   quantiles = bs.plot.quantiles[n-2,]

   pvalue.function = approxfun(x=quantiles, y=alphas, yleft=0.0, yright=1.0) 
   pval = pvalue.function(r)

   RVAL = list(statistic = c(correlation=r), p.value = pval, sample.size=n, 
              method = "Goodness of fit test for Birnbaum-Saunders distribution", 
              data.name = DNAME)
    class(RVAL) = "htest"
    return(RVAL)
}
## bs.test(1:3)

bs.test.pvalue <- function (r,n) { 
   if (n < 3L || n > 1000L) stop("sample size must be between 3 and 1000")
   if (r <= 0 || r > 1)     stop("sample correlation must be between 0 and 1")

   alphas = as.numeric( colnames(bs.plot.quantiles) )
   quantiles = bs.plot.quantiles[n-2,]

   pvalue.function = approxfun(x=quantiles, y=alphas, yleft=0.0, yright=1.0) 
   pval = pvalue.function(r)

   RVAL = list(statistic = c(correlation = r), p.value = pval,        
              method = "p-value for Goodness of fit test for Birnbaum-Saunders distribution", 
              data.name = NULL)
    class(RVAL) = "htest"
    return(RVAL)
}

## critical value for bs.test 
print.bs.test.critical.value <- function (x,...) { 
   cat("\n       ", x$title, "\n\n") 
   cat("significance level = ",x$alpha, ", ", 
       "sample size = ", x$sample.size,"\n\n",sep="" )
   cat("critical value =",   x$critical.value, "\n\n")
}

bs.test.critical <- function (alpha, n) { 
   if (n < 3L || n > 1000L) stop("sample size must be between 3 and 1000")
   if (alpha <  0 || alpha >  1) stop("Significance level must be between 0 and 1")
   AL = as.integer( round(alpha*1000) ) 

   alphas = as.numeric( colnames(bs.plot.quantiles) )
   quantiles = bs.plot.quantiles[n-2,AL+1]

   RVAL = list(sample.size=n, alpha=AL/1000, critical.value=quantiles, 
               title="Critical value for the Goodness of fit test for Birnbaum-Saunders distribution")
   class(RVAL) = "bs.test.critical.value"
   return(RVAL)
}

#####################################################################
bs.plot = function(x, plot.it=TRUE, a, col.line="black", lty.line=1, 
    xlim = NULL, ylim = NULL, main = NULL, sub = NULL, xlab = NULL, 
    ylab = "Probability",  ...) {
   x = sort(x[complete.cases(x)])
   n = length(x)
   if (missing(a))  {a= ifelse(n <= 10, 3/8, 1/2)}
   p = ppoints(n,a=a)
   w = sqrt(x)*qnorm(p)
   if (plot.it) { 
      plot(x, w, type="p", yaxt = "n", xlim = xlim, ylim = ylim, 
            main = main, sub = sub, xlab = xlab, ylab = ylab, ...)
      ticklabels=c(0.01, seq(0.1,0.9,by=0.1), seq(0.91,0.99,by=0.01) )
      qq = quantile(x, probs=ticklabels)
      ticksat= qnorm(ticklabels)* sqrt( qq )
      axis(2,at=ticksat,labels=ticklabels)
      LM = lm(w~x) 
      curve(predict(LM, newdata = data.frame(x)), add = TRUE, 
            lty = lty.line, col = col.line)
   }
   invisible(list(x = x, w = w))
}

