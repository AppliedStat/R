


###############################################################################
wp.plot <- function(x, plot.it=TRUE, a, col.line="black", lty.line=1, 
                    xlim = NULL, ylim = NULL, main = NULL, sub = NULL, xlab = NULL, ylab="Probability", ...) {
   x = sort(x[complete.cases(x)])
   n = length(x)
   if (missing(a)) {
       a = ifelse(n <= 10, 3/8, 1/2)
   }
   y = log(-log(1-ppoints(n,a=a)))
   if (plot.it) {
      plot(x,y, log="x", yaxt="n", xlim=xlim, ylim=ylim, main=main, sub=sub, xlab=xlab, ylab=ylab, ...)
      ticklabels=c( (1:5)/100, (1:9)/10)
      ticksat=log(-log(1-ticklabels))
      axis(2,at=ticksat,labels=ticklabels)
      LM = lm(y~log(x)) 
      curve( predict(LM,newdata=data.frame(x)), add=TRUE, lty=lty.line, col=col.line)
   }
   invisible(list(x=x, y=y))
}
###############################################################################

