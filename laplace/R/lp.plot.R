lp.plot <- function(x, plot.it=TRUE, a, col.line="black", lty.line=1,
                    xlim = NULL, ylim = NULL, main = NULL, sub = NULL, xlab = NULL, ylab="Probability", ...) {
   x = sort(x[complete.cases(x)])
   n = length(x)
   if (missing(a)) {
       a = ifelse(n <= 10, 3/8, 1/2)
   }
   p.minus.5 = ppoints(n) - 0.5
   d = sign(p.minus.5)
   y = -d*log(1-2*abs(p.minus.5))

   if (plot.it) {
      plot(x,y, yaxt="n", xlim=xlim, ylim=ylim, main=main, sub=sub, xlab=xlab, ylab=ylab, ...)
      ticklabels=c( (1:5)/100, (1:9)/10, 0.95, 0.99)
      ticklabels.minus.5 = ticklabels - 0.5
      dt = sign(ticklabels.minus.5)
      ticksat = -dt*log(1-2*abs( ticklabels.minus.5)) 
      axis(2,at=ticksat,labels=ticklabels)
      LM = lm(y~x)
      curve( predict(LM,newdata=data.frame(x)), add=TRUE, lty=lty.line, col=col.line)
   }
   invisible(list(x=x, y=y))
}

