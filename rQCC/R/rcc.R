#============================================
# Robust Control Chart
#--------------------------------------------
rcc <- function(x, location=c("mean", "median","HL1","HL2","HL3"),
               scale=c("sd", "range", "mad","shamos"), 
               correction=TRUE,sigma.factor=3) 
{
 location = match.arg(location)
    scale = match.arg(scale)
        n = ncol(x); m = nrow(x)

xbari = switch (location, 
           mean  = apply(x,MARGIN=1, FUN=mean), 
           median= apply(x,MARGIN=1, FUN=median), 
           HL1   = apply(x,MARGIN=1, function(x)HL(x,method="HL1")), 
           HL2   = apply(x,MARGIN=1, function(x)HL(x,method="HL2")),
           HL3   = apply(x,MARGIN=1, function(x)HL(x,method="HL3")) )
xbarbar = mean(xbari)

si = switch (scale,
              sd = apply(x,MARGIN=1, FUN=sd), 
           range = apply(x,MARGIN=1, function(x){diff(range(x))}), 
             mad = apply(x,MARGIN=1, FUN=mad),
          shamos = apply(x,MARGIN=1, FUN=shamos) )
sbar = mean(si)

  if (correction==TRUE) {
     SE = switch(scale, 
               sd = A3.factor(n,"sd",    sigma.factor=sigma.factor)*sbar, 
            range = A3.factor(n,"range", sigma.factor=sigma.factor)*sbar, 
              mad = A3.factor(n,"mad",   sigma.factor=sigma.factor)*sbar, 
           shamos = A3.factor(n,"shamos",sigma.factor=sigma.factor)*sbar )
  } else {
     SE = sigma.factor/sqrt(n)*sbar
  }

  UCL = xbarbar + SE
   CL = xbarbar
  LCL = xbarbar - SE 
  res <- list(LCL=LCL,CL=CL,UCL=UCL, data=x, n=n, m=m, xbari=xbari, xbarbar=xbarbar, sbar=sbar,
              si=si, SE=SE, loc=location, scale=scale, bc=correction, sig=sigma.factor)
  class(res) <- "rcc"
  res
}
#------------------
print.rcc <- function(x, ...) { 
   limits = c(x$LCL, x$CL, x$UCL)
   names(limits) = c("LCL", "CL", "UCL")
   print(limits, quote=FALSE)
   invisible(x)
}
#------------------
summary.rcc <- function(object, ...) { 
   z=object
   digits = getOption("digits") 

   TITLE = paste0(" Control chart using location (", z$loc, ")", " and ", 
                  "scale (", z$scale, ")", " with ", z$sig, "*sigma." )
   cat("\n", TITLE, "\n\n")

 # Basic info. Sample info. 
   cat("Basic Information:", "\n")
   if (z$bc) { 
       OUT = c(paste(c(z$n, z$m, round(A3.factor(z$n, z$scale, z$sig),digits))),  z$bc )
   } else {
       OUT = c(paste(c(z$n, z$m, round(z$sig/sqrt(z$n),digits))), z$bc)
   }
   names(OUT) = c("Sample_size", "Subgroup_size", "A3 factor", "Bias_correction")
   print(OUT, quote=FALSE, print.gap=1)
   cat("\n")

 # location
   cat("Location:", "\n")
   OUT = c(round(z$xbarbar,digits),  
       paste0(round(100*finite.breakdown(z$n,z$loc),3)), 
       paste0(round(100*RE(z$n,z$loc,correction=z$bc),3)) ) 
   names(OUT) = c("Estimate", "  Breakdown(%)", "RE(%)")
   print(OUT, quote=FALSE)
   cat("\n")

 # scale
   cat("Scale:", "\n")
   OUT = c(round(z$sbar,digits),
       paste0(round(100*finite.breakdown(z$n,z$scale),3)),  
       paste0(round(100*RE(z$n,z$scale,correction=z$bc),3)) )  
   names(OUT) = c("Estimate", "  Breakdown(%)", "RE(%)")
   print( noquote(OUT) )
   cat("\n")

 # LCL, CL, UCL
   cat("Control Limits:", "\n")
   limits = c(z$LCL, z$CL, z$UCL)
   names(limits) = c("LCL", "CL", "UCL")
   print(limits, quote=FALSE)
   cat("\n")
}
# summary(result)

#------------------
plot.rcc <- function(x, digits= getOption("digits")-2, cex.text=0.7, x.text=1, ...) {
   LCL = x$LCL; CL=x$CL;  UCL = x$UCL; SE=x$SE;
   xbar = apply(x$data,MARGIN=1, FUN=mean)
   ymin = min(xbar,LCL); ymax = max(xbar,UCL)

   m = length(xbar)
   subgroup = 1:m
   index1 = (xbar > LCL) & (xbar < UCL)  ## in-control
   index2 = !index1                      ## out-of-control 

   plot(subgroup[index1], xbar[index1], xlim=c(1, m), ylim=c(ymin-0.2*SE, ymax+0.2*SE),
        type="p", pch=1, frame=FALSE, axes=FALSE, ylab=NA, xlab="Subgroup" )
   xx = c(0.3,   0.3, m+0.5, m+0.6, 0.6)
   yy = c(LCL, UCL, UCL,          LCL,          LCL)
   polygon(xx,yy, col="#FFFFEF", border=FALSE)
   points( subgroup[index1], xbar[index1], pch=20, col="blue" )
   points( subgroup[index2], xbar[index2], pch=4,  col="red" )

   if (missing(cex.text)) cex.text = 0.7
   text(x.text, UCL+0.1*cex.text*SE, paste0("UCL=", round(UCL,digits)), adj=0, cex=cex.text) 
   text(x.text,  CL+0.1*cex.text*SE, paste0(" CL=",  round(CL,digits)), adj=0, cex=cex.text) 
   text(x.text, LCL-0.1*cex.text*SE, paste0("LCL=", round(LCL,digits)), adj=0, cex=cex.text) 

   lines(subgroup, xbar, lty=3)
   abline( h= CL, lty=2, lwd=0.5 )
   abline( h=c(LCL,UCL), lwd=1.0 )
   axis(1, at =subgroup, subgroup, lty=1, lwd=0.5)
   axis(2, lty=1)
}
#============================================
