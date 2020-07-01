#============================================
# Robust Control Chart
#--------------------------------------------
rcc <- function(x, location=c("mean", "median","HL1","HL2","HL3"),
       scale=c("sd", "range", "mad","shamos"), type=c("Xbar","S", "R"), sigma.factor=3) 
{
 location = match.arg(location)
    scale = match.arg(scale)
     type = match.arg(type)
        n = ncol(x); m = nrow(x)
   #--------------------------------------------

   if (type == "Xbar") {
   # Xbar chart
      xbari = switch (location, 
                 mean  = apply(x,MARGIN=1, FUN=mean), 
                 median= apply(x,MARGIN=1, FUN=median), 
                 HL1   = apply(x,MARGIN=1, function(x)HL(x,method="HL1")), 
                 HL2   = apply(x,MARGIN=1, function(x)HL(x,method="HL2")),
                 HL3   = apply(x,MARGIN=1, function(x)HL(x,method="HL3")) )
      xbarbar = mean(xbari)
      
        # unbiased scale estimate
        si.unbiased = switch (scale,
                    sd = apply(x,MARGIN=1, FUN=sd)/c4.factor(n,"sd"), 
                 range = apply(x,MARGIN=1, function(x){diff(range(x))})/c4.factor(n,"range"), 
                   mad = apply(x,MARGIN=1, FUN=mad.unbiased),
                shamos = apply(x,MARGIN=1, FUN=shamos.unbiased) )
       sigma.hat = mean(si.unbiased)
      
         SE = sigma.factor/sqrt(n)*sigma.hat
      
        UCL = xbarbar + SE
         CL = xbarbar
        LCL = xbarbar - SE 
        res <- list(LCL=LCL,CL=CL,UCL=UCL, data=x, n=n, m=m, 
                    xbari=xbari, xbarbar=xbarbar, sigma.hat=sigma.hat,
                    SE=SE, loc=location, scale=scale, type=type, sig=sigma.factor)
        class(res) <- "rcc"
        return(res)
   } else if (type == "S") {
   # S chart
        # unbiased scale estimate
        si.unbiased = switch (scale,
                    sd = apply(x,MARGIN=1, FUN=sd)/c4.factor(n,method="sd"), 
                 range = apply(x,MARGIN=1, function(x){diff(range(x))})/c4.factor(n,method="range"),
                   mad = apply(x,MARGIN=1, FUN=mad.unbiased),
                shamos = apply(x,MARGIN=1, FUN=shamos.unbiased) )
       sigma.hat = mean(si.unbiased)
   
        c4 = c4.factor(n,method="sd")
        se = sigma.factor*sqrt(1-c4^2)
        B1 = max(c4-se,0); B2 = c4+se
         
        UCL = B2 * sigma.hat
         CL = c4 * sigma.hat
        LCL = B1 * sigma.hat 
        res <- list(LCL=LCL,CL=CL,UCL=UCL, data=x, n=n, m=m, sigma.hat=sigma.hat,
                    SE=se*sigma.hat, scale=scale, type=type, sig=sigma.factor)
        class(res) <- "rcc"
        return(res)
   } else {
   # R chart
        d2 = c4.factor(n, method="range")
        d3 = factors.cc(n,factor="d3")

        # unbiased scale estimate
        Ri.unbiased = switch (scale,
                    sd = apply(x,MARGIN=1, FUN=sd)/c4.factor(n,"sd"),
                 range = apply(x,MARGIN=1, function(x){diff(range(x))})/d2,
                   mad = apply(x,MARGIN=1, FUN=mad.unbiased),
                shamos = apply(x,MARGIN=1, FUN=shamos.unbiased) )
        sigma.hat = mean(Ri.unbiased)
        se = sigma.factor*d3
        D1 = max(d2-se,0); D2 = d2+se

        UCL = D2 * sigma.hat
         CL = d2 * sigma.hat
        LCL = D1 * sigma.hat
        res <- list(LCL=LCL,CL=CL,UCL=UCL, data=x, n=n, m=m, sigma.hat=sigma.hat,
                    SE=se*sigma.hat, scale=scale, type=type, sig=sigma.factor)
        class(res) <- "rcc"
        return(res)
   }
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

   if (z$type == "Xbar") {
   # Xbar chart
      TITLE = paste0(" Xbar Control chart using location (", z$loc, ")", " and ", 
                     "scale (unbiased ", z$scale, ")", " with ", z$sig, "*sigma." )
      cat("\n", TITLE, "\n\n")
   
      # Basic info. Sample info. 
      cat("Basic Information:", "\n")
      OUT = c(paste(c(z$n, z$m)))
      names(OUT) = c("Sample_size", "Subgroup_size" )
      print(OUT, quote=FALSE, print.gap=1)
      cat("\n")
   
      # location
      cat("Location:", "\n")
      OUT = c(round(z$xbarbar,digits),  
          paste0(round(100*finite.breakdown(z$n,z$loc),3)), 
          paste0(round(100*RE(z$n,z$loc,correction=TRUE),3)) ) 
      names(OUT) = c("Estimate", "Breakdown(%)", "RE(%)")
      print(OUT, quote=FALSE)
      cat("\n")
   
      # scale
      cat("Scale:", "\n")
      OUT = c(round(z$sigma.hat,digits),
          paste0(round(100*finite.breakdown(z$n,z$scale),3)),  
          paste0(round(100*RE(z$n,z$scale,correction=TRUE),3)) )  
      names(OUT) = c("Estimate", "Breakdown(%)", "RE(%)")
      print( noquote(OUT) )
      cat("\n")
   
      # LCL, CL, UCL
      cat("Control Limits:", "\n")
      limits = c(z$LCL, z$CL, z$UCL)
      names(limits) = c("LCL", "CL", "UCL")
      print(limits, quote=FALSE)
      cat("\n")
   } else if (z$type == "S") {
   # S chart
      TITLE = paste0(" S Control chart using scale (unbiased ", z$scale, ")", 
                     " with ", z$sig, "*sigma." )
      cat("\n", TITLE, "\n\n")
    
      # Basic info. Sample info.  
      cat("Basic Information:", "\n")
      OUT = c(paste(c(z$n, z$m)))
      names(OUT) = c("Sample_size", "Subgroup_size" )
      print(OUT, quote=FALSE, print.gap=1)
      cat("\n")
    
      # scale
      cat("Scale:", "\n")
      OUT = c(round(z$sigma.hat,digits),
          paste0(round(100*finite.breakdown(z$n,z$scale),3)),  
          paste0(round(100*RE(z$n,z$scale,correction=TRUE),3)) )
      names(OUT) = c("Estimate", "Breakdown(%)", "RE(%)")
      print( noquote(OUT) )
      cat("\n")
    
      # LCL, CL, UCL
      cat("Control Limits:", "\n")
      limits = c(z$LCL, z$CL, z$UCL) 
      names(limits) = c("LCL", "CL", "UCL")
      print(limits, quote=FALSE)
      cat("\n")
   } else {
   # R chart
      TITLE = paste0(" R Control chart using scale (unbiased ", z$scale, ")",
                     " with ", z$sig, "*sigma." )
      cat("\n", TITLE, "\n\n")
      
      # Basic info. Sample info.  
      cat("Basic Information:", "\n")
      OUT = c(paste(c(z$n, z$m)))
      names(OUT) = c("Sample_size", "Subgroup_size" )
      print(OUT, quote=FALSE, print.gap=1)
      cat("\n")
      
      # scale
      cat("Scale:", "\n")
      OUT = c(round(z$sigma.hat,digits),
          paste0(round(100*finite.breakdown(z$n,z$scale),3)),  
          paste0(round(100*RE(z$n,z$scale,correction=TRUE),3)) )
      names(OUT) = c("Estimate", "Breakdown(%)", "RE(%)")
      print( noquote(OUT) )
      cat("\n")
      
      # LCL, CL, UCL
      cat("Control Limits:", "\n")
      limits = c(z$LCL, z$CL, z$UCL) 
      names(limits) = c("LCL", "CL", "UCL")
      print(limits, quote=FALSE)
      cat("\n")
   }
}
   
#------------------
plot.rcc <- function(x, digits= getOption("digits")-2, cex.text=0.7, x.text=1, ...) {
   LCL = x$LCL; CL=x$CL;  UCL = x$UCL; SE=x$SE;

   if (x$type == "Xbar") {
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
   
      text(x.text, UCL+0.1*cex.text*SE, paste0("UCL=", round(UCL,digits)), adj=0, cex=cex.text) 
      text(x.text,  CL+0.1*cex.text*SE, paste0(" CL=",  round(CL,digits)), adj=0, cex=cex.text) 
      text(x.text, LCL-0.1*cex.text*SE, paste0("LCL=", round(LCL,digits)), adj=0, cex=cex.text) 
   
      lines(subgroup, xbar, lty=3)
      abline( h= CL, lty=2, lwd=0.5 )
      abline( h=c(LCL,UCL), lwd=1.0 )
      axis(1, at =subgroup, subgroup, lty=1, lwd=0.5)
      axis(2, lty=1)
   } else if (x$type == "S") {
      si = apply(x$data,MARGIN=1, FUN=sd)
      ymin = min(si,LCL); ymax = max(si,UCL)
   
      m = length(si)
      subgroup = 1:m
      index1 = (si > LCL) & (si < UCL)  ## in-control
      index2 = !index1                  ## out-of-control 
   
      plot(subgroup[index1], si[index1], xlim=c(1, m), ylim=c(ymin-0.2*SE, ymax+0.2*SE),
           type="p", pch=1, frame=FALSE, axes=FALSE, ylab=NA, xlab="Subgroup" )
      xx = c(0.3,   0.3, m+0.5, m+0.6, 0.6)
      yy = c(LCL, UCL, UCL,          LCL,          LCL)
      polygon(xx,yy, col="#FFFFEF", border=FALSE)
      points( subgroup[index1], si[index1], pch=20, col="blue" )
      points( subgroup[index2], si[index2], pch=4,  col="red" )
   
      text(x.text, UCL+0.1*cex.text*SE, paste0("UCL=", round(UCL,digits)), adj=0, cex=cex.text)
      text(x.text,  CL+0.1*cex.text*SE, paste0(" CL=",  round(CL,digits)), adj=0, cex=cex.text)
      text(x.text, LCL-0.1*cex.text*SE, paste0("LCL=", round(LCL,digits)), adj=0, cex=cex.text)
   
      lines(subgroup, si, lty=3)
      abline( h= CL, lty=2, lwd=0.5 )
      abline( h=c(LCL,UCL), lwd=1.0 )
      axis(1, at =subgroup, subgroup, lty=1, lwd=0.5)
      axis(2, lty=1)
   } else {
      Ri = apply(x$data,MARGIN=1, FUN=function(x){diff(range(x))}) 
      ymin = min(Ri,LCL); ymax = max(Ri,UCL)
   
      m = length(Ri)
      subgroup = 1:m
      index1 = (Ri > LCL) & (Ri < UCL)  ## in-control
      index2 = !index1                  ## out-of-control 
   
      plot(subgroup[index1], Ri[index1], xlim=c(1, m), ylim=c(ymin-0.2*SE, ymax+0.2*SE),
           type="p", pch=1, frame=FALSE, axes=FALSE, ylab=NA, xlab="Subgroup" )
      xx = c(0.3,   0.3, m+0.5, m+0.6, 0.6)
      yy = c(LCL, UCL, UCL,          LCL,          LCL)
      polygon(xx,yy, col="#FFFFEF", border=FALSE)
      points( subgroup[index1], Ri[index1], pch=20, col="blue" )
      points( subgroup[index2], Ri[index2], pch=4,  col="red" )
   
      text(x.text, UCL+0.1*cex.text*SE, paste0("UCL=", round(UCL,digits)), adj=0, cex=cex.text)
      text(x.text,  CL+0.1*cex.text*SE, paste0(" CL=",  round(CL,digits)), adj=0, cex=cex.text)
      text(x.text, LCL-0.1*cex.text*SE, paste0("LCL=", round(LCL,digits)), adj=0, cex=cex.text)
   
      lines(subgroup, Ri, lty=3)
      abline( h= CL, lty=2, lwd=0.5 )
      abline( h=c(LCL,UCL), lwd=1.0 )
      axis(1, at =subgroup, subgroup, lty=1, lwd=0.5)
      axis(2, lty=1)
   }
}
#============================================
