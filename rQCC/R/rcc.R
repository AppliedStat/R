#============================================
# Robust Control Chart
#--------------------------------------------
rcc <- function(x, location=c("mean", "median","HL1","HL2","HL3"),
       scale=c("sd", "range", "mad","shamos"), type=c("Xbar","S", "R"), 
       poolLoc=c("A", "B", "C"), poolScale=c("A", "B", "C"), sigmaFactor=3, nk) 
{
  if (is.matrix(x)) x <- unname( split(x, row(x)) )
 location = match.arg(location)
    scale = match.arg(scale)
     type = match.arg(type)
  poolLoc = match.arg(poolLoc)
poolScale = match.arg(poolScale)

        m = length(x)
   nsubgr = sapply(x, length)
   is.balanced = (length(unique(nsubgr))==1L)

  LCLsubgr <- CLsubgr <- UCLsubgr <- numeric(m) 

  if (missing(nk)) nk = length(unlist(x))/m

  #--------------------------------------------
  if (type == "Xbar") {
  # Xbar chart
        xbarbar  = pooledEstimator(x, estimator=location, poolType=poolLoc) 
       sigma.hat = pooledEstimator(x, estimator=scale,    poolType=poolScale)

        SE = 1/sqrt(nk)*sigma.hat
       UCL = xbarbar + sigmaFactor*SE
        CL = xbarbar
       LCL = xbarbar - sigmaFactor*SE 

       SEsubgr  = 1/sqrt(nsubgr)*sigma.hat
       UCLsubgr = xbarbar + sigmaFactor*SEsubgr 
       LCLsubgr = xbarbar - sigmaFactor*SEsubgr

       res <- list(LCL=LCL, CL=CL, UCL=UCL, 
                   LCLsubgr=LCLsubgr, UCLsubgr=UCLsubgr,
                   data=x, nk=nk, m=m, nsubgr=nsubgr, is.balanced=is.balanced, 
                   xbarbar=xbarbar, sigma.hat=sigma.hat, SE=SE, loc=location, scale=scale, 
                   poolLoc=poolLoc, poolScale=poolScale, type=type, sig=sigmaFactor)
       class(res) <- "rcc"
       return(res)
  } else if (type == "S") {
  # S chart
       # unbiased scale estimator
       sigma.hat = pooledEstimator(x, estimator=scale,    poolType=poolScale)

       c4 = factors.cc(nk,factor="c4")   # c4.factor(nk,estimator="sd")
       se = sqrt(1-c4^2)
       B1 = max(c4-sigmaFactor*se,0); B2 = c4 + sigmaFactor*se

       UCL = B2 * sigma.hat
        CL = c4 * sigma.hat
       LCL = B1 * sigma.hat

       for ( i in seq_len(m) ) {
           c4 = factors.cc(nsubgr[i],factor="c4")   # c4.factor(nk,estimator="sd")
           se = sqrt(1-c4^2)
           B1 = max(c4-sigmaFactor*se,0); B2 = c4 + sigmaFactor*se
        
           UCLsubgr[i] = B2 * sigma.hat
            CLsubgr[i] = c4 * sigma.hat
           LCLsubgr[i] = B1 * sigma.hat
       }
       res <- list(LCL=LCL, CL=CL, UCL=UCL, 
                   LCLsubgr=LCLsubgr, UCLsubgr=UCLsubgr,
                   data=x, nk=nk, m=m, nsubgr=nsubgr, is.balanced=is.balanced,
                   sigma.hat=sigma.hat, SE=se*sigma.hat, scale=scale,
                   poolLoc=poolLoc, poolScale=poolScale, type=type, sig=sigmaFactor)
       class(res) <- "rcc"
       return(res)
  } else {
  # R chart
       # unbiased scale estimator
       sigma.hat = pooledEstimator(x, estimator=scale,    poolType=poolScale)

       d2 = factors.cc(nk,factor="d2")  # c4.factor(n,estimator="range")
       d3 = factors.cc(nk,factor="d3")

       se = d3
       D1 = max(d2-sigmaFactor*se,0); D2 = d2 + sigmaFactor*se

       UCL = D2 * sigma.hat
        CL = d2 * sigma.hat
       LCL = D1 * sigma.hat

       for ( i in seq_len(m) ) {
           d2 = factors.cc(nsubgr[i],factor="d2")  # c4.factor(n,estimator="range")
           d3 = factors.cc(nsubgr[i],factor="d3")

           se = d3
           D1 = max(d2-sigmaFactor*se,0); D2 = d2 + sigmaFactor*se
    
           UCLsubgr[i] = D2 * sigma.hat
            CLsubgr[i] = d2 * sigma.hat
           LCLsubgr[i] = D1 * sigma.hat

       }
       res <- list(LCL=LCL, CL=CL, UCL=UCL, 
                   LCLsubgr=LCLsubgr, CLsubgr=CLsubgr, UCLsubgr=UCLsubgr,
                   data=x, nk=nk, m=m, nsubgr=nsubgr, is.balanced=is.balanced, 
                   sigma.hat=sigma.hat, SE=se*sigma.hat, scale=scale,
                   poolLoc=poolLoc, poolScale=poolScale, type=type, sig=sigmaFactor)
       class(res) <- "rcc"
       return(res)
  }
}
#============================================
print.rcc <- function(x, ...) { 
   limits = c(x$LCL, x$CL, x$UCL)
   names(limits) = c("LCL", "CL", "UCL")
   print(limits, quote=FALSE)
   invisible(x)
}
#============================================
summary.rcc <- function(object, ...) {
   z=object
   digits = getOption("digits") 
   subgr.index = seq_len(z$m)

   if (z$type == "Xbar") {
   # Xbar 
      TITLE = paste0(" Xbar Control chart using location (", z$loc, ")", " and ",
                     "scale (unbiased ", z$scale, ")", " with ", z$sig, "*sigma." )
   ## cat("\n", "\033[1m", TITLE, "\033(B\033[m", "\n\n")
      cat("\n", TITLE, "\n\n")

      # Basic info. Sample info. 
      cat("=================", "\n")
      cat("Basic Information ", "\n")
      cat("-----------------", "\n")
      if (z$is.balanced) {
         cat("* Number of subgroups =", z$m, "\n")
         cat("* Sample size =",         z$nk,"\n")
      } else {
         cat("* Number of subgroups =", z$m, "\n")
         cat("* Subgroup sample sizes =", z$nsubgr, "\n")
      }
      cat("* Standard Error (SE) =", z$SE, "\n")
      cat("\n")

      # location  
      cat("=================", "\n")
      cat( paste0("Location (", z$loc, ")"), "\n")
      cat("-----------------", "\n")
      estimate = round(z$xbarbar,digits)
      if (z$is.balanced) {
         REpooled = suppressMessages(RE(n=z$nk, estimator=z$loc, correction=TRUE))
         OUT = c( estimate, paste0(round(100*finite.breakdown(z$nk,z$loc),2)), REpooled)
         names(OUT) = c("Estimate", "Breakdown(%)", "RE(%)")
         print(OUT, quote=FALSE)
      } else {
         REpooled = suppressMessages(RE(n=z$nsubgr, estimator=z$loc, poolType=z$poolLoc, correction=TRUE))
         subsample.sizes = sort(unique(z$nsubgr)) 
         REvalues <- BDvalues <- numeric(length(subsample.sizes))
         for ( i in seq_along(subsample.sizes) ) {
              BDvalues[i] = finite.breakdown(subsample.sizes[i],z$loc) 
              REvalues[i] = suppressMessages(RE(n=subsample.sizes[i],estimator=z$loc,correction=TRUE))
         }

         cat( paste0("* Pooled estimate of ", z$loc,  " (Type=",z$poolLoc,"): ", estimate), "\n")
         cat( paste0("* Minimum of subgroup breakdown points: ", round(100*min(BDvalues),2),"%"), "\n")
         cat( paste0("* RE of pooled ", z$loc, 
                     " (Type=",z$poolLoc,")"," with respect to mean", " (Type=",z$poolLoc,"): ", 
                       round(100*REpooled,2), "%"), "\n")
         cat("* Breakdown and RE for each subgroup sample size\n")
         OUT = rbind( round(100*BDvalues,2),  round(100*REvalues,2) ) 
         dimnames(OUT) = list( c("Breakdown (%)", "RE (%)"), paste0( "n=",subsample.sizes) )
         print(OUT, quote=FALSE)
      }
      cat("\n")
      
      # scale
      cat("=================", "\n")
      cat( paste0("Scale (", z$scale, ")"), "\n")
      cat("-----------------", "\n")
      estimate = round(z$sigma.hat,digits) 
      if (z$is.balanced) {
         REpooled = round(100*suppressMessages(RE(n=z$nk, estimator=z$scale, correction=TRUE)),2)
         OUT = c(estimate, paste0(round(100*finite.breakdown(z$nk,z$scale),2)), REpooled)
         names(OUT) = c("Estimate", "Breakdown(%)", "RE(%)")
         print(OUT, quote=FALSE)
      } else {
         REpooled = suppressMessages(RE(n=z$nsubgr, estimator=z$scale, poolType=z$poolScale, correction=TRUE))
         subsample.sizes = sort(unique(z$nsubgr))
         REvalues <- BDvalues <- numeric(length(subsample.sizes))
         for ( i in seq_along(subsample.sizes) ) {
              BDvalues[i] = finite.breakdown(subsample.sizes[i],z$scale)
              REvalues[i] = suppressMessages(RE(n=subsample.sizes[i],estimator=z$scale,correction=TRUE))
         }

         cat( paste0("* Pooled estimate of ", z$scale,  " (Type=",z$poolScale,"): ", estimate), "\n")
         cat( paste0("* Minimum of subgroup breakdown points: ", round(100*min(BDvalues),2),"%"), "\n")
         cat( paste0("* RE of pooled ", z$scale,
                     " (Type=",z$poolScale,")"," with respect to mean", " (Type=",z$poolScale,"): ",
                       round(100*REpooled,2), "%"), "\n")
         cat("* Breakdown and RE for each subgroup sample size\n")
         OUT = rbind( round(100*BDvalues,2),  round(100*REvalues,2) )
         dimnames(OUT) = list( c("Breakdown (%)", "RE (%)"), paste0( "n=",subsample.sizes) )
         print(OUT, quote=FALSE)
      }
      cat("\n")

      # LCL, CL, UCL
      cat("=================", "\n")
      cat( paste0("Control Limits (with sample size = ", z$nk, "):"), "\n" )
      cat("-----------------", "\n")
      limits = c(z$LCL, z$CL, z$UCL)
      names(limits) = c("LCL", "CL", "UCL")
      print(limits, quote=FALSE)
      cat("\n")
   } else if (z$type == "S") {
   # S chart
      TITLE = paste0(" S Control chart using scale (unbiased ", z$scale, ")",
                     " with ", z$sig, "*sigma." )
   ## cat("\n", "\033[1m", TITLE, "\033(B\033[m", "\n\n")
      cat("\n", TITLE, "\n\n")

      # Basic info. Sample info.  
      cat("=================", "\n")
      cat("Basic Information ", "\n")
      cat("-----------------", "\n")
      if (z$is.balanced) {
         cat("* Number of subgroups =", z$m, "\n")
         cat("* Sample size =",         z$nk,"\n")
      } else {
         cat("* Number of subgroups =", z$m, "\n")
         cat("* Subgroup sample sizes =", z$nsubgr, "\n")
      }
      cat("* Standard Error (SE) =", z$SE, "\n")
      cat("\n")

      # scale
      cat("=================", "\n")
      cat( paste0("Scale (", z$scale, ")"), "\n")
      cat("-----------------", "\n")
      estimate = round(z$sigma.hat,digits)
      if (z$is.balanced) {
         REpooled = round(100*suppressMessages(RE(n=z$nk, estimator=z$scale, correction=TRUE)),2)
         OUT = c(estimate, paste0(round(100*finite.breakdown(z$nk,z$scale),2)), REpooled)
         names(OUT) = c("Estimate", "Breakdown(%)", "RE(%)")
         print(OUT, quote=FALSE)
      } else {
         REpooled = suppressMessages(RE(n=z$nsubgr, estimator=z$scale, poolType=z$poolScale, correction=TRUE))
         subsample.sizes = sort(unique(z$nsubgr))
         REvalues <- BDvalues <- numeric(length(subsample.sizes))
         for ( i in seq_along(subsample.sizes) ) {
              BDvalues[i] = finite.breakdown(subsample.sizes[i],z$scale)
              REvalues[i] = suppressMessages(RE(n=subsample.sizes[i],estimator=z$scale,correction=TRUE))
         }

         cat( paste0("* Pooled estimate of ", z$scale,  " (Type=",z$poolScale,"): ", estimate), "\n")
         cat( paste0("* Minimum of subgroup breakdown points: ", round(100*min(BDvalues),2),"%"), "\n")
         cat( paste0("* RE of pooled ", z$scale,
                     " (Type=",z$poolScale,")"," with respect to mean", " (Type=",z$poolScale,"): ",
                       round(100*REpooled,2), "%"), "\n")
         cat("* Breakdown and RE for each subgroup sample size\n")
         OUT = rbind( round(100*BDvalues,2),  round(100*REvalues,2) )
         dimnames(OUT) = list( c("Breakdown (%)", "RE (%)"), paste0( "n=",subsample.sizes) )
         print(OUT, quote=FALSE)
      }
      cat("\n")

      # LCL, CL, UCL
      cat("=================", "\n")
      cat( paste0("Control Limits (with sample size = ", z$nk, "):"), "\n" )
      cat("-----------------", "\n")
      limits = c(z$LCL, z$CL, z$UCL)
      names(limits) = c("LCL", "CL", "UCL")
      print(limits, quote=FALSE)
      cat("\n")
   } else {
   # R chart
      TITLE = paste0(" R Control chart using scale (unbiased ", z$scale, ")",
                     " with ", z$sig, "*sigma." )
   ## cat("\n", "\033[1m", TITLE, "\033(B\033[m", "\n\n")
      cat("\n", TITLE, "\n\n")

      # Basic info. Sample info.  
      cat("=================", "\n")
      cat("Basic Information ", "\n")
      cat("-----------------", "\n")
      if (z$is.balanced) {
         cat("* Number of subgroups =", z$m, "\n")
         cat("* Sample size =",         z$nk,"\n")
      } else {
         cat("* Number of subgroups =", z$m, "\n")
         cat("* Subgroup sample sizes =", z$nsubgr, "\n")
      }
      cat("* Standard Error (SE) =", z$SE, "\n")
      cat("\n")

      # scale
      cat("=================", "\n")
      cat( paste0("Scale (", z$scale, ")"), "\n")
      cat("-----------------", "\n")
      estimate = round(z$sigma.hat,digits)
      if (z$is.balanced) {
         REpooled = round(100*suppressMessages(RE(n=z$nk, estimator=z$scale, correction=TRUE)),2)
         OUT = c(estimate, paste0(round(100*finite.breakdown(z$nk,z$scale),2)), REpooled)
         names(OUT) = c("Estimate", "Breakdown(%)", "RE(%)")
         print(OUT, quote=FALSE)
      } else {
         REpooled = suppressMessages(RE(n=z$nsubgr, estimator=z$scale, poolType=z$poolScale, correction=TRUE))
         subsample.sizes = sort(unique(z$nsubgr))
         REvalues <- BDvalues <- numeric(length(subsample.sizes))
         for ( i in seq_along(subsample.sizes) ) {
              BDvalues[i] = finite.breakdown(subsample.sizes[i],z$scale)
              REvalues[i] = suppressMessages(RE(n=subsample.sizes[i],estimator=z$scale,correction=TRUE))
         }

         cat( paste0("* Pooled estimate of ", z$scale,  " (Type=",z$poolScale,"): ", estimate), "\n")
         cat( paste0("* Minimum of subgroup breakdown points: ", round(100*min(BDvalues),2),"%"), "\n")
         cat( paste0("* RE of pooled ", z$scale,
                     " (Type=",z$poolScale,")"," with respect to mean", " (Type=",z$poolScale,"): ",
                       round(100*REpooled,2), "%"), "\n")
         cat("* Breakdown and RE of each subgroup sample size\n")
         OUT = rbind( round(100*BDvalues,2),  round(100*REvalues,2) )
         dimnames(OUT) = list( c("Breakdown (%)", "RE (%)"), paste0( "n=",subsample.sizes) )
         print(OUT, quote=FALSE)
      }
      cat("\n")

      # LCL, CL, UCL
      cat("=================", "\n")
      cat( paste0("Control Limits (with sample size = ", z$nk, "):"), "\n" )
      cat("-----------------", "\n")
      limits = c(z$LCL, z$CL, z$UCL)
      names(limits) = c("LCL", "CL", "UCL")
      print(limits, quote=FALSE)
      cat("\n")
   }
}
#============================================

#============================================
plot.rcc <- function(x, digits= getOption("digits")-2, 
            col.line="cyan4", col.background="ivory", col.boundary="ivory4", 
            cex.text=0.7, x.text, text.offset=0, LCL, CL, UCL, ...)
{  
   drawpolygon <- function(xrange, yrange, delta, col=col.background, border=FALSE) {
      polygon(c(xrange[1]-delta,xrange[1]-delta, xrange[2]+delta, xrange[2]+delta, xrange[1]-delta),
              c(yrange[1],yrange[2],yrange[2],yrange[1],yrange[1]), col=col, border=border)
   }
   if (missing(LCL)) LCL = x$LCL
   if (missing(CL) )  CL = x$CL
   if (missing(UCL)) UCL = x$UCL
   SE=x$SE;
   m = x$m
   subgr.index = seq_len(m)

   half = 0.5
   if (missing(x.text)) x.text = rep(m+half, 3) 

   par( mar=c(5,3,2,5) ) ## To write LCL, CL, UCL on the right side. (default: mar=c(5,4,4,2))
   # Basic setup for unbalanced case
   if (!x$is.balanced) {
       xx1 = c(1-half/4, rep(2:m,rep(2,m-1))-half, m+half/4)
       xx2 = rev(xx1)
       yy1 = rep(x$UCLsubgr, rep(2,m) )
       yy2 = rev( rep(x$LCLsubgr, rep(2,m)) )
       xxx = c(xx1,xx2,xx1[1])
       yyy = c(yy1,yy2,yy1[1])
   }

   # Control charts
   if (x$type == "Xbar") {  
      ## Xbar chart
      ## xbari = apply(x$data,MARGIN=1, FUN=mean)
      xbari = sapply(x$data, mean)  

      if (x$is.balanced) {
          ymin = min(xbari,LCL); ymax = max(xbari,UCL)

          index1 = (xbari >= LCL) & (xbari <= UCL)  ## in-control
          index2 = !index1                        ## out-of-control 

          plot(NA,NA, xlim=c(1, m), ylim=c(ymin-0.2*SE, ymax+0.2*SE),
               type="p", pch=1, frame=FALSE, axes=FALSE, ylab=NA, xlab="Subgroup" )

          drawpolygon(xrange=c(1,m), yrange=c(LCL, UCL), delta=half/4)
      } else {
          ymin = min(xbari,LCL,x$LCLsubgr); ymax = max(xbari,UCL,x$UCLsubgr)

          ## index1 = (xbari >= LCL) & (xbari <= UCL)             ## in-control
          index1 = (xbari >= x$LCLsubgr) & (xbari <= x$UCLsubgr)  ## in-control
          index2 = !index1                                      ## out-of-control 

          plot(NA,NA, xlim=c(1, m), ylim=c(ymin-0.2*SE, ymax+0.2*SE),
               type="p", pch=1, frame=FALSE, axes=FALSE, ylab=NA, xlab="Subgroup" )
          #for ( i in subgr.index ) {
          #    drawpolygon(xrange=c(i,i), yrange=c(x$LCLsubgr[i], x$UCLsubgr[i]))
          #    lines( c(i-0.5,i+0.5), c(x$LCLsubgr[i],x$LCLsubgr[i]), lwd=2, col=col.boundary )
          #    lines( c(i-0.5,i+0.5), c(x$UCLsubgr[i],x$UCLsubgr[i]), lwd=2, col=col.boundary )
          #}
          polygon(xxx,yyy, col=col.background, border=FALSE)
          lines(xx1,yy1, lwd=2, col=col.boundary)
          lines(xx2,yy2, lwd=2, col=col.boundary)
      }

      lines(  subgr.index, xbari, lty=1, col=col.line)
      points( subgr.index[index1], xbari[index1], pch=20, col="blue" )
      points( subgr.index[index2], xbari[index2], pch=4,  col="red" )

   } else if (x$type == "S") { 
      ## S chart
      ## si = apply(x$data,MARGIN=1, FUN=sd)
      si = sapply(x$data, sd)   # Don't call from res.

      if (x$is.balanced) {
          ymin = min(si,LCL); ymax = max(si,UCL)

          index1 = (si >= LCL) & (si <= UCL)  ## in-control
          index2 = !index1                    ## out-of-control 

          plot(NA,NA, xlim=c(1, m), ylim=c(ymin-0.2*SE, ymax+0.2*SE),
               type="p", pch=1, frame=FALSE, axes=FALSE, ylab=NA, xlab="Subgroup" )

          drawpolygon(xrange=c(1,m), yrange=c(LCL, UCL), delta=half/4)
      } else {
          ymin = min(si,LCL,x$LCLsubgr); ymax = max(si,UCL, x$UCLsubgr)
          
          ## index1 = (si >= LCL) & (si <= UCL)             ## in-control
          index1 = (si >= x$LCLsubgr) & (si <= x$UCLsubgr)  ## in-control
          index2 = !index1                                  ## out-of-control 
          
          plot(NA, NA, xlim=c(1, m), ylim=c(ymin-0.2*SE, ymax+0.2*SE),
               type="p", pch=1, frame=FALSE, axes=FALSE, ylab=NA, xlab="Subgroup" )
          
          #for ( i in subgr.index ) {
          #    drawpolygon(xrange=c(i,i), yrange=c(x$LCLsubgr[i], x$UCLsubgr[i]))
          #    lines( c(i-0.5,i+0.5), c(x$LCLsubgr[i],x$LCLsubgr[i]), lwd=2, col=col.boundary )
          #    lines( c(i-0.5,i+0.5), c(x$UCLsubgr[i],x$UCLsubgr[i]), lwd=2, col=col.boundary )
          #}
          polygon(xxx,yyy, col=col.background, border=FALSE)
          lines(xx1,yy1, lwd=2, col=col.boundary)
          lines(xx2,yy2, lwd=2, col=col.boundary)
      }

      lines(  subgr.index, si, lty=1, col=col.line)
      points( subgr.index[index1], si[index1], pch=20, col="blue" )
      points( subgr.index[index2], si[index2], pch=4,  col="red" )

   } else { 
      ## R chart
      ##Ri = apply(x$data,MARGIN=1, FUN=function(x){diff(range(x))}) 
      Ri = apply(sapply(x$data,range),2,diff)  # Don't call from res.

      if (x$is.balanced) {
          ymin = min(Ri,LCL); ymax = max(Ri,UCL)

          index1 = (Ri >= LCL) & (Ri <= UCL)              ## in-control
          index2 = !index1                                ## out-of-control 

          plot(NA, NA, xlim=c(1, m), ylim=c(ymin-0.2*SE, ymax+0.3*SE),
               type="p", pch=1, frame=FALSE, axes=FALSE, ylab=NA, xlab="Subgroup" )

          drawpolygon(xrange=c(1,m), yrange=c(LCL, UCL), delta=half/4)
      } else {
          ymin = min(Ri,LCL,x$LCLsubgr); ymax = max(Ri,UCL,x$UCLsubgr)

          index1 = (Ri >= x$LCLsubgr) & (Ri <= x$UCLsubgr)  ## in-control
          index2 = !index1                                  ## out-of-control 

          plot(NA,NA, xlim=c(1, m), ylim=c(ymin-0.2*SE, ymax+0.2*SE),
               type="p", pch=1, frame=FALSE, axes=FALSE, ylab=NA, xlab="Subgroup" )

          #for ( i in subgr.index ) {
          #    drawpolygon(xrange=c(i,i), yrange=c(x$LCLsubgr[i], x$UCLsubgr[i]))
          #    lines( c(i-0.5,i+0.5), c(x$LCLsubgr[i],x$LCLsubgr[i]), lwd=2, col=col.boundary )
          #    lines( c(i-0.5,i+0.5), c(x$UCLsubgr[i],x$UCLsubgr[i]), lwd=2, col=col.boundary )
          #}
          polygon(xxx,yyy, col=col.background, border=FALSE)
          lines(xx1,yy1, lwd=2, col=col.boundary)
          lines(xx2,yy2, lwd=2, col=col.boundary)
      }

      lines(  subgr.index, Ri, lty=1, col=col.line)
      points( subgr.index[index1], Ri[index1], pch=20, col="blue" )
      points( subgr.index[index2], Ri[index2], pch=4,  col="red" )
   }

   offset = text.offset*cex.text*SE
   text(x.text[1], LCL-offset, paste0("LCL=", round(LCL,digits)), xpd=NA, adj=0, cex=cex.text)
   text(x.text[2],  CL+offset, paste0(" CL=", round(CL,digits)),  xpd=NA, adj=0, cex=cex.text)
   text(x.text[3], UCL+offset, paste0("UCL=", round(UCL,digits)), xpd=NA, adj=0, cex=cex.text)

   ## abline( h= CL, lty=4, lwd=0.5 )
   ## abline( h=c(LCL,UCL), lty=2, lwd=1.0 )
   lines( c(1-half/4,m+half/4), c( CL, CL),  lty=1, lwd=1.0, col="gold4" )
   lines( c(1-half/4,m+half/4), c(LCL,LCL),  lty=2, lwd=1.0 )
   lines( c(1-half/4,m+half/4), c(UCL,UCL),  lty=2, lwd=1.0 )

   axis(1, at =subgr.index, subgr.index, lty=1, lwd=0.5)
   axis(2, lty=1)
}
#============================================

