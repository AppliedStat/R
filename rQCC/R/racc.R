#============================================
# Attribute Control Chart
#--------------------------------------------
racc <- function(x, gamma, type=c("g","h","t"), parameter, 
                gEstimator=c("cdf","MM"), tModel=c("E","W"),
                location.shift=0, sigmaFactor=3, nk) 
{
  if (is.matrix(x)) x <- unname( split(x, row(x)) )
        type = match.arg(type)
  gEstimator = match.arg(gEstimator)
      tModel = match.arg(tModel)

  m = length(x)

  g  = sigmaFactor
  g2 = g^2 

  ## parameter.missing = missing(parameter)
  parameter.missing = missing(parameter)

  TINY = sqrt(.Machine$double.eps)
  #--------------------------------------------
  a = location.shift
  N = length(unlist(x)) 
  if (missing(nk)) nk = N/m
  ni = sapply(x, length)
  is.balanced = (length(unique(ni))==1L)

  if (parameter.missing) {
     if ( (type=="g") | (type=="h") ) {
         xunlist =  unlist(x)
         xbarbar = mean(xunlist) 
         if (gEstimator == "cdf") { # cdf method
            d1 = max(0,quantile(xunlist, gamma/2))
            d2 = max(a,quantile(xunlist, gamma) - d1 + a - 1)
            F1  = sum( xunlist < d1+TINY )
            F2  = sum( xunlist < d2+TINY )
            F12 = sum( xunlist < d1+d2+1-a+TINY )
            phat = 1 - ((F12-F1)/F2)^(1/(d1+1-a))
            if ( phat < TINY )   phat=TINY
            if ( phat > 1-TINY ) phat=1-TINY
            parameter = phat
         } else {                  # MM method
            d = max(a, quantile(xunlist, gamma))
            y = xunlist[xunlist <= d]
            n = length(y)
            ybar = mean(y)
            s2 = var(y) * (n-1)/n
            if (ybar >= (a+d)/2) d = floor(2*ybar-a) + 1
          # phat = (d-a+2-2*(ybar-a+1)) / ((d-a+1)*(ybar-a+1)-(ybar-a+1)^2 - s2)
            phat = (d+a-2*ybar) / ((d-ybar)*(ybar-a+1)-s2)
            if ( phat < TINY )   phat=TINY
            if ( phat > 1-TINY ) phat=1-TINY
            parameter = phat
         }
     }
  } 
           
  if (type == "g") {  
  # g chart
     p = parameter 
     SE = sqrt(nk*(1-p)/(p^2))
     CL  = nk*(1/p-1+a)
     UCL = CL + g*SE
     LCL = max(CL - g*SE, a*nk)

     SEsub = sqrt(ni*(1-p)/(p^2))
     CLsubgr  = ni*(1/p-1+a)
     UCLsubgr = CLsubgr + g*SEsub
     LCLsubgr = pmax(CLsubgr - g*SEsub, a*ni)

     res <- list(LCL=LCL, CL=CL, UCL=UCL, LCLsubgr=LCLsubgr, CLsubgr=CLsubgr, UCLsubgr=UCLsubgr,
                 data=x, m=m, nk=nk, is.balanced=is.balanced, SE=SE,
                 type=type, gEstimator=gEstimator, parameter=parameter, 
                 parameter.missing=parameter.missing, sig=sigmaFactor)
     class(res) <- "racc"
     return(res)
  } else if (type == "h") { 
  # h chart
     p = parameter
     SE = sqrt((1-p)/(nk*p^2))
     CL  = 1/p-1+a
     UCL = CL + g*SE
     LCL = max(CL - g*SE, a)

     SEsub = sqrt((1-p)/(p^2)/ni)
     CLsubgr  = 1/p-1+a
     UCLsubgr = CLsubgr + g*SEsub
     LCLsubgr = pmax(CLsubgr - g*SEsub, a)

     res <- list(LCL=LCL, CL=CL, UCL=UCL, LCLsubgr=LCLsubgr, CLsubgr=CLsubgr, UCLsubgr=UCLsubgr,
                 data=x, m=m, nk=nk, is.balanced=is.balanced, SE=SE,
                 type=type, gEstimator=gEstimator, parameter=parameter, 
                 parameter.missing=parameter.missing, sig=sigmaFactor)
     class(res) <- "racc"
     return(res)
  } else if (type == "t")  {
  # t chart (exponential or Weibull)   
     is.balanced = TRUE
     a.over.2 = pnorm(sigmaFactor,lower.tail=FALSE)

     #------------------------------------------------------
##############################################################################################
     weibull.rm <- function(x, a){   ## Repeated median (added on 2022-07-10)
       x = sort( log(unique(x)) )
       n = length(x)
       if (missing(a)) { a = ifelse(n <= 10, 3/8, 1/2) }
       y = log(-log(1-ppoints(n, a=a)))
       DX = outer(x,x,"-"); DY = outer(y,y,"-")
       diag(DX) = NA
       beta.hat = median( apply(DY/DX,2, median, na.rm=TRUE) )  # slope 
       intercept = median(y-beta.hat*x) # This intercept is better than the below. 
       ## intercept = median( apply((outer(x,y,"*")-outer(y,x,"*"))/DX, 2, median, na.rm=TRUE) )
       theta.hat = exp(-intercept/beta.hat)
       return( c(beta.hat, theta.hat) )
     }
#------------------------------
     if (parameter.missing) {
        if (tModel == "E") {
           ## parameter = mean(x)                             # not robust
           parameter = median( -sort(x)/log(1-ppoints(m)) )   # robust version
        } else if ( tModel == "W") {
       ## parameter = weibull.MLE(x)                          # MLE version  
          parameter = weibull.rm(x)                           # robust version
        } else {
          stop("tModel should be one of \"E\" or \"W\".")
        }
     } 
     # LCL and UCL for t chart 
     if (tModel == "E") {
         CL = qexp(0.5, rate=1/parameter)
         SE = parameter
        UCL = qexp(a.over.2, rate=1/parameter, lower.tail=FALSE)
        LCL = qexp(a.over.2, rate=1/parameter, lower.tail=TRUE)
     } else  {
         CL = qweibull(0.5, shape=parameter[1], scale=parameter[2])
         SE = parameter[2] * ( base::gamma(1+2/parameter[1]) - base::gamma(1+1/parameter[1])^2 )
        UCL = qweibull(a.over.2, shape=parameter[1], scale=parameter[2], lower.tail=FALSE)
        LCL = qweibull(a.over.2, shape=parameter[1], scale=parameter[2], lower.tail=TRUE)
     }
     res <- list(LCL=LCL, CL=CL, UCL=UCL, data=x, m=m, is.balanced=is.balanced, SE=SE,
                 type=type, parameter=parameter, tModel=tModel,
                 parameter.missing=parameter.missing, sig=sigmaFactor)
     class(res) <- "racc"
     return(res)
  } else {
     stop("Type should be either one of \"g\", \"h\", and \"t\".")
  }
}
#============================================
print.racc <- function(x, ...) { 
   limits = c(x$LCL, x$CL, x$UCL)
   names(limits) = c("LCL", "CL", "UCL")
   print(limits, quote=FALSE)
   invisible(x)
}
#============================================
summary.racc <- function(object, ...) {
   z=object
   digits = getOption("digits") 
   subgr.index = seq_len(z$m)

   if (z$type == "g") {
   # g chart
      TITLE = paste0("The robust g Control Chart with ", z$sig, "*sigma", " based on ", z$gEstimator, " method.")
   ## cat("\n", "\033[1m", TITLE, "\033(B\033[m", "\n\n")
      cat("\n", TITLE, "\n\n")

      # Basic info. Sample info. 
      cat("=================", "\n")
      cat("Basic Information ", "\n")
      cat("-----------------", "\n")

      OUT = c(paste(c(z$m, z$nk)))
      names(OUT) = c("Number of subgroups", "Sample size")
      print(OUT, quote=FALSE, print.gap=1)
      tmp = ifelse(z$parameter.missing, "* Parameter estimate =", "* Parameter value =")
      cat(tmp, z$parameter, "\n")
      cat("* Standard error (SE) =", z$SE, "\n")
      cat("\n")

      # LCL, CL, UCL
      cat("=================", "\n")
      cat( paste0("Control Limits (with sample size = ", z$nk, "):"), "\n" )
      cat("-----------------", "\n")
      limits = c(z$LCL, z$CL, z$UCL)
      names(limits) = c("LCL", "CL", "UCL")
      print(limits, quote=FALSE)
      cat("\n")
   } else if (z$type == "h") {
   # h chart
      TITLE = paste0("The robust h Control Chart with ", z$sig, "*sigma", " based on ", z$gEstimator, " method.")
   ## cat("\n", "\033[1m", TITLE, "\033(B\033[m", "\n\n")
      cat("\n", TITLE, "\n\n")

      # Basic info. Sample info. 
      cat("=================", "\n")
      cat("Basic Information ", "\n")
      cat("-----------------", "\n")

      OUT = c(paste(c(z$m, z$nk)))
      names(OUT) = c("Number of subgroups", "Sample size")
      print(OUT, quote=FALSE, print.gap=1)
      tmp = ifelse(z$parameter.missing, "* Parameter estimate =", "* Parameter value =")
      cat(tmp, z$parameter, "\n")
      cat("* Standard error (SE) =", z$SE, "\n")
      cat("\n")

      # LCL, CL, UCL
      cat("=================", "\n")
      cat( paste0("Control Limits (with sample size = ", z$nk, "):"), "\n" )
      cat("-----------------", "\n")
      limits = c(z$LCL, z$CL, z$UCL)
      names(limits) = c("LCL", "CL", "UCL")
      print(limits, quote=FALSE)
      cat("\n")
   } else if ( (z$type=="t")&&(z$tModel=="E") ) {   # 
   # Exponential t chart
      TITLE = paste0("The robust exponential t Control Chart with ", z$sig, "*sigma." )
   ## cat("\n", "\033[1m", TITLE, "\033(B\033[m", "\n\n")
      cat("\n", TITLE, "\n\n")

      # Basic info. Sample info.  
      cat("=================",  "\n")
      cat("Basic Information ", "\n")
      cat("-----------------",  "\n")
      cat("* Number of observations =", z$m, "\n")
      tmp = ifelse(z$parameter.missing, "* Parameter estimate (mean) =", "* Parameter value (mean) =")
      cat(tmp, z$parameter, "\n")
      cat("\n")

      # LCL, CL, UCL
      cat("=================", "\n")
      cat( "Control Limits:",  "\n")
      cat("-----------------", "\n")
      limits = c(z$LCL, z$CL, z$UCL)
      names(limits) = c("LCL", "CL", "UCL")
      print(limits, quote=FALSE)
      cat("\n")
   } else if ( (z$type=="t")&&(z$tModel=="W") ) {   # 
   # Weibull t chart
      TITLE = paste0("The robust Weibull t Control Chart with ", z$sig, "*sigma." )
   ## cat("\n", "\033[1m", TITLE, "\033(B\033[m", "\n\n")
      cat("\n", TITLE, "\n\n")

      # Basic info. Sample info.  
      cat("=================",  "\n")
      cat("Basic Information ", "\n")
      cat("-----------------",  "\n")
      cat("* Number of observations =", z$m, "\n")

      if(z$parameter.missing) {
         cat("* Parameter estimate:", "shape =", z$parameter[1], " scale =", z$parameter[2])
      } else {
         cat("* Parameter value:", "shape =", z$parameter[1], " scale =", z$parameter[2])
      }
      cat("\n")

      # LCL, CL, UCL
      cat("=================", "\n")
      cat("Control Limits:",  "\n" )
      cat("-----------------", "\n")
      limits = c(z$LCL, z$CL, z$UCL)
      names(limits) = c("LCL", "CL", "UCL")
      print(limits, quote=FALSE)
      cat("\n")
   } else {
     stop("Type should be either one of \"g\", \"h\", and \"t\".")
   }
}
#============================================


#============================================
plot.racc <- function(x, digits= getOption("digits")-2, 
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
   SE=x$SE; SEsubgr=x$SEsubgr
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

   # Control Charts
   if (x$type == "g") { 
      ## g chart
      xbari = sapply(x$data, sum)
      if (x$is.balanced) {
          ymin = min(xbari,LCL); ymax = max(xbari,UCL)

          index1 = (xbari >= LCL) & (xbari <= UCL)  ## in-control
          index2 = !index1                          ## out-of-control 

          plot(NA,NA, xlim=c(1, m), ylim=c(ymin-0.2*SE, ymax+0.2*SE),
               type="p", pch=1, frame=FALSE, axes=FALSE, ylab=NA, xlab="Subgroup" )

          drawpolygon(xrange=c(1,m), yrange=c(LCL, UCL), delta=half/4)
      } else {
          ymin = min(xbari,LCL,x$LCLsubgr); ymax = max(xbari,UCL,x$UCLsubgr)

          ## index1 = (xbari >= LCL) & (xbari <= UCL)             ## in-control
          index1 = (xbari >= x$LCLsubgr) & (xbari <= x$UCLsubgr)  ## in-control
          index2 = !index1                                        ## out-of-control 

          plot(NA,NA, xlim=c(1, m), ylim=c(ymin-0.2*SE, ymax+0.2*SE),
               type="p", pch=1, frame=FALSE, axes=FALSE, ylab=NA, xlab="Subgroup" )
          polygon(xxx,yyy, col=col.background, border=FALSE)
          lines(xx1,yy1, lwd=2, col=col.boundary)
          lines(xx2,yy2, lwd=2, col=col.boundary)
      }

      lines(  subgr.index, xbari, lty=1, col=col.line)
      points( subgr.index[index1], xbari[index1], pch=20, col="blue" )
      points( subgr.index[index2], xbari[index2], pch=4,  col="red" )
   } else if (x$type == "h") { 
      ## h chart
      xbari = sapply(x$data, mean)
      if (x$is.balanced) {
          ymin = min(xbari,LCL); ymax = max(xbari,UCL)
          
          index1 = (xbari >= LCL) & (xbari <= UCL)  ## in-control
          index2 = !index1                          ## out-of-control 
          
          plot(NA,NA, xlim=c(1, m), ylim=c(ymin-0.2*SE, ymax+0.2*SE),
               type="p", pch=1, frame=FALSE, axes=FALSE, ylab=NA, xlab="Subgroup" )
          
          drawpolygon(xrange=c(1,m), yrange=c(LCL, UCL), delta=half/4)
      } else {
          ymin = min(xbari,LCL,x$LCLsubgr); ymax = max(xbari,UCL,x$UCLsubgr)

          ## index1 = (xbari >= LCL) & (xbari <= UCL)             ## in-control
          index1 = (xbari >= x$LCLsubgr) & (xbari <= x$UCLsubgr)  ## in-control
          index2 = !index1                                        ## out-of-control 

          plot(NA,NA, xlim=c(1, m), ylim=c(ymin-0.2*SE, ymax+0.2*SE),
               type="p", pch=1, frame=FALSE, axes=FALSE, ylab=NA, xlab="Subgroup" )
          polygon(xxx,yyy, col=col.background, border=FALSE)
          lines(xx1,yy1, lwd=2, col=col.boundary)
          lines(xx2,yy2, lwd=2, col=col.boundary)
      }

      lines(  subgr.index, xbari, lty=1, col=col.line)
      points( subgr.index[index1], xbari[index1], pch=20, col="blue" )
      points( subgr.index[index2], xbari[index2], pch=4,  col="red" )
   } else if (x$type == "t") {
      ## t chart
      xi = x$data
      ymin = min(xi,LCL)
      ymax = max(xi,UCL)

      index1 = (xi >= LCL) & (xi <= UCL)  ## in-control
      index2 = !index1                    ## out-of-control 

      plot(NA,NA, xlim=c(1, m), ylim=c( max(ymin,0), min(ymax) ),
           type="p", pch=1, frame=FALSE, axes=FALSE, ylab=NA, xlab="Observations" )

      drawpolygon(xrange=c(1,m), yrange=c(LCL, UCL), delta=half/4)

      lines(  subgr.index, xi, lty=1, col=col.line)
      points( subgr.index[index1], xi[index1], pch=20, col="blue" )
      points( subgr.index[index2], xi[index2], pch=4,  col="red" )
   } else { 
      stop("Type should be either one of \"g\", \"h\", and \"t\".")
   }

   ## Write CL, LCL, UCL on the right margin
   offset = text.offset*cex.text*SE*2 
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

