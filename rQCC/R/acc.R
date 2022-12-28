#============================================
# Attribute Control Chart
#--------------------------------------------
acc <- function(x, n, type=c("p", "np", "c", "u" ,"g", "h", "t"), parameter, 
                pEstimator=c("Wald","Wilson"), gEstimator=c("ML","MVU"), tModel=c("E","W"), 
                location.shift=0, sigmaFactor=3, nk) 
{
  if (is.matrix(x)) x <- unname( split(x, row(x)) )
        type = match.arg(type)
  pEstimator = match.arg(pEstimator)
  gEstimator = match.arg(gEstimator)
      tModel = match.arg(tModel)

  m = length(x)

  if (!missing(n))  {
     is.balanced = (length(unique(n))==1L)
     if (length(n)==1L) n = rep(n, length.out=m)
     # if (length(n) != m) stop("The number of sample sizes does not match.")
     N = sum(n)
     if (missing(nk))  nk = N/m
  } else {
     if ( any(type %in% c("p","np", "u")) ) stop("n is missing.")
  } 

  g  = sigmaFactor
  g2 = g^2 

  parameter.missing = missing(parameter)
  #--------------------------------------------
  if (type == "p") {
  #  p chart
     ## pbar= ifelse (parameter.missing, sum(x)/N, parameter)
     if (parameter.missing) {
        pbar = sum(x)/N
        parameter = pbar
     } else {
        pbar = parameter
     }
     CL  = switch(pEstimator, Wald = pbar, Wilson = (pbar+g2/(2*N))/(1+g2/N) )
     SE  = switch(pEstimator,
                  Wald  = sqrt(pbar*(1-pbar)/nk),
                  Wilson= 1/(1+g2/N)*sqrt( pbar*(1-pbar)/nk + g2/(4*nk*N) ) )
     UCL = min( CL + g*SE, 1)
     LCL = max( CL - g*SE, 0)

     # Subgroup
      SEsubgr = switch(pEstimator, 
                 Wald  = sqrt(pbar*(1-pbar)/n), 
                 Wilson= 1/(1+g2/N)*sqrt( pbar*(1-pbar)/n + g2/(4*n*N) ) )
     UCLsubgr = pmin( CL + g*SEsubgr, 1)
     LCLsubgr = pmax( CL - g*SEsubgr, 0)
     res <- list(LCL=LCL, CL=CL, UCL=UCL, LCLsubgr=LCLsubgr, UCLsubgr=UCLsubgr, 
                 data=x, n=n, m=m, nk=nk, is.balanced=is.balanced, SE=SE, 
                 type=type, pEstimator=pEstimator, parameter=parameter, 
                 parameter.missing=parameter.missing, sig=sigmaFactor)
     class(res) <- "acc"
     return(res)
  } else if (type == "np") {
  #  np chart
     if (!is.balanced) stop("The np chart is not appropriate with unbalanced samples.")
     ## pbar= ifelse (missing(parameter), sum(x)/N, parameter)
     if (parameter.missing) {
        pbar = sum(x)/N
        parameter = pbar
     } else {
        pbar = parameter
     }
     CL  = switch(pEstimator, Wald = pbar, Wilson = (pbar+g2/(2*N))/(1+g2/N) )
     SE  = switch(pEstimator, 
                  Wald  = sqrt(pbar*(1-pbar)/nk), 
                  Wilson= 1/(1+g2/N)*sqrt( pbar*(1-pbar)/nk + g2/(4*nk*N) ) )
     UCL = nk*min(CL+g*SE, 1)  
     LCL = nk*max(CL-g*SE, 0)  
     CL  = nk*CL
    
     res <- list(LCL=LCL, CL=CL, UCL=UCL, 
                 data=x, n=n, m=m, nk=nk, is.balanced=is.balanced, SE=SE, 
                 type=type, pEstimator=pEstimator, parameter=parameter, 
                 parameter.missing=parameter.missing, sig=sigmaFactor)

     class(res) <- "acc"
     return(res)
  } else if (type == "c") {
  # c chart
     is.balanced = TRUE
     # lam= ifelse (missing(parameter), mean(x), parameter)
     if (parameter.missing) {
        lam = mean(x) 
        parameter = lam
     } else {
        lam = parameter
     }

     CL  = lam
     SE  = sqrt(lam) 
     UCL = lam + g*SE 
     LCL = max(lam - g*SE, 0)

     res <- list(LCL=LCL, CL=CL, UCL=UCL, data=x, m=m, is.balanced=is.balanced, SE=SE,
                 type=type, parameter=parameter,
                 parameter.missing=parameter.missing, sig=sigmaFactor)
     class(res) <- "acc"
     return(res)
  } else if (type == "u") {
  # u chart
     ## lam= ifelse (missing(parameter), sum(x)/N, parameter)
     if (parameter.missing) {
        lam = sum(x)/N 
        parameter = lam
     } else {
        lam = parameter
     }

     CL  = lam 
     SE = sqrt(lam/nk) 
     UCL = lam + g*SE 
     LCL = max(lam - g*SE, 0) 
    
     # Subgroup
      SEsubgr = sqrt(lam/n) 
     UCLsubgr = CL + g*SEsubgr 
     LCLsubgr = pmax( CL - g*SEsubgr, 0)
     res <- list(LCL=LCL, CL=CL, UCL=UCL, LCLsubgr=LCLsubgr, UCLsubgr=UCLsubgr,
                 data=x, n=n, m=m, nk=nk, is.balanced=is.balanced, SE=SE,
                 type=type, pEstimator=pEstimator, parameter=parameter, 
                 parameter.missing=parameter.missing, sig=sigmaFactor)
     class(res) <- "acc"
     return(res)
  } else if (type == "g") {  
  # g chart
     a = location.shift
     N = length(unlist(x)) 
     if (missing(nk)) nk = N/m
     ni = sapply(x, length)
     is.balanced = (length(unique(ni))==1L)

     if (parameter.missing) {
        xbarbar = mean(unlist(x)) 
        parameter = switch(gEstimator, ML = 1/(xbarbar-a+1), MVU = (1-1/N)/(xbarbar-a+1-1/N) )
        SE0 = sqrt(nk*(xbarbar-a)*(xbarbar-a+1))
        SE  = switch(gEstimator, ML = SE0, MVU = SE0*sqrt(N/(N+1)) )
        CL  = nk*xbarbar
        UCL = CL + g*SE
        LCL = max(CL - g*SE, a*nk)
       
        SEsub0 = sqrt(ni*(xbarbar-a)*(xbarbar-a+1))
        SEsub  = switch(gEstimator, ML = SEsub0, MVU = SEsub0*sqrt(N/(N+1)) )
        CLsubgr  = ni*xbarbar
        UCLsubgr = CLsubgr + g*SEsub 
        LCLsubgr = pmax(CLsubgr - g*SEsub, a*ni)
     } else {
        p = parameter 
        SE = sqrt(nk*(1-p)/(p^2))
        CL  = nk*(1/p-1+a)
        UCL = CL + g*SE
        LCL = max(CL - g*SE, a*nk)

        SEsub = sqrt(ni*(1-p)/(p^2))
        CLsubgr  = ni*(1/p-1+a)
        UCLsubgr = CLsubgr + g*SEsub
        LCLsubgr = pmax(CLsubgr - g*SEsub, a*ni)
     } 
     res <- list(LCL=LCL, CL=CL, UCL=UCL, LCLsubgr=LCLsubgr, CLsubgr=CLsubgr, UCLsubgr=UCLsubgr,
                 data=x, m=m, nk=nk, is.balanced=is.balanced, SE=SE,
                 type=type, gEstimator=gEstimator, parameter=parameter, 
                 parameter.missing=parameter.missing, sig=sigmaFactor)
     class(res) <- "acc"
     return(res)
  } else if (type == "h") { 
  # h chart
     a = location.shift
     N = length(unlist(x))
     if (missing(nk)) nk = N/m
     ni = sapply(x, length)
     is.balanced = (length(unique(ni))==1L)

     if (parameter.missing) {
        xbarbar = mean(unlist(x))
        parameter = switch(gEstimator, ML = 1/(xbarbar-a+1), MVU = (1-1/N)/(xbarbar-a+1-1/N) )
        SE0 = sqrt((xbarbar-a)*(xbarbar-a+1)/nk)
        SE  = switch(gEstimator, ML = SE0, MVU = SE0*sqrt(N/(N+1)) )
        CL  = xbarbar
        UCL = CL + g*SE
        LCL = max(CL - g*SE, a)

        SEsub0 = sqrt((xbarbar-a)*(xbarbar-a+1)/ni)
        SEsub  = switch(gEstimator, ML = SEsub0, MVU = SEsub0*sqrt(N/(N+1)) )
        CLsubgr  = xbarbar
        UCLsubgr = CLsubgr + g*SEsub
        LCLsubgr = pmax(CLsubgr - g*SEsub, a)
     } else {
        p = parameter
        SE = sqrt((1-p)/(nk*p^2))
        CL  = 1/p-1+a
        UCL = CL + g*SE
        LCL = max(CL - g*SE, a)

        SEsub = sqrt((1-p)/(p^2)/ni)
        CLsubgr  = 1/p-1+a
        UCLsubgr = CLsubgr + g*SEsub
        LCLsubgr = pmax(CLsubgr - g*SEsub, a)
     } 
     res <- list(LCL=LCL, CL=CL, UCL=UCL, LCLsubgr=LCLsubgr, CLsubgr=CLsubgr, UCLsubgr=UCLsubgr,
                 data=x, m=m, nk=nk, is.balanced=is.balanced, SE=SE,
                 type=type, gEstimator=gEstimator, parameter=parameter, 
                 parameter.missing=parameter.missing, sig=sigmaFactor)
     class(res) <- "acc"
     return(res)
  } else if (type == "t")  { 
  # t chart (exponential or Weibull)
     is.balanced = TRUE
     a.over.2 = pnorm(sigmaFactor,lower.tail=FALSE)

     # MLE of Weibull
     weibull.MLE <- function (x, tol=.Machine$double.eps^0.25, maxiter=1000, trace=0) {   
       # Setup for interval 
         meanlog = mean(log(x))
         lower = 1/(log(max(x)) - meanlog)
         upper = sum((x^lower) * log(x))/sum(x^lower) - meanlog
         interval = c(lower, 1/upper)
       # EE equation 
         EEweibull = function(alpha, x) {
             xalpha = x^alpha
             sum(log(x) * (xalpha))/sum(xalpha) - 1/alpha - mean(log(x))
         }
         tmp = uniroot(EEweibull, interval = interval, x = x, tol = tol,
             maxiter = maxiter, trace = trace)
         alpha = tmp$root
         beta = mean(x^alpha)^(1/alpha)
         return( c(alpha,beta) )
     } # END of MLE of Weibull

     if (parameter.missing) {
        if (tModel == "E") {
           parameter = mean(x) 
        } else if ( tModel == "W") {
           parameter = weibull.MLE(x)
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
         SE = parameter[2] * ( gamma(1+2/parameter[1]) - gamma(1+1/parameter[1])^2 )
        UCL = qweibull(a.over.2, shape=parameter[1], scale=parameter[2], lower.tail=FALSE)
        LCL = qweibull(a.over.2, shape=parameter[1], scale=parameter[2], lower.tail=TRUE)
     }
     res <- list(LCL=LCL, CL=CL, UCL=UCL, data=x, m=m, is.balanced=is.balanced, SE=SE,
                 type=type, parameter=parameter, tModel=tModel,
                 parameter.missing=parameter.missing, sig=sigmaFactor)
     class(res) <- "acc"
     return(res)
  }  else {
     stop("Type should be either one of \"p\", \"np\", \"c\", \"u\", \"g\", \"h\", and \"t\".")
  }
}
#============================================
print.acc <- function(x, ...) { 
   limits = c(x$LCL, x$CL, x$UCL)
   names(limits) = c("LCL", "CL", "UCL")
   print(limits, quote=FALSE)
   invisible(x)
}
#============================================
summary.acc <- function(object, ...) {
   z=object
   digits = getOption("digits") 
   subgr.index = seq_len(z$m)

   if (z$type == "p") {
   # p chart
      if (z$pEstimator == "Wald") {
          TITLE = paste0("The conventional p Control Chart with ", z$sig, "*sigma." )
      } else {
          TITLE = paste0("The p Control Chart based on the Wilson CI with ", z$sig, "*sigma." )
      } 
   ## cat("\n", "\033[1m", TITLE, "\033(B\033[m", "\n\n")
      cat("\n", TITLE, "\n\n")

      # Basic info. Sample info. 
      cat("=================", "\n")
      cat("Basic Information ", "\n")
      cat("-----------------", "\n")
      if (z$is.balanced) {
         OUT = c(paste(c(z$m, z$nk)))
         names(OUT) = c("Number of subgroups", "Sample size")
         print(OUT, quote=FALSE, print.gap=1)
      } else {
         cat("* Number of subgroups =", z$m, "\n")
         cat("* Subgroup sample sizes =", z$n, "\n")
      }
      tmp = ifelse(z$parameter.missing, "* Estimated parameter =", "* Parameter value =")
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
   } else if (z$type == "np") {
   # np chart
      if (!z$is.balanced) stop("The np chart is not appropriate with unbalanced samples.")
      if (z$pEstimator == "Wald") {
          TITLE = paste0("The conventional np Control Chart with ", z$sig, "*sigma." )
      } else {
          TITLE = paste0("The np Control Chart based on the Wilson CI with ", z$sig, "*sigma." )
      }
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
   } else if (z$type == "c") {
   # c chart
      TITLE = paste0("The c Control Chart with ", z$sig, "*sigma." )
   ## cat("\n", "\033[1m", TITLE, "\033(B\033[m", "\n\n")
      cat("\n", TITLE, "\n\n")

      # Basic info. Sample info.  
      cat("=================", "\n")
      cat("Basic Information ", "\n")
      cat("-----------------", "\n")
      cat("* Number of subgroups =", z$m, "\n")
      tmp = ifelse(z$parameter.missing, "* Parameter estimate =", "* Parameter value =")
      cat(tmp, z$parameter, "\n")
      cat("* Standard error (SE) =", z$SE, "\n")
      cat("\n")

      # LCL, CL, UCL
      cat("=================", "\n")
      cat( "Control Limits ",  "\n" )
      cat("-----------------", "\n")
      limits = c(z$LCL, z$CL, z$UCL)
      names(limits) = c("LCL", "CL", "UCL")
      print(limits, quote=FALSE)
      cat("\n")
   } else if (z$type == "u") {
   # u chart
      TITLE = paste0("The u Control Chart with ", z$sig, "*sigma." )
   ## cat("\n", "\033[1m", TITLE, "\033(B\033[m", "\n\n")
      cat("\n", TITLE, "\n\n")

      # Basic info. Sample info.  
      cat("=================", "\n")
      cat("Basic Information ", "\n")
      cat("-----------------", "\n")
      if (z$is.balanced) {
         OUT = c(paste(c(z$m, z$nk)))
         names(OUT) = c("Number of subgroups", "Sample size")
         print(OUT, quote=FALSE, print.gap=1)
      } else {
         cat("* Number of subgroups =", z$m, "\n")
         cat("* Subgroup sample sizes =", z$n, "\n")
      }
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
   } else if (z$type == "g") {
   # g chart
      TITLE = paste0("The g Control Chart with ", z$sig, "*sigma", " based on ", z$gEstimator, " method.")
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
      TITLE = paste0("The h Control Chart with ", z$sig, "*sigma", " based on ", z$gEstimator, " method.")
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
   } else if ( (z$type=="t")&&(z$tModel=="E") ) {
   # Exponential t chart
      TITLE = paste0("The exponential t Control Chart with ", z$sig, "*sigma." )
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
   } else if ( (z$type=="t")&&(z$tModel=="W") ) {
   # Weibull t chart
      TITLE = paste0("The Weibull t Control Chart with ", z$sig, "*sigma." )
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
     stop("Type should be either one of \"p\", \"np\", \"c\", \"u\", \"g\", and \"h\".")
   }
}
#============================================


#============================================
plot.acc <- function(x, digits= getOption("digits")-2, 
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
   if (x$type == "p") {  
      ## p chart
      pbari = x$data / x$n   
      if (x$is.balanced) {
          ymin = min(pbari,LCL); ymax = max(pbari,UCL)

          index1 = (pbari >= LCL) & (pbari <= UCL)  ## in-control
          index2 = !index1                          ## out-of-control 

          plot(NA,NA, xlim=c(1, m), ylim=c( max(ymin-0.5*SE,0), min(ymax+0.5*SE,1) ),
               type="p", pch=1, frame=FALSE, axes=FALSE, ylab=NA, xlab="Subgroup" )

          drawpolygon(xrange=c(1,m), yrange=c(LCL, UCL), delta=half/4)
      } else {
          ymin = min(pbari,LCL,x$LCLsubgr); ymax = max(pbari,UCL,x$UCLsubgr)

          ## index1 = (pbari >= LCL) & (pbari <= UCL)             ## in-control
          index1 = (pbari >= x$LCLsubgr) & (pbari <= x$UCLsubgr)  ## in-control
          index2 = !index1                                        ## out-of-control 

          plot(NA,NA, xlim=c(1, m), ylim=c( max(ymin-0.5*SE,0), min(ymax+0.5*SE,1) ),
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

      lines(  subgr.index, pbari, lty=1, col=col.line)
      points( subgr.index[index1], pbari[index1], pch=20, col="blue" )
      points( subgr.index[index2], pbari[index2], pch=4,  col="red" )
   } else if (x$type == "np") { 
      ## np chart
      if (!x$is.balanced) stop("The np chart is not appropriate with unbalanced samples.")
      npbari = x$data 
      ymin = min(npbari,LCL)
      ymax = max(npbari,UCL)
  
      index1 = (npbari >= LCL) & (npbari <= UCL)  ## in-control
      index2 = !index1                            ## out-of-control 

      plot(NA,NA, xlim=c(1, m), ylim=c(max(ymin-0.5*SE,0), ymax+0.5*SE),
           type="p", pch=1, frame=FALSE, axes=FALSE, ylab=NA, xlab="Subgroup" )

      drawpolygon(xrange=c(1,m), yrange=c(LCL, UCL), delta=half/4)

      lines(  subgr.index, npbari, lty=1, col=col.line)
      points( subgr.index[index1], npbari[index1], pch=20, col="blue" )
      points( subgr.index[index2], npbari[index2], pch=4,  col="red" )
   } else if (x$type == "c") { 
      ## c chart
      lami = x$data 
      ymin = min(lami,LCL)
      ymax = max(lami,UCL)

      index1 = (lami >= LCL) & (lami <= UCL)  ## in-control
      index2 = !index1                        ## out-of-control 

      plot(NA,NA, xlim=c(1, m), ylim=c( max(ymin-0.5*SE,0), min(ymax+0.5*SE) ),
           type="p", pch=1, frame=FALSE, axes=FALSE, ylab=NA, xlab="Subgroup" )

      drawpolygon(xrange=c(1,m), yrange=c(LCL, UCL), delta=half/4)

      lines(  subgr.index, lami, lty=1, col=col.line)
      points( subgr.index[index1], lami[index1], pch=20, col="blue" )
      points( subgr.index[index2], lami[index2], pch=4,  col="red" )
   } else if (x$type == "u") { 
      ## u chart
      lami = x$data / x$n
      if (x$is.balanced) {
          ymin = min(lami,LCL); ymax = max(lami,UCL)

          index1 = (lami >= LCL) & (lami <= UCL)  ## in-control
          index2 = !index1                        ## out-of-control 

          plot(NA,NA, xlim=c(1, m), ylim=c( max(ymin-0.5*SE,0), min(ymax+0.5*SE) ),
               type="p", pch=1, frame=FALSE, axes=FALSE, ylab=NA, xlab="Subgroup" )

          drawpolygon(xrange=c(1,m), yrange=c(LCL, UCL), delta=half/4)
      } else {
          ymin = min(lami,LCL,x$LCLsubgr); ymax = max(lami,UCL,x$UCLsubgr)

          ## index1 = (lami >= LCL) & (lami <= UCL)             ## in-control
          index1 = (lami >= x$LCLsubgr) & (lami <= x$UCLsubgr)  ## in-control
          index2 = !index1                                        ## out-of-control 

          plot(NA,NA, xlim=c(1, m), ylim=c( max(ymin-0.5*SE,0), min(ymax+0.5*SE) ),
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

      lines(  subgr.index, lami, lty=1, col=col.line)
      points( subgr.index[index1], lami[index1], pch=20, col="blue" )
      points( subgr.index[index2], lami[index2], pch=4,  col="red" )
   } else if (x$type == "g") { 
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
      stop("Type should be either one of \"p\", \"np\", \"c\", \"u\", \"g\", \"h\", and \"t\".")
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

