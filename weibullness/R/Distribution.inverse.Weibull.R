# 
dinvweibull = function(x, shape, scale = 1, log = FALSE) {
    k <- max(lx<-length(x),lshape<-length(shape),lscale<-length(scale))
    if (lx     < k)     x <- rep(x,     length=k) 
    if (lshape < k) shape <- rep(shape, length=k) 
    if (lscale < k) scale <- rep(scale, length=k) 
    #---
    logd = numeric(k)
    id = (x > 0)
    logd[id]  = log(shape[id])-log(scale[id]) - (shape[id]+1)*log(x[id]/scale[id]) - 1/(x[id]/scale[id])^shape[id]
    logd[!id] = -Inf
    #---
    if(!is.null(Names <- names(x))) names(logd) = rep(Names,length=k)
    if (log) { 
       return(logd) 
    } else {
       return( exp(logd) )
    }
}

#
pinvweibull = function(q, shape, scale=1, lower.tail=TRUE, log.p=FALSE) {
    k <- max(lq<-length(q),lshape<-length(shape),lscale<-length(scale))
    if (lq     < k)     q <- rep(q,     length=k)
    if (lshape < k) shape <- rep(shape, length=k)
    if (lscale < k) scale <- rep(scale, length=k)
    #---
    id = which(q > 0)
    p = numeric(length(q))
    p[id] = pweibull(1/q[id], shape=shape[id], scale=1/scale[id], lower.tail=!lower.tail, log.p=log.p)
    if(!is.null(Names <- names(q))) names(p) = rep(Names,length=k)
    return(p)
}

#
qinvweibull = function(p, shape, scale=1, lower.tail = TRUE, log.p = FALSE) {
    k <- max(lp<-length(p),lshape<-length(shape),lscale<-length(scale))
    if (lp     < k)     p <- rep(p,     length=k)
    if (lshape < k) shape <- rep(shape, length=k)
    if (lscale < k) scale <- rep(scale, length=k)
    q = 1/qweibull(p = p, shape = shape, scale = 1/scale, lower.tail = !lower.tail, log.p = log.p)
    if(!is.null(Names<-names(p))) names(q)<-rep(Names,length=k)
    return(q)
}

#
rinvweibull = function(n, shape, scale=1) {
   lshape<-length(shape); lscale<-length(scale)
   if (lshape < n) shape <- rep(shape, length=n)
   if (lscale < n) scale <- rep(scale, length=n)
   1/rweibull(n, shape=shape, scale=1/scale)
}
#


