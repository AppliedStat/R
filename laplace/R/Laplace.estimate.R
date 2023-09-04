
#------------------------------
laplace.mle <- function(x) {
   loc  = median(x)
   scale= mean( abs(x-loc) )
   structure(list(location=loc, scale=scale),class="laplace.estimate")
}
#------------------------------
laplace.lp <- function(x,n, a){   ## It was a=0.5
  x = sort(x)
  r = length(x)
  if ( missing(n) ) n=r
  if (missing(a)) { a = ifelse(n <= 10, 3/8, 1/2) }
  if ( n < r ) stop("n cannot be smaller than the number of observations")
  R = ppoints(n,a=a)
  p = R[1:r]
  p.minus.5 = p-0.5 
  d = sign(p.minus.5)
  y = -d*log(1-2*abs(p.minus.5))
  LM=lm(y ~ x)
  B = as.numeric( coef (LM) )
  scale = 1/B[2]
  loc   = -scale*B[1]
  structure(list(location=loc, scale=scale), class="laplace.estimate")
}
#------------------------------
laplace.rm <- function(x, a){ 
    x = sort(unique(x))
    n = length(x)
    if (missing(a)) {
        a = ifelse(n <= 10, 3/8, 1/2)
    }
    p = ppoints(n, a = a)
    p.minus.5 = p-0.5 
    d = sign(p.minus.5)
    y = -d*log(1-2*abs(p.minus.5))
    DX = outer(x, x, "-")
    DY = outer(y, y, "-")
    diag(DX) = NA
    slope.hat = median(apply(DY/DX, 2, median, na.rm = TRUE))
    intercept = median(y - slope.hat * x)
    loc = - slope.hat * intercept
    structure(list(location = loc, scale = 1/slope.hat ), class = "laplace.estimate")
}
#------------------------------
print.laplace.estimate <-
function(x, digits = getOption("digits"), ...)
{ 
   ans = format(x, digits=digits)
   dn  = dimnames(ans)
   print(ans, quote=FALSE)
   invisible(x)
}
#------------------------------------------------------------------
