#------------------------------
# We need n when the data are right-censored at the max. obs
gumbel.gp <- function(x,n, a){   ## It was a=0.5
  x = sort(x)
  r = length(x)
  if ( missing(n) ) n=r
  if (missing(a)) { a = ifelse(n <= 10, 3/8, 1/2) }
  if ( n < r ) stop("n cannot be smaller than the number of observations")
  R = ppoints(n,a=a)
  LM=lm(-log(-log( R[1:r] ))~x)
  B = as.numeric( coef (LM) )
  beta = 1/B[2]
  mu   = -beta*B[1]
  structure(list(location=mu, scale=beta), class="gumbel.estimate")
}
#------------------------------
print.gumbel.estimate <-
function(x, digits = getOption("digits"), ...)
{ 
   ans = format(x, digits=digits)
   dn  = dimnames(ans)
   print(ans, quote=FALSE)
   invisible(x)
}
#------------------------------


