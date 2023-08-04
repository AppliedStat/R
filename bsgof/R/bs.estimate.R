#####################################################################


# =============================================
# MLE of Birnbaum-Saunders 
# See Engelhardt, Bain, Wright. Tech. V.23, pp.251-256
# ---------------
bs.mle <- function(x) {
   r = 1/ mean(1/x)
   s =  mean(x)
   # lexical scoping of R (so, this won't work in S/Splus). 
   EE.mle = function(beta) {
      k = 1 / mean( 1/(beta+x) )
      return (beta^2 - beta*(2*r+k) + r*(s+k))
   }
   tmp = uniroot(EE.mle, interval=c(r,s))
   beta = tmp$root
   alpha = sqrt ( s/beta + beta/r - 2 )
   # list(alpha=alpha, beta=beta)
    structure(list(alpha=alpha,beta=beta),class="bs.estimate")

}
#---------------------
# MM of bs. Ng et al (2002). CSDA 43, 283-298
# ---------------
bs.mme <- function(x) {
   r = 1/ mean(1/x)
   s =  mean(x)
   beta = sqrt(s*r)
   alpha = sqrt( 2*sqrt(s/r) - 2 )
    structure(list(alpha=alpha,beta=beta),class="bs.estimate")
}
#------------------------------
print.bs.estimate <-
function(x, digits = getOption("digits"), ...)
{  
   ans = format(x, digits=digits)
   dn  = dimnames(ans)
   print(ans, quote=FALSE)
   invisible(x)
}
#------------------------------------------------------------------


#####################################################################
# pdf of Birnbaum-Saunders distribution 
dbs <- function(x, alpha=1, beta=1, log=FALSE)  {
   y = sqrt(x/beta)
   if (log==FALSE) { 
      (0.5/x) * (y+1/y) * dnorm(y-1/y, sd=alpha)
   } else { 
      log(0.5)-log(x) + log(y+1/y) + dnorm(y-1/y, sd=alpha, log=TRUE)
   }
}
# CDF of Birnbaum-Saunders distribution 
pbs <- function(q, alpha=1, beta=1, lower.tail=TRUE, log.p=FALSE) {
   y = sqrt(q/beta)
   pnorm(y-1/y, mean=0, sd=alpha, lower.tail=lower.tail, log.p=log.p)
}
# random number of Birnbaum-Saunders distribution 
rbs <- function(n, alpha=1, beta=1)  {
   x = rnorm(n, sd=alpha/2)
   x2 = x*x
   beta*(1+2*x2+2*x*sqrt(1+x2))
}
# quantile of Birnbaum-Saunders distribution 
qbs <- function(p, alpha=1, beta=1, lower.tail=TRUE, log.p=FALSE) {
   az = alpha*qnorm(p, mean=0, sd=1, lower.tail=lower.tail, log.p=log.p)
   .25*beta*(az + sqrt(az*az+4))^2
}
# # test the above
# p0 = runif(5); a0 = runif(5)*10; b0=runif(5)*10
# q = qbs(p0, a0,b0)
# cbind(p0, pbs(q, a0,b0) )
# #
# q0 = runif(5)*20; a0 = runif(5)*10; b0=runif(5)*10
# p1 = pbs(q0, a0, b0)
# cbind(q0, qbs(p1,a0, b0) )
#====================================================================


