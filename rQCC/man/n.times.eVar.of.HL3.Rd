\name{n.times.eVar.of.HL3}
\alias{n.times.eVar.of.HL3}
\title{Empirical variances of the Hodges-Lehmann (HL3) estimator (times n)}

\description{
n times the empirical variances of the Hodges-Lehmann estimator (HL3) 
under the standard normal distribution, 
where \eqn{n}{n} is the sample size and \eqn{n}{n} is from 1 to 100.
For the HL3, the median is taken over all \eqn{(i,j)}{(i,j)} 
where \eqn{i,j=1,2,\ldots,n}{i, j=1,2,...,n},
that is,
\deqn{\mathrm{HL3} = \mathop{\mathrm{median}}_{\forall(i,j)} \Big( \frac{X_i+X_j}{2} \Big).}{HL3 = median of (Xi+Xj)/2 over all (i,j).}
They are obtained using the extensive Monte Carlo simulation
with 1E07 replicates. 
}

\usage{n.times.eVar.of.HL3}
\format{This data frame contains 100 values.}

\value{\code{n.times.eVar.of.HL3} returns a vector of 100 values.} 

\references{
Park, C., H. Kim, and M. Wang (2019).
Finite-sample properties of robust location and scale estimators. 
\emph{arXiv:1908.00462}.


Hodges, J. L. and E. L. Lehmann (1963).
Estimates of location based on rank tests.
\emph{Annals of Mathematical Statistics}, \bold{34}, 598--611. \cr 
}

\author{Chanseok Park and Min Wang}

\keyword{datasets}
\keyword{variance}

