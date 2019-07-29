\name{Hodges-Lehmann}
\alias{HL}
\title{Hodges-Lehmann estimate}
\usage{
HL(x, method = c("HL1", "HL2", "HL3"), na.rm = FALSE)
}
\arguments{
 \item{x}{a numeric vector of observations.}
 \item{method}{a character string specifying the estimator, must be 
one of "HL1" (default), "HL2" and "HL3."}
 \item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
}

\description{Calculates the Hodges-Lehmann estimate.}

\details{
The Hodges-Lehmann (HL1) is defined as 
\deqn{\mathrm{HL1} = \mathop{\mathrm{median}}_{i<j} \Big( \frac{X_i+X_j}{2} \Big).}{HL1 = median of 
(Xi+Xj)/2 over i<j}  
where \eqn{i,j=1,2,\ldots,n}{i, j=1,2,...,n}.

The Hodges-Lehmann (HL2) is defined as 
\deqn{\mathrm{HL2} = \mathop{\mathrm{median}}_{i \le j}\Big(\frac{X_i+X_j}{2} \Big).}{HL2 = median of
(Xi+Xj)/2 over i \le j.}

The Hodges-Lehmann (HL3) is defined as 
\deqn{\mathrm{HL3} = \mathop{\mathrm{median}}_{\forall(i,j)} \Big( \frac{X_i+X_j}{2} \Big).}{HL3 = 
median of (Xi+Xj)/2 over all (i,j).}
}


\references{
Park, C., H. Kim and M. Wang (2019).
Investigation of finite-sample properties of robust location
and scale estimators.  \emph{Submitted for Publication}.

Hodges, J. L. and E. L. Lehmann (1963).
Estimates of location based on rank tests.
\emph{Annals of Mathematical Statistics}, \bold{34}, 598--611.
}

\seealso{
\code{base::\link{mean}} for calculating sample mean and 
\code{stats::\link{median}} for calculating sample median.  \cr

\code{rQCC::finite.breakdown} for calculating the finite-sample breakdown point.
}

\examples{
x = c(0:10, 50)
HL(x, method="HL2")
}
%% \keyword{print}

\author{Chanseok Park and Min Wang}

\keyword{breakdown}
\keyword{robust}

