\name{finite.sample.breakdown}
\alias{finite.breakdown}
\title{Finite-sample breakdown point of various estimators}
\usage{
finite.breakdown(n, method=c("median","HL1","HL2","HL3","mad","shamos")) 
}
\arguments{
 \item{n}{sample size (\eqn{n \ge 1}{n \ge 1}).}
 \item{method}{a character string specifying the estimator, must be 
one of "median" (default), "HL1", "HL2", "HL3", "mad", and "shamos."}
}


\description{Calculates finite-sample breakdown point of 
median, Hodges-Lehmann estimators (HL1, HL2, HL3), 
MAD (median absolute deviation) and Shamos.}

\details{
The Hodges-Lehmann (HL1) is defined as 
\deqn{\mathrm{HL1} = \mathop{\mathrm{median}}_{i<j} \Big( \frac{X_i+X_j}{2} \Big).}{HL1 = median of (X
i+Xj)/2 over i<j}  
where \eqn{i,j=1,2,\ldots,n}{i, j=1,2,...,n}.

The Hodges-Lehmann (HL2) is defined as 
\deqn{\mathrm{HL2} = \mathop{\mathrm{median}}_{i \le j}\Big(\frac{X_i+X_j}{2} \Big).}{HL2 = median of (Xi+Xj)/2 over i \le j.}

The Hodges-Lehmann (HL3) is defined as 
\deqn{\mathrm{HL3} = \mathop{\mathrm{median}}_{\forall(i,j)} \Big( \frac{X_i+X_j}{2} \Big).}{HL3 = median of (Xi+Xj)/2 over all (i,j).}
}


\value{\code{finite.breakdown} gives the finite-sample breakdown point
of the specified estimator.}

\references{
Park, C., H. Kim and M. Wang (2019).
Investigation of finite-sample properties of robust location
and scale estimators.  \emph{Submitted for Publication}.

Hodges, Jr., J. L. (1967).
Efficiency in normal samples and tolerance of extreme values 
for some estimates of location.
\emph{Proceedings of the Fifth Berkeley Symposium on Mathematical
Statistics and Probability}, Vol. \bold{1}, 163--186, 
Berkeley. University of California Press.

Hampel, F. R., Ronchetti, E., Rousseeuw, P. J., and Stahel, W. A. (1986).
\emph{Robust Statistics: The Approach Based on Influence Functions},
Subsection 2.2a.  John Wiley & Sons, New York.
}


\seealso{
\code{rQCC::HL} for the Hodges-Lehmann estimate.
}

\examples{
# finite-sample breakdown point of the Hodges-Lehmann (HL1) with size n=10.
finite.breakdown(n=10, method="HL2")
}


\author{Chanseok Park and Min Wang}
\keyword{breakdown}
\keyword{robust}
