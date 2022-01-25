\encoding{utf-8}
\name{Hodges-Lehmann}
\alias{HL}
\title{Hodges-Lehmann estimator}

\description{Calculates the Hodges-Lehmann estimator.}

\usage{
HL(x, estimator = c("HL1", "HL2", "HL3"), na.rm = FALSE)
}

\arguments{
 \item{x}{a numeric vector of observations.}
 \item{estimator}{a character string specifying the estimator, must be 
one of \code{"HL1"} (default), \code{"HL2"} and \code{"HL3"}.}
 \item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
}


\details{
\code{HL} computes the Hodges-Lehmann estimators (one of \code{"HL1"}, \code{"HL2"}, \code{"HL3"}).

The Hodges-Lehmann (HL1) is defined as 
\deqn{\mathrm{HL1} = \mathop{\mathrm{median}}_{i<j} \Big( \frac{X_i+X_j}{2} \Big)}{HL1 = median of 
(Xi+Xj)/2 over i<j}  
where \eqn{i,j=1,2,\ldots,n}{i, j=1,2,...,n}.

The Hodges-Lehmann (HL2) is defined as 
\deqn{\mathrm{HL2} = \mathop{\mathrm{median}}_{i \le j}\Big(\frac{X_i+X_j}{2} \Big).}{HL2 = median of
(Xi+Xj)/2 over i \le j.}

The Hodges-Lehmann (HL3) is defined as 
\deqn{\mathrm{HL3} = \mathop{\mathrm{median}}_{\forall(i,j)} \Big( \frac{X_i+X_j}{2} \Big).}{HL3 = 
median of (Xi+Xj)/2 over all (i,j).}
}

\value{It returns a numeric value.}

\references{
Park, C., H. Kim, and M. Wang (2022).
Investigation of finite-sample properties of robust location and scale estimators.
\emph{Communications in Statistics - Simulation and Computation}, To appear.\cr
\doi{10.1080/03610918.2019.1699114}

Hodges, J. L. and E. L. Lehmann (1963).
Estimates of location based on rank tests.
\emph{Annals of Mathematical Statistics}, \bold{34}, 598--611.
}

\seealso{
\link{mean}\{base\} for calculating sample mean and 
\link{median}\{stats\} for calculating sample median. 

\code{finite.breakdown}\{rQCC\} for calculating the finite-sample breakdown point.
}

\examples{
x = c(0:10, 50)
HL(x, estimator="HL2")
}

\author{Chanseok Park and Min Wang}
\keyword{robust}
\keyword{univar}
