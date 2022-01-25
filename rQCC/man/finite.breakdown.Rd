\encoding{utf-8}
\name{finite.sample.breakdown}
\alias{finite.breakdown}
\title{Finite-sample breakdown point}
       
\description{Calculates the finite-sample breakdown point of the mean,
median, Hodges-Lehmann estimators (HL1, HL2, HL3), standard deviation, range, 
MAD (median absolute deviation) and Shamos estimators.
Note that for the case of the mean, standard deviation and  range, 
the finite-sample breakdown points
are always zero.}

\usage{
finite.breakdown (n, estimator=c("mean","median","HL1","HL2","HL3",
                                 "sd","range","mad","shamos") )
}
\arguments{
 \item{n}{a numeric vector of sample sizes.}
 \item{estimator}{a character string specifying the estimator, must be 
one of \code{"mean"} (default), \code{"median"}, \code{"HL1"}, \code{"HL2"}, \code{"HL3"}, 
\code{"sd"}, \code{"range"}, \code{"mad"}, and \code{"shamos"}.}
}

\details{
\code{finite.breakdown} gives the finite-sample breakdown point
of the specified estimator.

The Hodges-Lehmann (HL1) is defined as 
\deqn{\mathrm{HL1} = \mathop{\mathrm{median}}_{i<j} \Big( \frac{X_i+X_j}{2} \Big)}{HL1 = median of (X
i+Xj)/2 over i<j}  
where \eqn{i,j=1,2,\ldots,n}{i, j=1,2,...,n}.

The Hodges-Lehmann (HL2) is defined as 
\deqn{\mathrm{HL2} = \mathop{\mathrm{median}}_{i \le j}\Big(\frac{X_i+X_j}{2} \Big).}{HL2 = median of (Xi+Xj)/2 over i \le j.}

The Hodges-Lehmann (HL3) is defined as 
\deqn{\mathrm{HL3} = \mathop{\mathrm{median}}_{\forall(i,j)} \Big( \frac{X_i+X_j}{2} \Big).}{HL3 = median of (Xi+Xj)/2 over all (i,j).}
}

\value{It returns a numeric value.}

\references{
Park, C., H. Kim, and M. Wang (2022).
Investigation of finite-sample properties of robust location and scale estimators.
\emph{Communications in Statistics - Simulation and Computation}, To appear.\cr
\doi{10.1080/03610918.2019.1699114}

Hodges, Jr., J. L. (1967).
Efficiency in normal samples and tolerance of extreme values 
for some estimators of location.
\emph{Proceedings of the Fifth Berkeley Symposium on Mathematical
Statistics and Probability}, Vol. \bold{1}, 163--186.
University of California Press, Berkeley.

Hampel, F. R., Ronchetti, E., Rousseeuw, P. J., and Stahel, W. A. (1986).
\emph{Robust Statistics: The Approach Based on Influence Functions},
Subsection 2.2a.  John Wiley & Sons, New York.
}

\seealso{
\code{HL}\{rQCC\} for the Hodges-Lehmann estimator.
}

\examples{
# finite-sample breakdown point of the Hodges-Lehmann (HL1) with size n=10.
finite.breakdown(n=10, estimator="HL2")

# finite-sample breakdown points of the median with sizes n=4,5,6 
finite.breakdown(n=4:6, estimator="median")
}


\author{Chanseok Park and Min Wang}
\keyword{breakdown}
\keyword{robust}
