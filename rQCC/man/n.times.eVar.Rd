\encoding{utf-8}
\name{n.times.eVar}
\alias{n.times.eVar.of.HL1}
\alias{n.times.eVar.of.HL2}
\alias{n.times.eVar.of.HL3}
\alias{n.times.eVar.of.mad}
\alias{n.times.eVar.of.median}
\alias{n.times.eVar.of.shamos}

\title{Empirical variances (times n)}

\description{\eqn{n}{n} times the empirical variances of the Hodges-Lehmann (HL1, HL2, HL3), 
the median, the median absolute deviation (MAD), and the Shamos estimators.
}

\usage{
n.times.eVar.of.HL1

n.times.eVar.of.HL2

n.times.eVar.of.HL3

n.times.eVar.of.median

n.times.eVar.of.mad

n.times.eVar.of.shamos
}

\details{\eqn{n}{n} times the empirical variances of the Hodges-Lehmann (HL1, HL2, HL3), the median,
the median absolute deviation (MAD), and the Shamos estimators
under the standard normal distribution,
where \eqn{n}{n} is the sample size and \eqn{n}{n} is from 1 to 100.


For the MAD estimators, \link{mad}\{stats\} is used. 
For the Hodges-Lehmann, \code{HL}\{rQCC\} is used.
For the Shamos, the Fisher-consistent Shamos estimator, \code{shamos}\{rQCC\}, is used.

The empirical variances  are obtained using the Monte Carlo simulation
with 1E07 replicates. 
}

\value{They return a vector of 100 values.} 

\references{
Park, C., H. Kim, and M. Wang (2022).
Investigation of finite-sample properties of robust location and scale estimators.
\emph{Communications in Statistics - Simulation and Computation}, To appear.\cr
\url{https://doi.org/10.1080/03610918.2019.1699114}

Hodges, J. L. and E. L. Lehmann (1963).
Estimates of location based on rank tests.
\emph{Annals of Mathematical Statistics}, \bold{34}, 598--611.

Shamos, M. I. (1976).  Geometry and statistics: Problems at the interface.
In Traub, J. F., editor,
\emph{Algorithms and Complexity: New Directions and Recent Results},
pages 251--280. Academic Press, New York.

Lèvy-Leduc, C., Boistard, H., Moulines, E., Taqqu, M. S.,      
and Reisen, V. A.  (2011).
Large sample behaviour of some well-known robust estimators under
long-range dependence.
\emph{Statistics}, \bold{45}, 59--71.
}

\author{Chanseok Park and Min Wang}
\keyword{datasets}
\keyword{variance}

