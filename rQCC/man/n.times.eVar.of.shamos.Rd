\name{n.times.eVar.of.shamos}
\alias{n.times.eVar.of.shamos}
\encoding{utf-8}
\title{Empirical variances of the Shamos estimator (times n)}

\description{
n times the empirical variances of the Shamos estimator under the standard
normal distribution, 
where \eqn{n}{n} is the sample size and \eqn{n}{n} is from 1 to 100.
Note that the Fisher-consistent Shamos estimator is used. 
The empirical variances are obtained 
using the extensive Monte Carlo simulation with 1E07 replicates. 
}

\usage{n.times.eVar.of.shamos}
\format{This data frame contains 100 values.}

\value{\code{n.times.eVar.of.shamos} returns a vector of 100 values.}

\references{
Park, C., H. Kim, and M. Wang (2019).
Finite-sample properties of robust location and scale estimators. 
\emph{arXiv:1908.00462}.

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
