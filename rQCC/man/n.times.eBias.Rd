\encoding{utf-8}
\name{n.times.eBias}
\alias{n.times.eBias.of.mad}
\alias{n.times.eBias.of.shamos}
\title{Empirical biases (times n)}

\description{\eqn{n}{n} times the empirical biases of the median absolute deviation (MAD) 
and Shamos estimators.
}

\usage{
n.times.eBias.of.mad

n.times.eBias.of.shamos
}

\details{\eqn{n}{n} times the empirical biases of the median absolute deviation (MAD) 
and Shamos estimators  
under the standard normal distribution, where \eqn{n}{n} is the sample size
and \eqn{n}{n} is from 1 to 100.
For the MAD estimator, \link{mad}\{stats\} is used.  
For the Shamos estimator, the Fisher-consistent Shamos estimator, \code{shamos}\{rQCC\}, is used.

These estimators are not unbiased with a finite sample.
The empirical biases are obtained using the extensive Monte Carlo simulation
with 1E07 replicates.
}

\value{They return a vector of 100 values.}

\references{
Park, C., H. Kim, and M. Wang (2022).
Investigation of finite-sample properties of robust location and scale estimators.
\emph{Communications in Statistics - Simulation and Computation}, 
\bold{51}, 2619-2645. \cr 
\doi{10.1080/03610918.2019.1699114}

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
\keyword{bias}
