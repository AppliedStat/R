\name{n.times.eVar.of.mad}
\alias{n.times.eVar.of.mad}
\title{Empirical variances of the MAD (median absolute deviation) estimator (times n)}

\description{
n times the empirical variances of the median absolute deviation(MAD) estimator under the standard
normal distribution, 
where \eqn{n}{n} is the sample size and \eqn{n}{n} is from 1 to 100.
For the median absolute deviation(MAD) estimates, mad\{stats\} is used. 
They are obtained using the extensive Monte Carlo simulation
with 1E07 replicates. 
}

\usage{n.times.eVar.of.mad}
\format{This data frame contains 100 values.}

\value{\code{n.times.eVar.of.mad} returns a numeric value.}


\references{
Park, C., H. Kim and M. Wang (2019).
Investigation of finite-sample properties of robust location 
and scale estimators.  \emph{Submitted for Publication}.
}

\author{Chanseok Park and Min Wang}

\keyword{datasets}
\keyword{variance}
