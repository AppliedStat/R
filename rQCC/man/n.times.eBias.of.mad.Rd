\name{n.times.eBias.of.mad}
\alias{n.times.eBias.of.mad}
\title{Empirical biases of the MAD (median absolute deviation) estimator (times n)}

\description{
n times the empirical biases of the median absolute deviation(MAD) estimator under the standard
normal distribution, where n is the sample size and n is from 1 to 100.
For the median absolute deviation(MAD) estimates, mad\{stats\} is used. 
They are obtained using the extensive Monte Carlo simulation 
with 1E07 replicates. 
}

\usage{n.times.eBias.of.mad}
\format{This data frame contains 100 values.}

\value{\code{n.times.eBias.of.mad} returns a vector of 100 values.}

\references{
Park, C., H. Kim, and M. Wang (2019).
Finite-sample properties of robust location and scale estimators. 
\emph{arXiv:1908.00462}.
}



\author{Chanseok Park and Min Wang}
\keyword{datasets}
\keyword{bias}

