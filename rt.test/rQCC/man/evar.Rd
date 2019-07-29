\name{empirical.variance}
\alias{evar}
\title{Empirical variances of robust estimators}
\usage{
evar (n, method = c("median", "HL1", "HL2", "HL3", "mad", "shamos")) 
}
\arguments{
 \item{n}{sample size (\eqn{n \ge 1}{n \ge 1}).}
 \item{method}{a character string specifying the estimator, must be 
one of "median" (default), "HL1", "HL2", "HL3", "mad", and "shamos."}
}

\description{Provides the empirical variances of the median, Hodges-Lehmann (HL1, HL2, HL3) 
median absolute deviation (MAD) and Shamos estimators.
}

\details{
The empirical variances for \eqn{n=1,2,\ldots,10}{n = 1, 2, ... ,100} are obtained
using the extensive Monte Carlo simulation with 1E07 replicates.
For the case of \eqn{n > 100}{n > 100}, they are obtained using the method of Hayes (2014).
}

\value{\code{evar} computes the empirical varinace of a specific 
estimator (one of "median", "HL1", "HL2", "HL3", "mad", and "shamos") 
when a random sample is from the standard normal distribution.}

\seealso{
\code{rQCC::RE} for the relative efficiency. \cr 
\code{rQCC::n.times.eVar.of.HL1} for the empirical variance of the HL1 estimator (times n). \cr
\code{rQCC::n.times.eVar.of.HL2} for the empirical variance of the HL2 estimator (times n). \cr
\code{rQCC::n.times.eVar.of.HL3} for the empirical variance of the HL3 estimator (times n). \cr
\code{rQCC::n.times.eVar.of.mad} for the empirical variance of the MAD estimator (times n). \cr
\code{rQCC::n.times.eVar.of.median} for the empirical variance of the median estimator (times n). \cr
\code{rQCC::n.times.eVar.of.shamos} for the empirical variance of the Shamos estimator (times n). \cr
}

\references{
Park, C., H. Kim and M. Wang (2019).
Investigation of finite-sample properties of robust location
and scale estimators.  \emph{Submitted for Publication}.

Hayes, K. (2014).
Finite-sample bias-correction factors for the median absolute deviation.
\emph{Communications in Statistics: Simulation and Computation},
\bold{43}, 2205--2212.
}


\examples{
# Empirical variance of the Hodges-Lehmann estimator (HL2) under the standard normal distribution.
evar (n=10, method="HL2")
}

\keyword{variance}
\author{Chanseok Park and Min Wang}
