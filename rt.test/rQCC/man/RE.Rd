\name{relative.efficiency}
\alias{RE}
\title{Relative efficiency (RE)}
\usage{
RE (n, method = c("median", "HL1", "HL2", "HL3", "mad", "shamos")) 
}
\arguments{
 \item{n}{sample size (\eqn{n \ge 1}{n \ge 1}).}
 \item{method}{a character string specifying the estimator, must be 
one of "median" (default), "HL1", "HL2", "HL3", "mad", and "shamos."}
}

\description{Calculates the relative efficiency of 
the median and Hodges-Lehmann (HL1, HL2, HL3) estimators with respect to the sample mean
and that of the median absolute deviation (MAD) and Shamos estimators with respect 
to the sample standard deviation.

The relative efficiencies for \eqn{n=1,2,\ldots,10}{n = 1, 2, ... ,100} are obtained
using the extensive Monte Carlo simulation with 1E07 replicates.
For the case of \eqn{n > 100}{n > 100}, they are obtained based on  the method of Hayes (2014).
}

\details{
The relative efficiency (RE) of \eqn{\hat{\theta}_2}{\theta2_hat} with respect to 
\eqn{\hat{\theta}_1}{\theta1_hat} is defined as 
\deqn{\mathrm{RE}(\hat{\theta}_2 | \hat{\theta}_1) 
=\frac{\mathrm{Var}(\hat{\theta}_1)}{\mathrm{Var}(\hat{\theta}_2)}.}{RE(\theta2_hat|\theta1_hat) = Var(\theta1_hat) divided by Var(\theta2_hat).}
}

\value{\code{RE} calculates the relative efficiency
of a location estimator ("median", "HL1", "HL2", "HL3") with respect to the sample mean and 
that of a scale estimator ("mad", "shamos") with respect to the sample standard deviation.}

\seealso{
\code{rQCC::n.times.eVar.of.HL1} for the empirical variance of the HL1 estimator (times n).\cr
\code{rQCC::n.times.eVar.of.HL2} for the empirical variance of the HL2 estimator (times n).\cr
\code{rQCC::n.times.eVar.of.HL3} for the empirical variance of the HL3 estimator (times n).\cr
\code{rQCC::n.times.eVar.of.mad} for the empirical variance of the MAD estimator (times n).\cr
\code{rQCC::n.times.eVar.of.median} for the empirical variance of the median estimator (times n).\cr
\code{rQCC::n.times.eVar.of.shamos} for the empirical variance of the Shamos estimator (times n).\cr
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
# The relative efficiency of the Hodges-Lehmann (HL2) estimator
#     with respect to the sample standard deviation under the normal distribution.

RE(n=10, method="HL2")
}

\keyword{variance}
\author{Chanseok Park and Min Wang}
