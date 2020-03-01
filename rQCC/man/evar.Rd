\encoding{utf-8}
\name{empirical.variance}
\alias{evar}
\title{Empirical variances of robust estimators}
\usage{
evar (n, method = c("median", "HL1", "HL2", "HL3", "mad", "shamos")) 
}
\arguments{
 \item{n}{sample size (\eqn{n \ge 1}{n \ge 1}).}
 \item{method}{a character string specifying the estimator, must be 
one of \code{"median"} (default), \code{"HL1"}, \code{"HL2"}, \code{"HL3"}, \code{"mad"}, 
and \code{"shamos"}.}
}

\description{Provides the empirical variances of the median, Hodges-Lehmann (HL1, HL2, HL3) 
median absolute deviation (MAD) and Shamos estimators.
}

\details{
The empirical variances for \eqn{n=1,2,\ldots,100}{n = 1, 2, ... ,100} are obtained
using the extensive Monte Carlo simulation with 1E07 replicates.
For the case of \eqn{n > 100}{n > 100}, they are obtained using the method of Hayes (2014).

\code{evar} computes the empirical variance of a specific 
estimator (one of \code{"median"}, \code{"HL1"}, \code{"HL2"}, \code{"HL3"}, \code{"mad"}, 
and \code{"shamos"}) 
when a random sample is from the standard normal distribution.
}

\value{It returns a numeric value.}

\seealso{
\code{RE}\{rQCC\} for the relative efficiency. \cr 
\code{n.times.eVar.of.HL1}\{rQCC\} for the empirical variance of the HL1 estimator (times n). \cr
\code{n.times.eVar.of.HL2}\{rQCC\} for the empirical variance of the HL2 estimator (times n). \cr
\code{n.times.eVar.of.HL3}\{rQCC\} for the empirical variance of the HL3 estimator (times n). \cr
\code{n.times.eVar.of.mad}\{rQCC\} for the empirical variance of the MAD estimator (times n). \cr
\code{n.times.eVar.of.median}\{rQCC\} for the empirical variance of the median estimator (times n). \cr
\code{n.times.eVar.of.shamos}\{rQCC\} for the empirical variance of the Shamos estimator (times n). 
}

\references{
Park, C., H. Kim, and M. Wang (2020).
Investigation of finite-sample properties of robust location and scale estimators.
\emph{Communications in Statistics - Simulation and Computation}, To appear.\cr
\url{https://doi.org/10.1080/03610918.2019.1699114}

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
