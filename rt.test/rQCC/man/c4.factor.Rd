\encoding{utf-8}
\name{correction.factor}
\alias{c4.factor}

\title{Finite-sample correction factor for the standard deviation}

\description{
Finite-sample correction factor for the standard deviation (\eqn{\sigma}{\sigma})
under the normal distribution.
}

\usage{
c4.factor(n, method=c("sd","mad","shamos"))
}
\arguments{
 \item{n}{sample size (\eqn{n \ge 1}{n \ge 1}).}
 \item{method}{a character string specifying the estimator, must be
one of "sd" (default), "mad" and "shamos."}
}


\details{
The conventional sample standard deviation, median absolute deviation (MAD)
and Shamos estimators are Fisher-consistent under the normal distribution, 
but they are not unbiased with a sample of finite size.

An unbiased estimator of (\eqn{\sigma}{\sigma}) 
with the sample standard deviation (\eqn{s}) is given by \cr
\code{sd(x)/c4.factor(length(x), method="sd")} \cr
where \eqn{n}{n} is the sample size.

When the MAD (\code{stats::\link{mad}}) estimator is used, 
an unbiased estimator of (\eqn{\sigma}{\sigma}) is given by \cr
\code{mad(x)/c4.factor(length(x), method="mad")}.

Using the Shamos estimator (\code{rQCC::shamos}), 
an unbiased estimator of (\eqn{\sigma}{\sigma}) is given by \cr 
\code{shamos(x)/c4.factor(length(x), method="shamos")}.

Note that the formula for the correction factor \eqn{c_4(n)}{c4(n)} is given by 
\deqn{c_4(n) = \sqrt{\frac{2}{n-1}}\cdot\frac{\Gamma(n/2)}{\Gamma((n-1)/2)}.}{c4(n) = \sqrt\{2/(n-1)\}*\Gamma(n/2)/\Gamma((n-1)/2).}

The finite-sample correction factors for the MAD and Shamos estimators 
are obtained for \eqn{n=1,2,\ldots,10}{n = 1, 2, ... ,100} 
using the extensive Monte Carlo simulation with 1E07 replicates. 
For the case of \eqn{n > 100}{n > 100}, they are obtained 
using the method of Hayes (2014).
}

\value{
\code{c4.factor} computes the finite-sample correction factor when
the sample standard deviations, the sample MAD and Shamos estimates are used.
}

\author{Chanseok Park and Min Wang}


\references{
Park, C., H. Kim and M. Wang (2019).
Investigation of finite-sample properties of robust location
and scale estimators.  \emph{Submitted for Publication}.

Shamos, M. I. (1976).  Geometry and statistics: Problems at the interface.
In Traub, J. F., editor,
\emph{Algorithms and Complexity: New Directions and Recent Results},
pages 251--280. Academic Press, New York.

Hayes, K. (2014).
Finite-sample bias-correction factors for the median absolute deviation.
\emph{Communications in Statistics: Simulation and Computation},
\bold{43}, 2205--2212.
}

\seealso{
\code{stats::\link{mad}} for the Fisher-consistent MAD estimator 
for the standard deviation (\eqn{\sigma}{\sigma}) of a normal distribution. 

\code{rQCC::mad.unbiased} for finite-sample unbiased MAD estimator
for the standard deviation (\eqn{\sigma}{\sigma}) of a normal distribution. 

\code{rQCC::shamos} for the Fisher-consistent Shamos estimator 
for the standard deviation (\eqn{\sigma}{\sigma}) of a normal distribution. 

\code{rQCC::shamos.unbiased} for finite-sample unbiased Shamos estimator
for the standard deviation (\eqn{\sigma}{\sigma}) of a normal distribution. 

\code{rQCC::n.times.eBias.of.mad} 
for the values of the empirical biases of 
the MAD estimator under the standard normal distribution. 

\code{rQCC::n.times.eBias.of.shamos} 
for the values of the empirical biases of 
the Shamos estimator under the standard normal distribution. 
}


\examples{
c4.factor(n=10, method="sd")

c4.factor(n=10, method="mad")

c4.factor(n=10, method="shamos")
}

\keyword{correction factor}
