\encoding{utf-8}
\name{unbiasing.factor}
\alias{c4.factor}
\alias{w4.factor}

\title{Finite-sample unbiasing factor}

\description{
Finite-sample unbiasing factor for estimating the standard deviation (\eqn{\sigma}{\sigma})
and the variance (\eqn{\sigma^2}{\sigma^2}) under the normal distribution.
}

\usage{
c4.factor(n, method=c("sd","range", "mad","shamos"))
w4.factor(n, method=c("mad2","shamos2"))
}
\arguments{
 \item{n}{sample size (\eqn{n \ge 1}{n \ge 1}).}
 \item{method}{a character string specifying the estimator, must be
one of \code{"sd"} (default), \code{"range"}, \code{"mad"}, \code{"shamos"} 
for \code{c4.factor}, and one of \code{"mad2"} (default), \code{"shamos2"} for \code{w4.factor}.} 
}


\details{
The conventional sample standard deviation, range, median absolute deviation (MAD)
and Shamos estimators are Fisher-consistent under the normal distribution, 
but they are not unbiased with a sample of finite size.

Using the sample standard deviation,
an unbiased estimator of the standard deviation (\eqn{\sigma}{\sigma}) is calculated by \cr
\code{sd(x)/c4.factor(length(x), method="sd")} \cr

Using the range (maximum minus minimum), 
an unbiased estimator of \eqn{\sigma}{\sigma} is calculated by \cr
\code{diff(range(x))/c4.factor(length(x), method="range")}

Using the median absolute deviation (\link{mad}\{stats\}),
an unbiased estimator of \eqn{\sigma}{\sigma} is calculated by \cr
\code{mad(x)/c4.factor(length(x), method="mad")}

Using the Shamos estimator (\code{shamos}\{rQCC\}), 
an unbiased estimator of \eqn{\sigma}{\sigma} is calculated by \cr 
\code{shamos(x)/c4.factor(length(x), method="shamos")}

Note that the formula for the unbiasing factor \eqn{c_4(n)}{c4(n)} is given by 
\deqn{c_4(n) = \sqrt{\frac{2}{n-1}}\cdot\frac{\Gamma(n/2)}{\Gamma((n-1)/2)}.}{c4(n) = \sqrt\{2/(n-1)\}*\Gamma(n/2)/\Gamma((n-1)/2).}  \cr

The \emph{squared} MAD and \emph{squared} Shamos are Fisher-consistent for the variance
(\eqn{\sigma^2}{\sigma^2}) under the normal distribution,
but they are not unbiased with a sample of finite size.

An unbiased estimator of the variance (\eqn{\sigma^2}{\sigma^2})
is obtained using the finite-sample unbiasing factor (\code{w4.factor}).

Using the squared MAD, an unbiased estimator of \eqn{\sigma^2}{\sigma^2} is calculated by \cr
\code{mad(x)^2/w4.factor(length(x), method="mad2")}

Using the squared Shamos estimator,
an unbiased estimator of \eqn{\sigma^2}{\sigma^2} is calculated by \cr
\code{shamos(x)^2/w4.factor(length(x), method="shamos2")}

The finite-sample unbiasing factors for the median absolute deviation (MAD)  
and Shamos estimators 
are obtained for \eqn{n=1,2,\ldots,100}{n = 1, 2, ... ,100} 
using the extensive Monte Carlo simulation with 1E07 replicates. 
For the case of \eqn{n > 100}{n > 100}, they are obtained 
using the method of Hayes (2014).
}

\value{It returns a numeric value.}

\author{Chanseok Park and Min Wang}

\references{
Park, C., H. Kim, and M. Wang (2020).
Investigation of finite-sample properties of robust location and scale estimators.
\emph{Communications in Statistics - Simulation and Computation}, To appear.\cr
\url{https://doi.org/10.1080/03610918.2019.1699114}

Shamos, M. I. (1976).  Geometry and statistics: Problems at the interface.
In Traub, J. F., editor,
\emph{Algorithms and Complexity: New Directions and Recent Results},
251--280. Academic Press, New York.

Hayes, K. (2014).
Finite-sample bias-correction factors for the median absolute deviation.
\emph{Communications in Statistics: Simulation and Computation},
\bold{43}, 2205--2212.
}

\seealso{
\link{mad}\{stats\} for the Fisher-consistent median absolute deviation (MAD) estimator 
of the standard deviation (\eqn{\sigma}{\sigma}) under the normal distribution.

\code{mad.unbiased}\{rQCC\} for finite-sample unbiased median absolute deviation (MAD) estimator
of the standard deviation (\eqn{\sigma}{\sigma}) under the normal distribution.

\code{shamos}\{rQCC\} for the Fisher-consistent Shamos estimator 
of the standard deviation (\eqn{\sigma}{\sigma}) under the normal distribution.

\code{shamos.unbiased}\{rQCC\} for finite-sample unbiased Shamos estimator
of the standard deviation (\eqn{\sigma}{\sigma}) under the normal distribution.

\code{n.times.eBias.of.mad}\{rQCC\} for the values of the empirical bias of 
the median absolute deviation (MAD) estimator under the standard normal distribution. 

\code{n.times.eBias.of.shamos}\{rQCC\} for the values of the empirical bias of 
the Shamos estimator under the standard normal distribution.   \cr 

\code{mad2.unbiased}\{rQCC\} for finite-sample unbiased squared
MAD estimator of the variance (\eqn{\sigma^2}{\sigma^2}) under the normal distribution.

\code{shamos2.unbiased}\{rQCC\} for finite-sample unbiased squared Shamos estimator
of the variance (\eqn{\sigma^2}{\sigma^2}) under the normal distribution.

\code{n.times.evar.of.mad}\{rQCC\} for the values of the empirical variance of
the median absolute deviation (MAD) estimator under the standard normal distribution.

\code{n.times.evar.of.shamos}\{rQCC\} for the values of the empirical variance of
the Shamos estimator under the standard normal distribution.
}


\examples{
# unbiasing factor for estimating the standard deviation
c4.factor(n=10, method="sd")
c4.factor(n=10, method="mad")
c4.factor(n=10, method="shamos")

# Note: d2 notation is widely used for the bias-correction of the range.
d2 = c4.factor(n=10, method="range") 
d2


# unbiasing factor for estimating the variance
w4.factor(n=10, "mad2")
w4.factor(n=10, "shamos2")
}

\keyword{factor}
\keyword{variance}

