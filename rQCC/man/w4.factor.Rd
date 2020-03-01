\encoding{utf-8}
\name{factor.for.variance}
\alias{w4.factor}

\title{Finite-sample unbiasing factor for the variance}

\description{
Finite-sample unbiasing factor of the variance estimators
for the variance (\eqn{\sigma^2}{\sigma^2}) under the normal distribution.
}

\usage{
w4.factor(n, method=c("mad2","shamos2"))
}
\arguments{
 \item{n}{sample size (\eqn{n \ge 1}{n \ge 1}).}
 \item{method}{a character string specifying the estimator, must be
one of \code{"mad2"} (squared MAD) and "\code{shamos2}" (squared Shamos).}

}

\details{
The squared MAD and squared Shamos are Fisher-consistent for the variance 
(\eqn{\sigma^2}{\sigma^2}) under the normal distribution,
but they are not unbiased with a sample of finite size.

An unbiased estimator of the variance (\eqn{\sigma^2}{\sigma^2}) 
is obtained using the finite-sample unbiasing factor (\code{w4.factor}).

Using the MAD, an unbiased estimator of \eqn{\sigma^2}{\sigma^2} is calculated by \cr
\code{mad(x)^2/w4.factor(length(x), method="mad2")} 

Using the Shamos estimator, 
an unbiased estimator of \eqn{\sigma^2}{\sigma^2} is calculated by \cr
\code{shamos(x)^2/w4.factor(length(x), method="shamos2")}

The finite-sample unbiasing  factors for
\eqn{n=1,2,\ldots,100}{n = 1, 2, ... ,100} are obtained 
using the extensive Monte Carlo simulation with 1E07 replicates. 
For the case of \eqn{n > 100}{n > 100}, they are obtained 
using the method of Hayes (2014).
}

\value{It returns a numeric value.}

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

\seealso{
\code{mad2.unbiased}\{rQCC\} for finite-sample unbiased squared 
MAD estimator of the variance (\eqn{\sigma^2}{\sigma^2}) under the normal distribution. 

\code{shamos2.unbiased}\{rQCC\} for finite-sample unbiased squared Shamos estimator
of the variance (\eqn{\sigma^2}{\sigma^2}) under the normal distribution. 
}

\examples{
w4.factor(n=10, "mad2")

w4.factor(n=10, "shamos2")
}


\author{Chanseok Park and Min Wang}
\keyword{unbiasing factor}
\keyword{variance}

