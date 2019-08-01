\name{factor.for.squared.mad}
\alias{w5.for.mad2}
\encoding{utf-8}
\title{Finite-sample correction factor for the squared MAD (median absolute deviation) estimator
under the normal distribution}
\usage{
w5.for.mad2(n)
}
\arguments{
 \item{n}{sample size (\eqn{n \ge 1}{n \ge 1}).}
}

\description{Finite-sample correction factor 
for the squared median absolute deviation (MAD) estimator under the normal distribution
with a sample of size n.
Note that the conventional squared median absolute deviation (MAD) estimator is 
Fisher-consistent under the normal distribution, but 
it is not unbiased with a sample of finite size.

The finite-sample correction factors for
\eqn{n=1,2,\ldots,10}{n = 1, 2, ... ,100} are obtained 
using the extensive Monte Carlo simulation with 1E07 replicates. 
For the case of \eqn{n > 100}{n > 100}, they are obtained 
using the method of Hayes (2014).
}

\value{\code{w5.for.mad2} calculates the finite-sample correction factor 
for the variance (\eqn{\sigma^2}{\sigma^2}) under the normal distribution.}

\references{
Park, C., H. Kim and M. Wang (2019).
Investigation of finite-sample properties of robust location
and scale estimators.  \emph{Submitted for Publication}.

Hayes, K. (2014).
Finite-sample bias-correction factors for the median absolute deviation.
\emph{Communications in Statistics: Simulation and Computation},
\bold{43}, 2205--2212.
}

\seealso{
\code{rQCC::mad2.unbiased} for robust finite-sample unbiased squared 
median absolute deviation (MAD) estimator
for the variance (\eqn{\sigma^2}{\sigma^2}) of a normal distribution. \cr

\code{rQCC::shamos2.unbiased} for robust finite-sample unbiased squared Shamos estimator
for the variance (\eqn{\sigma^2}{\sigma^2}) of a normal distribution. \cr

\code{stats::\link{mad}} for calculating the sample median absolute deviation (MAD). \cr

\code{rQCC::shamos} for calculating the sample Shamos estimate.
}

\examples{
w5.for.mad2(n=10)
}


\author{Chanseok Park and Min Wang}
\keyword{correction factor}
\keyword{variance}

