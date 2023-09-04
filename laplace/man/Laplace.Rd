\name{Laplace}
\alias{Laplace}
\alias{laplace}
\alias{dlaplace}
\alias{plaplace}
\alias{qlaplace}
\alias{rlaplace}

\title{The Laplace distribution}

\description{Density, distribution function, quantile function and random generation for the Laplace distribution with
location and scale}

\usage{
dlaplace(x, location = 0, scale = 1, rate=1/scale, log = FALSE)
plaplace(q, location = 0, scale = 1, rate=1/scale, lower.tail = TRUE, log.p = FALSE)
qlaplace(p, location = 0, scale = 1, rate=1/scale, lower.tail = TRUE, log.p = FALSE)
rlaplace(n, location = 0, scale = 1, rate=1/scale)
}

\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{location}{location parameter.}
  \item{scale}{scale parameter.}
  \item{rate}{rate parameter.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.}

}

\details{
The probability density function of the Laplace distribution is given by 
\deqn{f(x) = \frac{1}{2\sigma} \exp\left( -\frac{|x-\mu|}{\sigma} \right)}{f(x)=1/(2s) exp(-|x-m|/s)}
}

\value{\code{dlaplace} gives the density, \code{plaplace} gives the distribution function, \code{qlaplace} gives the quantile function, 
and \code{rlaplace} generates random deviates.}

\author{Chanseok Park}

\examples{
dlaplace(1.5,  location=5, scale=2)
exp( dlaplace(1.5, location=5, scale=2, log=TRUE) )

plaplace(2.5,  location=5, scale=2)
1 - plaplace(2.5, location=5,scale=2, lower.tail = FALSE, log.p = FALSE)
1 - exp( plaplace(2.5, location=5,scale=2, lower.tail = FALSE, log.p = TRUE) )

qlaplace(0.1, location=5, scale=2)
qlaplace(0.9, location=5, scale=2, lower.tail = FALSE, log.p = FALSE)
qlaplace(log(0.1), location=5, scale=2, lower.tail = TRUE,  log.p = TRUE)
qlaplace(log(0.9), location=5, scale=2, lower.tail = FALSE, log.p = TRUE)

rlaplace(n=10, location=5, scale=2)
}

\keyword{distribution}
