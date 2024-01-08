\name{inverseWeibull}
\alias{inverseWeibull}
\alias{dinvweibull}
\alias{pinvweibull}
\alias{qinvweibull}
\alias{rinvweibull}
\title{The inverse Weibull distribution}
\description{
  Density function, distribution function, quantile function, random generation,
  raw moments and limited moments for the Inverse Weibull distribution
  with parameters \code{shape} and \code{scale}.
}
\usage{
dinvweibull(x, shape, scale = 1, log = FALSE)
pinvweibull(q, shape, scale = 1, lower.tail = TRUE, log.p = FALSE)
qinvweibull(p, shape, scale = 1, lower.tail = TRUE, log.p = FALSE)
rinvweibull(n, shape, scale = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is taken to be the number required.}
  \item{shape, scale}{parameters. Must be positive.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities or densities \eqn{p}{p} are given as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
}
\details{
The probability density function of the inverse Weibull distribution with parameters \code{shape} 
\eqn{=\beta}{=b} and \code{scale} \eqn{= \theta}{=s} is given by 
\deqn{f(x) = \frac{\beta (\theta/x)^\beta e^{-(\theta/x)^\beta}}{x}}{f(x) = b (s/x)^b exp(-(s/x)^b)/x}
where \eqn{x > 0}, \eqn{\beta > 0}{b > 0} and \eqn{\theta > 0}{s > 0}.

The cumulative distribution function is given by 
\deqn{F(X)=\exp(-(\theta/x)^\beta)}{F(x) = exp(-(s/x)^b)}
}
\value{
\code{dinvweibull} gives the density, \code{pinvweibull} gives the distribution function,
\code{qinvweibull} gives the quantile function, and \code{rinvweibull} generates random deviates.
}

\author{
  Chanseok Park
}
\examples{
x = (-1):2
names(x) = letters[1:4]
dinvweibull(x, shape=2) 
exp( dinvweibull(x, shape=2, log=TRUE) )


pinvweibull (1, shape=2)
exp(pinvweibull (1, shape=2,  log=TRUE))

q = c(-1,0,1,2)
qinvweibull ( pinvweibull (q, shape=2), shape=2 )
}
\keyword{distribution}
