\name{BS}
\alias{BS}
\alias{dbs}
\alias{pbs}
\alias{qbs}
\alias{rbs}

\title{The Birnbaum-Saunders distribution}

\description{Density, distribution function, quantile function and random generation for the Birnbaum-Saunders distribution with
alpha (shape) and beta (scale)}

\usage{
dbs(x, alpha = 1, beta = 1, log = FALSE)
pbs(q, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
qbs(p, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
rbs(n, alpha = 1, beta = 1)
}

\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{alpha}{shape parameter.}
  \item{beta}{scale parameter.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.}

}

\details{
The Birnbaum-Saunders distribution was proposed by Birnbaum and Saunders (1969) and its probability density function 
and cumulative distribution function are given by 
\deqn{f(x) = \frac{1}{\sqrt{2\pi}} \exp\left[-\frac{1}{2\alpha^{2}}
\left(\frac{x}{\beta}+\frac{\beta}{x}-2\right) \right]
\frac{x^{-\frac{3}{2}} (x+\beta)}{2\alpha\sqrt{\beta}}}
and
\deqn{F(x) = \Phi \Big[ \frac{1}{\alpha} \Big( \sqrt{\frac{x}{\beta}}-\sqrt{\frac{\beta}{x}} \Big) \Big],}
where \eqn{x>0}, \eqn{\alpha>0}, and \eqn{\beta>0}. 
}

\value{\code{dbs} gives the density, \code{pbs} gives the distribution function, \code{qbs} gives the quantile function, 
and \code{rbs} generates random deviates.}

\references{Birnbaum, Z. W. and Saunders, S. C. (1969). A new family of life distributions. J. Appl. Probab. 6(2): 637-652.}

\author{Chanseok Park}

\examples{
dbs(1.5,  alpha=0.5, beta=1.5)
exp( dbs(1.5, alpha=0.5, beta=1.5, log=TRUE) )

pbs(2.5,  alpha=0.5, beta=1.5)
1 - pbs(2.5, alpha=0.5,beta=1.5, lower.tail = FALSE, log.p = FALSE)
1 - exp( pbs(2.5, alpha=0.5,beta=1.5, lower.tail = FALSE, log.p = TRUE) )

qbs(0.1, alpha=0.5, beta=1.5)
qbs(0.9, alpha=0.5, beta=1.5, lower.tail = FALSE, log.p = FALSE)
qbs(log(0.1), alpha=0.5, beta=1.5, lower.tail = TRUE,  log.p = TRUE)
qbs(log(0.9), alpha=0.5, beta=1.5, lower.tail = FALSE, log.p = TRUE)

rbs(n=10, alpha=0.5, beta=1.5)
}

\keyword{distribution}
