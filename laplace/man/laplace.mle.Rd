\name{laplace.mle}
\alias{laplace.mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximum likelihood estimates of the Laplace distribution}
\description{Calculates the maximum likelihood estimates of the Laplace distribution.}
\usage{
laplace.mle(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector of observations.}
}
\details{
\code{laplace.mle} obtains the estimates of the location and scale
parameters using the maximum likelihood estimation method.
The probability density function of the Laplace distribution is given by 
\deqn{f(x) = \frac{1}{2\sigma} \exp\left( -\frac{|x-\mu|}{\sigma} \right)}{f(x)=1/(2s) exp(-|x-m|/s)}
Thus, the maximum likelihood estimate of the location parameter, \eqn{\hat{\mu}}, is given by the median
of a sample. The maximum likelihood estimate of the scale parameter is given by 
\deqn{\hat{\sigma} = \frac{1}{n} \sum_{i=1}^{n}|x_i - \hat{\mu}|.}
}

\value{An object of class \code{"laplace.estimate"}, a list with
two parameter estimates}

\author{Chanseok Park}

\seealso{
\code{\link{laplace.lp}} for the parameter estimation using the intercept and slope estimates from the
Laplace probability plot.
}
\examples{
x = c(6.96, 3.91, -5.01, 14.91, 13.65, -13.53, 3.4, -4.04, 55.29, 13.36)

laplace.mle(x)
}
\keyword{parameter estimate}

