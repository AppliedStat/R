\name{bs.mme}
\alias{bs.mme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Methof-of-moments estimates of Birnbaum-Saunders distribution}
\description{Calculates the methof-of-moments estimates of Birnbaum-Saunders distribution.
}
\usage{
bs.mme(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector of observations.}
}
\details{
The Birnbaum-Saunders distribution has the probability density function
\deqn{f(x) = \frac{1}{\sqrt{2\pi}} \exp\left[-\frac{1}{2\alpha^{2}}
\left(\frac{x}{\beta}+\frac{\beta}{x}-2\right) \right]
\frac{x^{-\frac{3}{2}} (x+\beta)}{2\alpha\sqrt{\beta}}}
where \eqn{x>0}, \eqn{\alpha>0}, and \eqn{\beta>0}.

The parameters are estimated using the methof-of-moments estimates method.
}
\value{An object of class \code{"bs.estimate"}, a list with parameter estimates.
}

\references{Birnbaum, Z. W. and Saunders, S. C. (1969). Estimation for a Family of Life Distributions with Applications to Fatigue.  J. Appl. Probab. 6(2): 328-347.}

\author{Chanseok Park}
%% ~Make other sections like Warning with \section{Warning }{\ldots} ~
\seealso{
\code{\link{bs.mle}} for the parameter estimation using the maximum likelihood method.
}
\examples{
# Aluminum-Coupons data set from Birnbaum and Saunders (1969).
data = c(0.37, 0.706, 0.716, 0.746, 0.785, 0.797, 0.844, 0.855, 0.858,
0.886, 0.886, 0.93, 0.96, 0.988, 0.99, 1, 1.01, 1.016, 1.018,
1.02, 1.055, 1.085, 1.102, 1.102, 1.108, 1.115, 1.12, 1.134,
1.14, 1.199, 1.2, 1.2, 1.203, 1.222, 1.235, 1.238, 1.252, 1.258,
1.262, 1.269, 1.27, 1.29, 1.293, 1.3, 1.31, 1.313, 1.315, 1.33,
1.355, 1.39, 1.416, 1.419, 1.42, 1.42, 1.45, 1.452, 1.475, 1.478,
1.481, 1.485, 1.502, 1.505, 1.513, 1.522, 1.522, 1.53, 1.54,
1.56, 1.567, 1.578, 1.594, 1.602, 1.604, 1.608, 1.63, 1.642,
1.674, 1.73, 1.75, 1.75, 1.763, 1.768, 1.781, 1.782, 1.792, 1.82,
1.868, 1.881, 1.89, 1.893, 1.895, 1.91, 1.923, 1.94, 1.945, 2.023,
2.1, 2.13, 2.215, 2.268, 2.44)

bs.mme(data)
}

\keyword{methof-of-moments}
\keyword{parameter}
\keyword{estimate}
