\name{weibull.ic}
\alias{weibull.ic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximum likelihood estimates with Interval Censoring}
\description{Calculates the maximum likelihood estimates with Interval Censoring Using the EM Algorithm.
}
\usage{
weibull.ic(X, start=c(1,1), maxits=10000, eps=1E-5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric matrix (n x 2) of observations.}
  \item{start}{a starting value.}
  \item{maxits}{the maximum number of iterations.}
  \item{eps}{the desired accuracy (convergence tolerance).}
}
\details{
The expectation-maximization(EM) algorithm is used for estimating the parameters with interval-censored data.
}
\value{Calculates the maximum likelihood estimates with interval-censored data}
\references{
Park, C. (2023). 
A Note on Weibull Parameter Estimation with Interval Censoring Using the EM Algorithm.
\emph{Mathematics}, \bold{11}(14), 3156. \cr
\doi{10.3390/math11143156}

Lawless, J. F. (2003).
Statistical Models and Methods for Lifetime Data, 2nd ed.; John Wiley & Sons: New York, NY.
}
\author{Chanseok Park}
%% ~Make other sections like Warning with \section{Warning }{\ldots} ~
\seealso{
\code{\link{weibull.wp}}  for the parameter estimation using the Weibull plot with full observations.
\code{\link{weibull.mle}} for the parameter estimation using the maximum likelihood method with full observations.
}
\examples{
library(weibullness)

attach(Wdata)
weibull.ic(radio.chemotherapy)

# Two-parameter Weibull with full observations
weibull.ic( cbind(bearings,bearings) )

# Two-parameter Weibull with full observations (using weibull.mle)
weibull.mle(bearings, threshold=0)
}

\keyword{Maximum likelihood estimator (MLE)}
\keyword{interval censored}
\keyword{parameter}
\keyword{estimate}
