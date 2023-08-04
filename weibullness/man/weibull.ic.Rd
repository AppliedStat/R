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
\doi{10.1155/2018/6056975}

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
## Data from Table 3.10 of Lawless (2003). 
X1 = c(
 8, 0, 24, 17, 17, 24, 16, 13, 11, 16,
18, 17, 32, 23, 44, 10, 0, 5, 12, 11,
33, 31, 13, 19, 34, 13, 16, 35, 15, 11,
22, 48, 30, 13, 10, 8, 4, 11, 14, 4,
34, 30, 18, 16, 35, 21, 11 ) 
X2 = c(
 12, 22, 31, 27, 23, 30, 24, Inf, 13, 20, 
 25, 26, Inf, Inf, 48, 35, 5, 8, 20, Inf, 
 40, Inf, 39, 32, Inf, Inf, 24, Inf, 22, 17, 
 32, Inf, 34, Inf, 17, 21, 9, Inf, 19, 8, 
 Inf, 36, 24, 60, 39, Inf, 20)
XX = cbind(X1,X2)
weibull.ic(XX)

# Two-parameter Weibull with full observations
data = c(355,725,884,462,1092,190,166,172,188,224,267,298,355,471,
        154,101,76,811,80,249,752,305,301,386,667,212,186,127,
        121,214,242,237,355,210,253,400,401,514,211,285)
weibull.ic( cbind(data,data) )

# Two-parameter Weibull with full observations (using weibull.mle)
weibull.mle(data, threshold=0)
}

\keyword{Maximum likelihood estimator (MLE)}
\keyword{interval censored}
\keyword{parameter}
\keyword{estimate}
