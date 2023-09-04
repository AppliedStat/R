\name{laplace.rm}
\alias{laplace.rm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Robust estimate of shape and scale parameters of the Laplace distribution 
using the repeated median method}
%%  ~~function to do ... ~~
\description{Calculates the estimates of the shape and scale
parameters.}
\usage{
laplace.rm(x, a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{a numeric vector of observations.}
 \item{a}{the offset fraction to be used; typically in (0,1). See \code{ppoints}.}
}
\details{ 
\code{laplace.rm} obtains the robust estimates of the location and scale
parameters using the intercept and slope estimates using the repeated median method from the 
Laplace probability plot.
}
\value{An object of class \code{"laplace.estimate"}, a list with
two parameter estimates}
\references{
Siegel, A. F. (1982). Robust Regression Using Repeated Medians. 
\emph{Biometrika}, \bold{69}, 242-244. 
}
\author{Chanseok Park}
\seealso{
\code{\link{laplace.mle}} for the parameter estimation using the maximum likelihood method. \cr 
\code{\link{laplace.lp}} for the parameter estimation using the Laplace probability plot. \cr 
}
\examples{
x = c(3.00, -0.66, 0.15, -2.04, 0.73, -0.02, 1.36, 1.14, 1.39, -1.87)
laplace.rm(x)
}
\keyword{Laplace probability plot}% use one of  RShowDoc("KEYWORDS")
\keyword{parameter estimate}
\keyword{robust}
