\name{weibull.rm}
\alias{weibull.rm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Robust estimate of shape and scale parameters of Weibull using the repeated median method}
%%  ~~function to do ... ~~
\description{Calculates the estimates of the shape and scale
parameters.}
\usage{
weibull.rm(x, a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{a numeric vector of observations.}
 \item{a}{the offset fraction to be used; typically in (0,1). See \code{ppoints}.}
}
\details{ 
\code{weibull.rm} obtains the robust estimates of the shape and scale
parameters using the intercept and slope estimates using the repeated median method from the 
Weibull plot.
}
\value{An object of class \code{"weibull.estimate"}, a list with
two parameter estimates}
\references{
Siegel, A. F. (1982). Robust Regression Using Repeated Medians. 
\emph{Biometrika}, \bold{69}, 242-244. 
}
\author{Chanseok Park}
\seealso{
\code{\link{weibull.mle}} for the parameter estimation using the maximum likelihood method.

\code{\link{weibull.wp}} for the parameter estimation using the Weibull plot.

\code{\link[MASS:fitdistr]{fitdistr}} for maximum-likelihood fitting of univariate distributions in package \CRANpkg{MASS}.
}
\examples{
library(weibullness)

data = c(355,725,884,462,1092,190,166,172,188,224,267,298,355,471,
        154,101,76,811,80,249,752,305,301,386,667,212,186,127,
        121,214,242,237,355,210,253,400,401,514,211,285)
weibull.rm(data)
}
\keyword{Weibull plot}% use one of  RShowDoc("KEYWORDS")
\keyword{parameter}
\keyword{estimate}
\keyword{robust}
