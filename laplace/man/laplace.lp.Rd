\name{laplace.lp}
\alias{laplace.lp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate of location and scale parameters of the Laplace distribution}
\description{Calculates the estimates of the location and scale
parameters fromn the Laplace probability plot.}
\usage{
laplace.lp(x, n, a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{a numeric vector of observations.}
 \item{n}{The number of observations is needed if there is right-censoring.}
 \item{a}{the offset fraction to be used; typically in (0,1). See \code{ppoints}.}
}
\details{ 
\code{laplace.lp} obtains the estimates of the location and scale
parameters using the intercept and slope estimates from the 
Laplace probability plot.
}
\value{An object of class \code{"laplace.estimate"}, a list with
two parameter estimates}

\author{Chanseok Park}
\seealso{
\code{\link{laplace.mle}} for the maximum likelihood estimate. 
}
\examples{
x = c(6.96, 3.91, -5.01, 14.91, 13.65, -13.53, 3.4, -4.04, 55.29, 13.36)

laplace.lp(x)
}
\keyword{Laplace probability plot}
\keyword{parameter estimate}
