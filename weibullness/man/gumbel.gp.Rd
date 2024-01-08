\name{gumbel.gp}
\alias{gumbel.gp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate of location and scale parameters of the Gumbel distribution}
\description{Calculates the estimates of the location and scale
parameters fromn the Gumbel probability plot.}
\usage{
gumbel.gp(x, n, a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{a numeric vector of observations.}
 \item{n}{The number of observations is needed if there is right-censoring.}
 \item{a}{the offset fraction to be used; typically in (0,1). See \code{ppoints}.}
}
\details{ 
\code{gumbel.gp} obtains the estimates of the location and scale
parameters using the intercept and slope estimates from the 
Gumbel probability plot.
}
\value{An object of class \code{"gumbel.estimate"}, a list with
two parameter estimates}

\author{Chanseok Park}
\seealso{
\code{\link{weibull.wp}} for the parameter estimation using the Weibull plot.  \cr 
}
\examples{
x = c(-2.73, 11.69, 34.85, 7.97, -0.86, 17.46, 9.69, -7.14, 14.86, 16.85)
gumbel.gp(x)
}
\keyword{Gumbel probability plot}
\keyword{parameter estimate}
