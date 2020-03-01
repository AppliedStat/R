\encoding{utf-8}
\name{factor.for.control.limits}
\alias{A3.factor}
\title{Factor for control limits}

\description{
The conventional control limits of Xbar-S type chart are given by
\deqn{\bar{\bar{X}} \pm A_3 \bar{S}.}{"X.bar.bar + A3*S.bar"  and  "X.bar.bar - A3*S.bar."}

\code{A3.factor} calculates the factor for control limits based on
the specified scale estimator.
}

\usage{
A3.factor(n, method=c("sd", "range", "mad","shamos"), sigma.factor=3)
}

\arguments{
 \item{n}{sample size (\eqn{n \ge 1}{n \ge 1}).}
 \item{method}{a character string specifying the estimator, must be
one of \code{"sd"} (default), \code{"range"}, \code{"mad"} and \code{"shamos"}.}
 \item{sigma.factor}{a factor for the standard deviation (\eqn{\sigma}{\sigma}).  
The American Standard uses "\emph{3*sigma}" limits (0.27\% false alarm rate),
while the British Standard uses "\emph{3.09*sigma}" limits (0.20\% false alarm rate).}
}

\details{\code{A3.factor} calculates the factor for control limits. 
For example, \code{A3.factor(n,"sd")} is used for Xbar-S type chart
and \code{A3.factor(n,"range")} is used for  Xbar-R type chart.

Also, \code{A3.factor(n, "mad")} and \code{A3.factor(n, "shamos")} 
are used when the median absolute deviation (MAD) and Shamos estimators are considered.
}

\value{It returns a numeric value.}

\examples{
# The American 3*sigma limits with the sample standard deviation.
A3.factor(n=10, method="sd")  # Xbar-S chart 

# The American 3*sigma limits with the sample range 
# Note: A2 notation is widely used for the sample range.
A2 = A3.factor(n=10, method="range")  # Xbar-R chart
A2 

# The British 3.09*sigma limits with the Shamos estimator.
A3.factor(n=10, method="shamos", sigma.factor=3.09)
}

\author{Chanseok Park and Min Wang}

\references{
ASTM (2018). 
\emph{Manual on Presentation of Data and Control Chart Analysis} (9th edition).
American Society for Testing and Materials, West Conshohocken, PA.
}

\keyword{control chart}
