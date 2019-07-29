\encoding{utf-8}
\name{control.limits}
\alias{A3.factor}
\title{Control limits factor for X-bar chart}

\usage{
A3.factor(n, method=c("sd","mad","shamos"), sigma.factor=3)
}

\arguments{
 \item{n}{sample size (\eqn{n \ge 1}{n \ge 1}).}
 \item{method}{a character string specifying the estimator, must be
one of "sd" (default), "mad" and "shamos."}
 \item{sigma.factor}{a factor for the standard deviation (\eqn{\sigma}{\sigma}).  
The American Standard uses "\emph{3*sigma}" limits (0.27\% false alarm rate),
while the British Standard uses "\emph{3.09*sigma}" limits (0.20\% false alarm rate).}
}

\description{
The conventional control limits of X-bar chart are given by 
\deqn{\bar{\bar{X}} \pm A_3 \bar{S}.}{"X.bar.bar + A3*S.bar"  and  "X.bar.bar - A3*S.bar."}

\code{A3.factor(n, "mad")} and \code{A3.factor(n, "shamos")} 
are used when the median absolute deviation (MAD) and Shamos estimators are considered.
}

\value{\code{A3.factor} calculates the control limits factor for X-bar chart.}

\examples{
# The American 3*sigma limits with the sample standard deviation.
A3.factor(n=10, method="sd")

# The British 3.09*sigma limits with the Shamos estimator.
A3.factor(n=10, method="shamos", sigma.factor=3.09)
}

\author{Chanseok Park and Min Wang}
\keyword{correction factor}
