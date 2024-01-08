\name{ep.test.critical}
\alias{ep.test.critical}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Critical Value for the Exponential Goodness-of-Fit test}
\description{
Calculates the critical value for the Weibullness test}
\usage{
ep.test.critical(alpha, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{the significance level.}
  \item{n}{the sample size.}
}
\details{
This function calculates the critical value for the Weibullness test
which is  constructed using the sample correlation 
from the  associated Weibull plot.
The critical value is then looked up in Weibull.Plot.Quantiles.
There is \code{\link{print}} method for class \code{"ep.test.critical"}. 
} 
\value{
A list with class \code{"ep.test.critical"} containing the following components:
   \item{sample.size}{sample size (missing observations are deleted).}
   \item{alpha}{significance level.}
   \item{critical.value}{critical value.}
   \item{data.name}{a character string giving the name(s) of the data.}
}

\author{Chanseok Park}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ks.test}} for performing the Kolmogorov-Smirnov test 
for the goodness of fit test of two samples. \cr 
\code{\link{shapiro.test}} for performing the Shapiro-Wilk test for normality. \cr
\code{\link{wp.test}} for performing the Weibullness test. 
}
\examples{
# Critical value with alpha (significance level) and n (sample size).
ep.test.critical(alpha=0.0982, n=7)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{critical value}
\keyword{exponential probability plot}
