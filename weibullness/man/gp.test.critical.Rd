\name{gp.test.critical}
\alias{gp.test.critical}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Critical Value for the Gumbel Goodness-of-Fit Test}
\description{
Calculates the critical value for the Gumbel goodness-of-fit test}
\usage{
gp.test.critical(alpha, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{the significance level.}
  \item{n}{the sample size.}
}
\details{
This function calculates the critical value for the Gumbel Goodness-of-Fit test
which is  constructed using the sample correlation from 
the associated Gumbel probability plot.
The critical value is then looked up in Gumbel.Plot.Quantiles.
There is \code{\link{print}} method for class \code{"gp.test.critical"}. 
} 
\value{
A list with class \code{"gp.test.critical"} containing the following components:
   \item{sample.size}{sample size (missing observations are deleted).}
   \item{alpha}{significance level.}
   \item{critical.value}{critical value.}
   \item{data.name}{a character string giving the name(s) of the data.}
}
\references{
Kinnison, R. (1989). 
Correlation Coefficient Goodness-of-Fit Test for the Extreme-Value Distribution. 
\emph{The American Statistician}, 
\bold{43}(2), 98-100. 

Vogel, R. M. and C. N. Kroll (1989).
Low-Flow Frequency Analysis Using Probability-Plot Correlation Coefficients.
\emph{Journal of Water Resources Planning and Management},
\bold{115}, 338-357.
}
\author{Chanseok Park}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ks.test}} for performing the Kolmogorov-Smirnov test 
for the goodness of fit test of two samples. \cr 
\code{\link{wp.test}} for performing the Weibullness test. \cr 
\code{\link{shapiro.test}} for  performing the Shapiro-Wilk test for normality.
}
\examples{
# Critical value with alpha (significance level) and n (sample size).
gp.test.critical(alpha=0.01, n=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{critical value}
\keyword{Gumbel probability plot}
