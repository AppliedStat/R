\name{ep.test}
\alias{ep.test}
%- \alias{print.ep.test.critical}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The Exponential Goodness-of-Fit  Test from the Exponential Probability Plot}
\description{
Performs Goodness-of-fit test for the exponential distribution
using the sample correlation from the exponential probability plot.}
\usage{
ep.test(x, a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector of data values. Missing values are allowed,
   but the number of non-missing values must be between 3 and 1000.  } 
  \item{a}{the offset fraction to be used; typically in (0,1). See ppoints().}
}
\details{The exponential goodness-of-fit test is constructed using the sample correlation 
which is calculated using the associated exponential probability plot. 
The critical value is then looked up in Exponential.Plot.Quantiles.
There is \code{\link{print}} method for class \code{"htest"}.
} 
\value{
A list with class "htest" containing the following components:
   \item{statistic}{the value of the test statistic (sample correlation from the exponential probability plot)}
   \item{p.value}{the p-value for the test.}
   \item{sample.size}{sample size (missing observations are deleted).}
   \item{method}{a character string indicating the exponential goodness-of-fit test.}
   \item{data.name}{a character string giving the name(s) of the data.}
}
\references{
Shapiro, S. S. and M. B. Wilk (1972). 
An Analysis of Variance Test for the Exponential Distribution (Complete Samples).
\emph{Technometrics}, \bold{14}(2), 355-370.
}
\author{Chanseok Park}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ks.test}} for performing the Kolmogorov-Smirnov test for the goodness of fit test of two samples.  \cr 
\code{\link{shapiro.test}} for performing the Shapiro-Wilk test for normality. \cr 
\code{\link{wp.test}} for performing the Weibullness test based on the Weibull probability plot.
}
\examples{
# For Exponential GOF. 
# Dataset from Section 2.5 of Shapiro and Wilk (1972).
x = c(6, 1, -4, 8, -2, 5, 0)
ep.test(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{exponential probability plot} %% RShowDoc("KEYWORDS")
\keyword{hypothesis test}
