\name{gp.test}
\alias{gp.test}
%- \alias{print.gp.test.critical}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Gumbel Goodness-of-Fit Test from a Gumbel Probability Plot}
\description{
Performs the statistical goodness-of-fit test for the Gumbel distribution
using the sample correlation from the Gumbel probability plot.}
\usage{
gp.test(x, a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector of data values. Missing values are allowed,
   but the number of non-missing values must be between 3 and 1000.  } 
  \item{a}{the offset fraction to be used; typically in (0,1). See ppoints().}
}
\details{The Gumbel goodness-of-fit test is constructed using the sample correlation 
which is calculated using the associated Gumbel probability plot. 
The critical value is then looked up in Gumbel.Plot.Quantiles.
There is \code{\link{print}} method for class \code{"htest"}.
} 
\value{
A list with class "htest" containing the following components:
   \item{statistic}{the value of the test statistic (sample correlation from the Gumbel probability plot)}
   \item{p.value}{the p-value for the test.}
   \item{sample.size}{sample size (missing observations are deleted).}
   \item{method}{a character string indicating the Gumbel goodness-of-fit test.}
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
  for the goodness of fit test of two samples.  \cr 
\code{\link{shapiro.test}} for  performing the Shapiro-Wilk test for normality. \cr 
\code{\link{wp.test}} for performing the Weibullness test.
}
\examples{
# Gumbel goodness-of-fit test. 
x = c(-3.16, -3.07, -2.24, -1.8, -1.48, -0.92, -0.87, -0.41, -0.06, 1.15)
gp.test(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Gumbel probability plot} %% RShowDoc("KEYWORDS")
\keyword{hypothesis test}
