\name{wp.test}
\alias{wp.test}
\alias{wp.test.pvalue}
\alias{wp.test.critical}
\alias{Weibullness}
\alias{weibullness}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The Weibullness Test from the Weibull Plot}
\description{
Performs the statistical test of Weibullness 
(Goodness-of-fit test for the Weibull distribution) 
using the sample correlation from the Weibull plot.}
\usage{
wp.test(x, a)
wp.test.pvalue(r, n)
wp.test.critical(alpha, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector of data values. Missing values are allowed,
   but the number of non-missing values must be between 3 and 1000.  } 
  \item{a}{the offset fraction to be used; typically in (0,1). See ppoints().}
  \item{r}{the sample correlation coefficient from the Weibull plot; 
           r is in (0,1).}
  \item{n}{the sample size.}
  \item{alpha}{the significance level.}
}
\details{The Weibullness test is constructed using the sample correlation 
which is calculated using the associated Weibull plot. 
The critical value is then looked up in Weibull.Plot.Quantiles.} 
\value{
A list with class "htest" containing the following components:
   \item{statistic}{the value of the test statistic (sample correlation from the Weibull plot)}
   \item{p.value}{the p-value for the test.}
   \item{sample.size}{sample size (missing observations are deleted).}
   \item{method}{a character string indicating the Weibullness test.}
   \item{data.name}{a character string giving the name(s) of the data.}
}
\references{
C. Park (2017).
Weibullness test and parameter estimation of the three-parameter
Weibull model using the sample correlation coefficient. 
\emph{International Journal of Industrial Engineering - Theory,
Applications and Practice},  
\bold{24}(4), 376-391. \cr 
\url{http://journals.sfu.ca/ijietap/index.php/ijie/article/view/2848}


Vogel, R. M. and C. N. Kroll (1989). 
Low-FLow Frequency Analysis Using Probability-Plot Correlation Coefficients.
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
 \code{\link{shapiro.test}} for  performing the Shapiro-Wilk test for normality.
}
\examples{
library(weibullness)
# For Weibullness hypothesis test. 
x = rweibull(10, shape=1)
wp.test(x)

# p.value with r (sample correlation from the Weibull plot) and n (sample size).
wp.test.pvalue(r=0.6, n=10)

# Critical value with alpha (significance level) and n (sample size).
wp.test.critical(alpha=0.01, n=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
