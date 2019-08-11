\name{wp.test.critical}
\alias{wp.test.critical}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Critical value for the Weibullness test}
\description{
Calculates the critical value for the Weibullness test}
\usage{
wp.test.critical(alpha, n)
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
There is \code{\link{print}} method for class \code{"wp.test.critical"}. 
} 
\value{
A list with class \code{"wp.test.critical"} containing the following components:
   \item{sample.size}{sample size (missing observations are deleted).}
   \item{alpha}{significance level.}
   \item{critical.value}{critical value.}
   \item{data.name}{a character string giving the name(s) of the data.}
}
\references{
Park, C. (2017).
Weibullness test and parameter estimation of the three-parameter
Weibull model using the sample correlation coefficient. 
\emph{International Journal of Industrial Engineering - Theory,
Applications and Practice},  
\bold{24}(4), 376-391. \cr 
\url{http://journals.sfu.ca/ijietap/index.php/ijie/article/view/2848}

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

\code{\link{shapiro.test}} for  performing the Shapiro-Wilk test for normality.
}
\examples{
library(weibullness)

# Critical value with alpha (significance level) and n (sample size).
wp.test.critical(alpha=0.01, n=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{critical value}
\keyword{Weibull plot}
