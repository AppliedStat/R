\name{wp.test.pvalue}
\alias{wp.test.pvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The p-value for the Weibullness test}
\description{
Calculates the p-value for the Weibullness test which is based on
the sample correlation from the Weibull plot.}
\usage{
wp.test.pvalue(r, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{the sample correlation coefficient from the Weibull plot; 
           r is in (0,1).}
  \item{n}{the sample size.}
}
\details{The p-value for the Weibullness test which is based on
the sample correlation from the Weibull plot.
There is \code{\link{print}}  method for class \code{"htest"}. 
}
\value{
A list with class "htest" containing the following components:
   \item{statistic}{the value of the test statistic (sample correlation from the Weibull plot)}
   \item{p.value}{the p-value for the test.}
   \item{method}{a character string indicating the Weibullness test.}
}
\references{
Park, C. (2017).
Weibullness test and parameter estimation of the three-parameter
Weibull model using the sample correlation coefficient. 
\emph{International Journal of Industrial Engineering - Theory,
Applications and Practice},  
\bold{24}(4), 376-391. \cr 
\doi{10.23055/ijietap.2017.24.4.2848}

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
  for the goodness of fit test of two samples. 

\code{\link{shapiro.test}} for  performing the Shapiro-Wilk test for normality.
}
\examples{
library(weibullness)

# p.value with r (sample correlation from the Weibull plot) and n (sample size).
wp.test.pvalue(r=0.6, n=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{p-value}
\keyword{Weibull plot}

%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
