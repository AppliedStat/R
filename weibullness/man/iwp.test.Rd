\name{iwp.test}
\alias{iwp.test}
%- \alias{print.iwp.test.critical}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weibullness Test from a inverse Weibull Plot}
\description{
Performs the statistical test of inverse Weibullness 
(Goodness-of-fit test for the inverse Weibull distribution) 
using the sample correlation from the inverse Weibull plot.}
\usage{
iwp.test(x, a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector of data values. Missing values are allowed,
   but the number of non-missing values must be between 3 and 1000.  } 
  \item{a}{the offset fraction to be used; typically in (0,1). See ppoints().}
}
\details{The inverse Weibullness test is constructed using the sample correlation 
which is calculated using the associated inverse Weibull plot. 
The critical value is then looked up in IW.Plot.Quantiles.
There is \code{\link{print}} method for class \code{"htest"}.
} 
\value{
A list with class "htest" containing the following components:
   \item{statistic}{the value of the test statistic (sample correlation from the inverse Weibull plot)}
   \item{p.value}{the p-value for the test.}
   \item{sample.size}{sample size (missing observations are deleted).}
   \item{method}{a character string indicating the inverse Weibullness test.}
   \item{data.name}{a character string giving the name(s) of the data.}
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
\code{\link{wp.test}} for performing the Weibullness test. \cr 
\code{\link{ks.test}} for performing the Kolmogorov-Smirnov test 
  for the goodness of fit test of two samples.  \cr 
\code{\link{shapiro.test}} for performing the Shapiro-Wilk test for normality.
}
\examples{
# For inverse Weibullness hypothesis test. 
attach(Wdata)
iwp.test(urinary)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{inverse Weibull plot} %% RShowDoc("KEYWORDS")
\keyword{hypothesis test}
