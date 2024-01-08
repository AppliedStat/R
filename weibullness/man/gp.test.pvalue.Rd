\name{gp.test.pvalue}
\alias{gp.test.pvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The p-value for the Gumbel goodness-of-Fit Test}
\description{
Calculates the p-value for the Gumbel goodness-of-fit test which is based on
the sample correlation from the Gumbel probability plot.}
\usage{
gp.test.pvalue(r, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{the sample correlation coefficient from the Gumbel probability plot; 
           r is in (0,1).}
  \item{n}{the sample size.}
}
\details{The p-value for the Gumbel goodness-of-fit test which is based on
the sample correlation from the Gumbel probability plot.
There is \code{\link{print}}  method for class \code{"htest"}. 
}
\value{
A list with class "htest" containing the following components:
   \item{statistic}{the value of the test statistic (sample correlation from the Gumbel probability plot)}
   \item{p.value}{the p-value for the test.}
   \item{method}{a character string indicating the Gumbel goodness-of-fit test.}
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
\code{\link{wp.test}} for performing the Weibullness test.
\code{\link{shapiro.test}} for performing the Shapiro-Wilk test for normality.
}
\examples{
# p.value with r (sample correlation from the Gumbel probability plot) and n (sample size).
gp.test.pvalue(r=0.98504, n=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{p-value}
\keyword{Gumbel probability plot}

%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
