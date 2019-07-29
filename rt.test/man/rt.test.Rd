\name{Robustified-t-test} 
\alias{rt.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Robustified t-test}
\description{Performs robustified one-sample t-test on a vector of data.}
\usage{
rt.test(x, alternative = c("two.sided", "less", "greater"), 
  mu = 0, test.stat = c("TA", "TB"), conf.level = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of quantiles.}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". 
You can specify just the initial letter.}
  \item{mu}{a number indicating the true value of the mean.}
  \item{test.stat}{a character string specifying the test statistic.}
  \item{conf.level}{confidence level of the interval.}
}
\details{Based on the empirical distributions of the TA statistic (based on median and MAD) 
and the TB statistic (based on Hodges-Lehmann and Shamos),
this function performs one-sample robustified t-test.}
\value{
A list with class "htest" containing the following components:
   \item{statistic}{the value of the test statistic.}
   \item{parameter}{sample size (non-missing observations in the sample).}
   \item{p.value}{the p-value for the test.}
   \item{conf.int}{a confidence interval for the mean appropriate to the specified alternative hypothesis.} 
   \item{estimate}{the specified hypothesized value of the median (TA) or the Hodges-Lehmann (TB).}
   \item{sample.size}{numeric scalar containing the number of non-missing observations in the sample used for the hypothesis test}
   \item{null.value}{the specified hypothesized value of the true mean.}
   \item{alternative}{a character string describing the alternative hypothesis.} 
   \item{method}{a character string indicating which statistic (TA or TB) is used.}
   \item{data.name}{a character string giving the name(s) of the data.}
}


\references{
Park, C. and M. Wang (2018).
Empirical distributions of the robustified \emph{t}-test statistics.
\emph{ArXiv e-prints}, 1807.02215. \cr
\url{https://arxiv.org/abs/1807.02215}

Jeong, R., S. B. Son, H. J. Lee, and H. Kim (2018).
On the robustification of the z-test statistic.
Presented at KIIE Conference, Gyeongju, Korea. April 6, 2018.

Park, C. (2018). 
Note on the robustification of the Student \emph{t}-test statistic using the median 
 and the median absolute deviation.  
\emph{ArXiv e-prints}, 1805.12256. \cr 
\url{https://arxiv.org/abs/1805.12256} 
}
\author{Chanseok Park and Min Wang}
% \note{
% %%  ~~further notes~~
% }
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{t.test}} for performing the Student t-test. \cr 
 \code{\link{prop.test}} for testing the proportion.
}
\examples{
# For robustified t-test (two-sided) using median and MAD (TA).
#    test.stat="TA" (default)
x = rnorm(10) 
rt.test(x)

# For robustified t-test (two-sided) using Hodges-Lehmann and Shamos (TB).
x = rnorm(10)
rt.test(x, test.stat="TB")

# 90% CI (two sides).
x = rnorm(10)
rt.test(x, conf.level=0.9)

# 90% CI (one side).
x = rnorm(10)
rt.test(x, alternative="less", conf.level=0.9)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
