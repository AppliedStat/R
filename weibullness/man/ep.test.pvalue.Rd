\name{ep.test.pvalue}
\alias{ep.test.pvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The p-value for the Exponential Goodness-of-Fit Test}
\description{
Calculates the p-value for the exponential goodness-of-fit test which is based on
the sample correlation from the exponential probability plot.}
\usage{
ep.test.pvalue(w, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{w}{the ANOVA test statistic; w is in (0,1).}
  \item{n}{the sample size.}
}
\details{The p-value for the exponential goodness-of-fit test which is based on
the sample correlation from the exponential probability plot.
There is \code{\link{print}}  method for class \code{"htest"}. 
}
\value{
A list with class "htest" containing the following components:
   \item{statistic}{the value of the test statistic (sample correlation from the exponential probability plot)}
   \item{p.value}{the p-value for the test.}
   \item{method}{a character string indicating the exponential probability test.}
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
\code{\link{ks.test}} for performing the Kolmogorov-Smirnov test 
  for the goodness of fit test of two samples.  \cr 
\code{\link{shapiro.test}} for performing the Shapiro-Wilk test for normality. \cr 
\code{\link{wp.test}} for performing the Weibullness test.
}
\examples{
# p.value with w (ANOVA statistic) and n (sample size).
ep.test.pvalue(w=0.35593, n=7)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{p-value}
%% \keyword{exponential probability plot}

%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
