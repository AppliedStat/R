\name{bs.test}
\alias{bs.test}
\alias{bs.test.pvalue}
\alias{bs.test.critical}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The Birnbaum-Saunders goodness-of-fit test from the probability plot}
\description{
Performs goodness-of-fit test for the Birnbaum-Saunders distribution}
\usage{
bs.test(x, a)
bs.test.pvalue(r, n)
bs.test.critical(alpha, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector of data values. Missing values are allowed,
   but the number of non-missing values must be between 3 and 1000.  } 
  \item{a}{the offset fraction to be used; typically in (0,1). See \code{\link{ppoints}}.}
  \item{r}{the sample correlation coefficient 
           from the Birnbaum-Saunders probability plot; r is in (0,1).}
  \item{n}{the sample size.}
  \item{alpha}{the significance level.}
}
\details{Using the sample correlation coefficeint 
         from the Birnbaum-Saunders probability plot, 
         it performs the goodness fit test for the Birnbaum-Saunders distribution.}
\value{A list with class "htest" containing the following components:
   \item{statistic}{the value of the test statistic 
         (sample correlation from the Birnbaum-Saunders probability plot)}
   \item{p.value}{the p-value for the test.}
   \item{sample.size}{sample size (missing observations are deleted).}
   \item{method}{a character string indicating the goodness fit test for the Birnbaum-Saunders distribution.}
   \item{data.name}{a character string giving the name(s) of the data.}
}
%% \references{ }
\author{Chanseok Park}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{ks.test}} for performing a one- or two-sample Kolmogorov-Smirnov test. \cr 
 \code{\link{shapiro.test}} for  performing the Shapiro-Wilk test of normality. \cr
 \code{wp.test}\{weibullness\} for performing the Weibullness test. 
}
\examples{
# For the goodness of fit test
x = c(1.2, 2.0, 3.3)
bs.test(x)

# p.value with r (sample correlation from the probability plot) and n (sample size) are given
bs.test.pvalue(r=0.6, n=10)

# Critical value with alpha (significane level) and n (sample size).
bs.test.critical(alpha=0.01, n=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
