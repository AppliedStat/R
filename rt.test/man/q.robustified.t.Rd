\name{q.robustified.t}
\alias{q.robustified.t}
\title{Lower quantiles of  TA or TB}
\description{Calculates the quantiles of the robustified t-test statistic (TA or TB).}
\usage{
  q.robustified.t(p, n, test.stat=c("TA", "TB"), lower.tail=TRUE)
}
\arguments{
 \item{p}{vector of probabilities.}
 \item{n}{the sample size}
 \item{test.stat}{a character string specifying the test statistic.}
 \item{lower.tail}{logical; if TRUE (default), probabilities are p=P[X <= x],
            otherwise, p=P[X > x].}
}

\details{Using the empirical distributions of TA and TB statistics, it calculates the quantile.}

\references{
Park, C., M. Wang, and W.-Y. Hwang (2022).
Empirical distributions of the robustified \emph{t}-test statistics.
\emph{Industrial Engineering & Management Systems}, \bold{21}(3), 432-439. \cr 
\doi{10.7232/iems.2022.21.3.432}

Park, C., M. Wang, and W.-Y. Hwang (2020).
A study on robustness of the paired sample tests. 
Empirical distributions of the robustified \emph{t}-test statistics.
\emph{Industrial Engineering & Management Systems}, \bold{19}(2), 386-397. \cr
\doi{10.7232/iems.2020.19.2.386} 
}

\author{Chanseok Park}

\seealso{
 \code{\link{qt}} for obtaining quantile value of Student t-distribution.
}

\examples{
# quantile value of TA (using median and MAD) statistic
q.robustified.t(p=0.01, n=10, test.stat="TA")

# quantile value of TB (using Hodges-Lehmann and Shamos) statistic
q.robustified.t(p=0.01, n=10, test.stat="TB")
}
%% \keyword{print}


