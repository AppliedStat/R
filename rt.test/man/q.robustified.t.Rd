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
Park, C. and M. Wang (2018).
Empirical distributions of the robustified \emph{t}-test statistics.
\emph{ArXiv e-prints}, 1807.02215. \cr
\url{https://arxiv.org/abs/1807.02215}
}

\author{Chanseok Park and Min Wang}

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


