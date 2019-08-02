\name{unbiased.shamos}
\alias{shamos.unbiased}
\encoding{utf-8}
\title{Unbiased Shamos estimate}
\usage{
shamos.unbiased(x, constant=1.048358, na.rm = FALSE,  IncludeEqual=FALSE)
}
\arguments{
 \item{x}{a numeric vector of observations.}
 \item{constant}{Correction factor for the Fisher-consistency under the standard normal distribution}
 \item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
 \item{IncludeEqual}{\code{FALSE} (default) calculates median of 
\eqn{|X_i-X_j|}{|Xi-Xj|} with \eqn{i < j}{i < j}, 
while \code{TRUE} calculates median of \eqn{|X_i-X_j|}{|Xi-Xj|} with \eqn{i \le j}{i \le j}.} 
}

\description{Calculates the unbiased Shamos estimate under the normal distribution
with the adjustment by a 
the Fisher-consistency factor and a finite-sample correction factor.}

\details{
The unbiased Shamos is defined as 
\deqn{\textrm{Shamos}=\code{constant}\times\mathop{\mathrm{median}}_{i < j} \big(|X_i-X_j|\big)/{c_6(n)}}{Shamos = constant * median of (Xi+Xj)/2 over i<j divided by c6(n)} 
where \eqn{i,j=1,2,\ldots,n}{i, j=1,2,...,n} and the unbiasing factor with a sample of finite size,
\eqn{c_6(n)}{c6(n)}, is calculated by \code{rQCC::c4.factor} with \code{method="shamos"} option.

The default value (\code{constant=1.048358}) ensures the Fisher-consistency under the standard normal.
Note that 
\eqn{\code{constant}=1/\{\sqrt{2}\,\Phi^{-1}(3/4)\}\approx 1.048358}{constant=1/(\sqrt(2)*\Phi^(-1)(3/4)) = 1.048358 (approximately)}.
}

\value{\code{hamos.unbiased} returns a numeric value.}


\references{
Park, C., H. Kim, and M. Wang (2019).
Finite-sample properties of robust location and scale estimators. 
\emph{arXiv:1908.00462}.

Shamos, M. I. (1976).  Geometry and statistics: Problems at the interface.
In Traub, J. F., editor,
\emph{Algorithms and Complexity: New Directions and Recent Results},
pages 251--280. Academic Press, New York.

Lèvy-Leduc, C., Boistard, H., Moulines, E., Taqqu, M. S.,
and Reisen, V. A.  (2011).
Large sample behaviour of some well-known robust estimators under
long-range dependence.
\emph{Statistics}, \bold{45}, 59--71.
}


\seealso{
\code{rQCC::shamos} for robust Fisher-consistent estimator
for the standard deviation of a normal distribution. \cr

\code{rQCC::c4.factor} for finite-sample correction factor for the standard deviation
under the normal distribution. \cr

\code{stats::\link{mad}} for calculating the sample median absolute deviation (MAD). \cr

\code{rQCC::finite.breakdown} for calculating the finite-sample breakdown point.
}

\examples{
x = c(0:10, 50)

# Unbiased Shamos
shamos.unbiased(x)

# Fisher-consistent Shamos
shamos(x)
}
%% \keyword{print}


\author{Chanseok Park and Min Wang}

\keyword{breakdown}
\keyword{robust}
\keyword{unbiased}

