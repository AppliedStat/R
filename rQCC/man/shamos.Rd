\name{Shamos}
\alias{shamos}
\encoding{utf-8}
\title{Shamos estimate (Fisher consistent)}
\usage{
shamos(x, constant=1.048358, na.rm = FALSE,  IncludeEqual=FALSE)
}
\arguments{
 \item{x}{a numeric vector of observations.}
 \item{constant}{Correction factor for the Fisher-consistency under the standard normal distribution}
 \item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
 \item{IncludeEqual}{\code{FALSE} (default) calculates median of 
\eqn{|X_i-X_j|}{|Xi-Xj|} with \eqn{i < j}{i < j}, 
while \code{TRUE} calculates median of \eqn{|X_i-X_j|}{|Xi-Xj|} with \eqn{i \le j}{i \le j}.} 
}

\description{Calculates the Shamos estimate with the adjustment by a correction
factor for the Fisher-consistency under the standard normal distribution.}

\details{
The Shamos estimator is defined as 
\deqn{\textrm{Shamos}=\code{constant}\times\mathop{\mathrm{median}}_{i < j} \big(|X_i-X_j|\big)}{Shamos = constant * median of (Xi+Xj)/2 over i<j} 
where \eqn{i,j=1,2,\ldots,n}{i, j=1,2,...,n}. 
The default value (\code{constant=1.048358}) ensures the Fisher-consistency under the standard normal.
Note that 
\eqn{\code{constant}=1/\{\sqrt{2}\,\Phi^{-1}(3/4)\}\approx 1.048358}{constant=1/(\sqrt(2)*\Phi^(-1)(3/4)) = 1.048358 (approximately)}.
}

\value{\code{shamos} returns a numeric value.}

\references{
Park, C., H. Kim and M. Wang (2019).
Investigation of finite-sample properties of robust location
and scale estimators.  \emph{Submitted for Publication}.

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
\code{rQCC::shamos.unbiased} for robust finite-sample unbiased estimator
for the standard deviation of a normal distribution.  \cr

\code{stats::\link{mad}} for calculating the sample median absolute deviation (MAD).\cr

\code{rQCC::finite.breakdown} for calculating the finite-sample breakdown point.
}

\examples{
x = c(0:10, 50)
shamos(x)
}
%% \keyword{print}


\author{Chanseok Park and Min Wang}

\keyword{breakdown}
\keyword{robust}

