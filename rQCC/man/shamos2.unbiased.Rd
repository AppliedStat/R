\name{unbiased.squared.shamos}
\alias{shamos2.unbiased}
\encoding{utf-8}
\title{Unbiased squared Shamos estimate}

\usage{
shamos2.unbiased(x, constant=1.048358, na.rm = FALSE, IncludeEqual=FALSE)
}

\arguments{
 \item{x}{a numeric vector of observations.}
 \item{constant}{Correction factor for the Fisher-consistency under the standard normal distribution}
 \item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
 \item{IncludeEqual}{\code{FALSE} (default) calculates median of
\eqn{|X_i-X_j|}{|Xi-Xj|} with \eqn{i < j}{i < j},
while \code{TRUE} calculates median of \eqn{|X_i-X_j|}{|Xi-Xj|} with \eqn{i \le j}{i \le j}.}
}

\description{Calculates the unbiased squared Shamos
estimate under the normal distribution which is adjusted by 
the Fisher-consistency and finite-sample correction factors.}

\details{The unbiased squared Shamos estimator is defined as the 
squared \code{rQCC::shamos} divided by \eqn{w_6(n)}{w6(n)} 
where the finite-sample correction factor  \eqn{w_6(n)}{w6(n)} 
is calculated by \code{rQCC::w6.for.shamos2}.
The default value (\code{constant=1.048358}) ensures the Fisher-consistency under the standard normal.
Note that the square of the conventional Shamos estimator is
Fisher-consistent for the variance (\eqn{\sigma^2}{\sigma^2}) under the normal distribution, but
it is not unbiased with a sample of finite size.
}

\value{\code{shamos2.unbiased} returns a numeric value.}


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
\code{rQCC::w6.for.shamos2} for finite-sample correction factor for the squared Shamos estimator
for the variance (\eqn{\sigma^2}{\sigma^2}) under the normal distribution. \cr

\code{rQCC::shamos.unbiased} for robust finite-sample unbiased Shamos  estimator
for the standard deviation (\eqn{\sigma}{\sigma}) of a normal distribution. \cr

\code{rQCC::finite.breakdown} for calculating the finite-sample breakdown point.
}

\examples{
x = c(0:10, 50)

# Unbiased squared Shamos
shamos2.unbiased(x)

# Fisher-consistent squared Shamos
shamos(x)^2
}
%% \keyword{print}


\author{Chanseok Park and Min Wang}

\keyword{breakdown}
\keyword{robust}
\keyword{variance}
\keyword{unbiased}

