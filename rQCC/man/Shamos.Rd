\encoding{utf-8}
\name{Shamos}
\alias{shamos}
\alias{shamos.unbiased}
\alias{shamos2.unbiased}

\title{Shamos estimator}

\description{Calculates the conventional Shamos, unbiased Shamos and unbiased squared 
Shamos estimators.
The conventional Shamos is calculated by \code{shamos}
which is Fisher-consistent under the normal distribution.
Note that it is not unbiased with a sample of finite size.
The unbiased Shamos estimator under the normal distribution is
calculated by \code{shamos.unbiased} with a finite-sample unbiasing factor.
The unbiased \emph{squared} Shamos estimator under the normal distribution is
calculated by \code{shamos2.unbiased} with a finite-sample unbiasing factor.
}

\usage{
shamos(x, constant=1.048358, na.rm = FALSE,  IncludeEqual=FALSE)

shamos.unbiased(x, constant=1.048358, na.rm = FALSE,  IncludeEqual=FALSE)

shamos2.unbiased(x, constant=1.048358, na.rm = FALSE, IncludeEqual=FALSE)
}
\arguments{
 \item{x}{a numeric vector of observations.}
 \item{constant}{Correction factor for the Fisher-consistency under the normal distribution}
 \item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
 \item{IncludeEqual}{\code{FALSE} (default) calculates median of 
\eqn{|X_i-X_j|}{|Xi-Xj|} with \eqn{i < j}{i < j}, 
while \code{TRUE} calculates median of \eqn{|X_i-X_j|}{|Xi-Xj|} with \eqn{i \le j}{i \le j}.} 
}

\details{
The Shamos estimator is defined as 
\deqn{\textrm{Shamos}=\code{constant}\times\mathop{\mathrm{median}}_{i < j} \big(|X_i-X_j|\big)}{Shamos = constant * median of |Xi-Xj| over i<j} 
where \eqn{i,j=1,2,\ldots,n}{i, j=1,2,...,n}. 
The default value (\code{constant=1.048358}) ensures the Fisher-consistency under the normal distribution.
Note that 
\eqn{\code{constant}=1/\{\sqrt{2}\,\Phi^{-1}(3/4)\}\approx 1.048358}{constant=1/(\sqrt(2)*\Phi^(-1)(3/4)) = 1.048358 (approximately)}.

The unbiased Shamos is defined as
\deqn{\textrm{Shamos}=\code{constant}\times\mathop{\mathrm{median}}_{i < j} \big(|X_i-X_j|\big)/{c_6(n)}}{Shamos = constant * median of |Xi-Xj| over i<j divided by c6(n)}
for \eqn{i,j=1,2,\ldots,n}{i,j=1,2,...,n}, where 
\eqn{c_6(n)}{c6(n)} is the finite-sample unbiasing factor.
Note that \eqn{c_6(n)}{c6(n)} notation is used in Park et. al (2022), and 
\eqn{c_6(n)}{c6(n)} is calculated using the function \code{c4.factor}\{rQCC\} with \code{estimator="shamos"} option.

The unbiased squared Shamos is defined as the
squared \code{shamos}\{rQCC\} divided by  \eqn{w_6(n)}{w6(n)} where 
\eqn{w_6(n)}{w6(n)} is the finite-sample unbiasing factor.
Note that \eqn{w_6(n)}{w6(n)} notation is used in Park et. al (2022), and 
\eqn{w_6(n)}{w6(n)} is calculated using the function \code{w4.factor}\{rQCC\} 
with \code{estimator="shamos2"} option.
Note that the square of the conventional Shamos estimator is
Fisher-consistent for the variance (\eqn{\sigma^2}{\sigma^2}) under the normal distribution, but
it is not unbiased with a sample of finite size.
}

\value{They return a numeric value.}

\references{
Park, C., H. Kim, and M. Wang (2022).
Investigation of finite-sample properties of robust location and scale estimators.
\emph{Communications in Statistics - Simulation and Computation}, 
\bold{51}, 2619-2645. \cr 
\doi{10.1080/03610918.2019.1699114}

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
\code{mad.unbiased}\{rQCC\} for calculating the \emph{unbiased} sample MAD.

\link{mad}\{stats\} for calculating the Fisher-consistent sample MAD.

\code{c4.factor}\{rQCC\} for finite-sample unbiasing 
factor for the standard deviation (\eqn{\sigma}{\sigma}) under the normal distribution.

\code{w4.factor}\{rQCC\} for finite-sample unbiasing  factor for the squared Shamos estimator
of the variance (\eqn{\sigma^2}{\sigma^2}) under the normal distribution.

\code{finite.breakdown}\{rQCC\} for calculating the finite-sample breakdown point.
}

\examples{
x = c(0:10, 50)

# Fisher-consistent Shamos, but not unbiased with a finite sample. 
shamos(x)

# Unbiased Shamos. 
shamos.unbiased(x)

# Fisher-consistent squared Shamos, but not unbiased with a finite sample. 
shamos(x)^2 

# Unbiased squared Shamos. 
shamos2.unbiased(x)
}

\author{Chanseok Park and Min Wang}
\keyword{robust}
\keyword{univar}
