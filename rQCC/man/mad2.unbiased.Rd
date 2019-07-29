\name{unbiased.squared.mad}
\alias{mad2.unbiased}
\encoding{utf-8}
\title{Unbiased squared MAD estimate}

\usage{
mad2.unbiased(x, center = median(x), constant=1.4826, na.rm = FALSE)
}

\arguments{
 \item{x}{a numeric vector of observations.}
 \item{center}{Optionally, the centre: defaults to the median.}
 \item{constant}{Correction factor for the Fisher-consistency under the standard normal distribution}
 \item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
}

\description{Calculates the unbiased squared median absolute deviation(MAD) 
estimate under the normal distribution which is adjusted by 
the Fisher-consistency and finite-sample correction factors.}

\details{The unbiased squared MAD is defined as the 
squared \code{stats::\link{mad}} divided by \eqn{w_5(n)}{w5(n)} 
where the finite-sample correction factor  \eqn{w_5(n)}{w5(n)} 
is calculated by \code{rQCC::w5.for.mad2}.
The default value (\code{constant=1.4826}) ensures the Fisher-consistency under the standard normal.
Note that the square of the conventional MAD estimator is
Fisher-consistent for the variance (\eqn{\sigma^2}{\sigma^2}) under the normal distribution, but
it is not unbiased with a sample of finite size.
}

\references{
Park, C., H. Kim and M. Wang (2019).
Investigation of finite-sample properties of robust location
and scale estimators.  \emph{Submitted for Publication}.

Hampel, F. R. (1974).
The influence curve and its role in robust estimation.
\emph{Journal of the American Statistical Association}, \bold{69}, 383--393.
}


\seealso{
\code{rQCC::w5.for.mad2} for finite-sample correction factor for the squared MAD estimator
for the variance (\eqn{\sigma^2}{\sigma^2}) under the normal distribution. \cr

\code{rQCC::mad.unbiased} for robust finite-sample unbiased MAD  estimator
for the standard deviation  (\eqn{\sigma}{\sigma}) of a normal distribution. \cr

\code{rQCC::finite.breakdown} for calculating the finite-sample breakdown point.
}

\examples{
x = c(0:10, 50)

# Unbiased squared MAD
mad2.unbiased(x)

# Fisher-consistent squared MAD
mad(x)^2
}
%% \keyword{print}


\author{Chanseok Park and Min Wang}

\keyword{breakdown}
\keyword{robust}
\keyword{variance}
\keyword{unbiased}

