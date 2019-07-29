\name{unbiased.mad}
\alias{mad.unbiased}
\encoding{utf-8}
\title{Unbiased MAD estimate}

\usage{
mad.unbiased(x, center = median(x), constant=1.4826, na.rm = FALSE)
}

\arguments{
 \item{x}{a numeric vector of observations.}
 \item{center}{Optionally, the centre: defaults to the median.}
 \item{constant}{Correction factor for the Fisher-consistency under the standard normal distribution}
 \item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
}

\description{Calculates the unbiased median absolute deviation(MAD) estimate under the normal distribution
which is adjusted by 
the Fisher-consistency and finite-sample correction factors.}

\details{
The unbiased MAD is defined as the \code{stats::\link{mad}} divided by \eqn{c_5(n)}{c5(n)},
where the bias-correction factor with a sample of finite size,
\eqn{c_5(n)}{c5(n)}, is calculated by \code{rQCC::c4.factor} with \code{method="mad"} option.
The default value (\code{constant=1.4826}) ensures the Fisher-consistency under the standard normal.
Note that the original MAD was proposed by Hampel (1974).
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
\code{rQCC::c4.factor} for finite-sample correction factor for the standard deviation
under the normal distribution. \cr

\code{rQCC::shamos} for robust Fisher-consistent estimator
for the standard deviation of a normal distribution. \cr

\code{rQCC::shamos.unbiased} for robust finite-sample unbiased estimator
for the standard deviation of a normal distribution. \cr

\code{stats::\link{mad}} for calculating the sample median absolute deviation (MAD). \cr

\code{rQCC::finite.breakdown} for calculating the finite-sample breakdown point.
}

\examples{
x = c(0:10, 50)

# Unbiased MAD
mad.unbiased(x)

# Fisher-consistent MAD
mad(x)
}
%% \keyword{print}


\author{Chanseok Park and Min Wang}

\keyword{breakdown}
\keyword{robust}
\keyword{unbiased}

