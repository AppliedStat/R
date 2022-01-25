\encoding{utf-8}
\name{MAD}
\alias{mad.unbiased}
\alias{mad2.unbiased}
\title{Median absolute deviation (MAD)}

\description{Calculates the unbiased median absolute deviation (MAD)
estimator and the unbiased squared MAD under the normal distribution which are adjusted by 
the Fisher-consistency and finite-sample unbiasing factors.}

\usage{
mad.unbiased (x, center = median(x), constant=1.4826, na.rm = FALSE)

mad2.unbiased(x, center = median(x), constant=1.4826, na.rm = FALSE)
}

\arguments{
 \item{x}{a numeric vector of observations.}
 \item{center}{pptionally, the center: defaults to the median.}
 \item{constant}{correction factor for the Fisher-consistency under the normal distribution}
 \item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
}

\details{
The unbiased MAD (\code{mad.unbiased})
is defined as the \link{mad}\{stats\} divided by \eqn{c_5(n)}{c5(n)},
where  \eqn{c_5(n)}{c5(n)} is the finite-sample unbiasing factor.
Note that \eqn{c_5(n)}{c5(n)} notation is used in Park et. al (2022),
and \eqn{c_5(n)}{c5(n)} is calculated 
using the function \code{c4.factor}\{rQCC\} with \code{estimator="mad"} option.
The default value (\code{constant=1.4826}) ensures the Fisher-consistency 
under the normal distribution. 
Note that the original MAD was proposed by Hampel (1974).

The unbiased squared MAD (\code{mad2.unbiased}) is defined as the
squared \link{mad}\{stats\} divided by \eqn{w_5(n)}{w5(n)} where 
 \eqn{w_5(n)}{w5(n)} is the finite-sample unbiasing  factor. 
Note that \eqn{w_5(n)}{w5(n)} notation is used in Park et. al (2022),
and \eqn{w_5(n)}{w5(n)} is calculated 
using the function \code{w4.factor}\{rQCC\} with \code{estimator="mad2"} option.
The default value (\code{constant=1.4826}) ensures the Fisher-consistency
 under the normal distribution.
Note that the square of the conventional MAD is
Fisher-consistent for the variance (\eqn{\sigma^2}{\sigma^2}) under the normal distribution, but
it is not unbiased with a sample of finite size.
}

\value{They return a numeric value.}

\references{
Park, C., H. Kim, and M. Wang (2022).
Investigation of finite-sample properties of robust location and scale estimators.
\emph{Communications in Statistics - Simulation and Computation}, To appear.\cr
\doi{10.1080/03610918.2019.1699114}

Hampel, F. R. (1974).
The influence curve and its role in robust estimation.
\emph{Journal of the American Statistical Association}, \bold{69}, 383--393.
}

\seealso{
\code{c4.factor}\{rQCC\} for finite-sample unbiasing  factor for the standard deviation
under the normal distribution. 

\code{w4.factor}\{rQCC\} for finite-sample unbiasing  factor
for the variance under the normal distribution.

\code{shamos}\{rQCC\} for robust Fisher-consistent estimator
of the standard deviation under the normal distribution. 

\code{shamos.unbiased}\{rQCC\} for robust finite-sample unbiased estimator
of the standard deviation under the normal distribution. 

\link{mad}\{stats\} for calculating the sample MAD. 

\code{finite.breakdown}\{rQCC\} for calculating the finite-sample breakdown point.
}

\examples{
x = c(0:10, 50)

# Fisher-consistent MAD, but not unbiased with a finite sample.
mad(x)

# Unbiased MAD.
mad.unbiased(x)

# Fisher-consistent squared MAD, but not unbiased.
mad(x)^2

# Unbiased squared MAD.
mad2.unbiased(x)
}

\author{Chanseok Park and Min Wang}
\keyword{robust}
\keyword{univar}
