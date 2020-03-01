\encoding{utf-8}
\name{relative.efficiency}
\alias{RE}
\title{Relative efficiency (RE)}
\usage{
RE (n, method = c("mean", "median","HL1","HL2","HL3", "sd", "range","mad","shamos"),
       correction = FALSE)
}
\arguments{
 \item{n}{sample size (\eqn{n \ge 1}{n \ge 1}).}
 \item{correction}{A finite-sample bias correction.
\code{TRUE} adjusts a finite-sample bias correction using \code{A3.factor} function.}
 \item{method}{a character string specifying the estimator, must be 
one of \code{"mean"} (default), \code{"median"}, \code{"HL1"}, \code{"HL2"}, \code{"HL3"}, 
\code{"range"}, \code{"mad"}, and \code{"shamos"}.}
}

\description{Calculates the relative efficiency value of 
the mean, median and Hodges-Lehmann (HL1, HL2, HL3) estimators with respect to the sample mean
and that of the standard deviation, range, median absolute deviation (MAD) and 
Shamos estimators with respect to the sample standard deviation.
Note that the relative efficiency value of the mean and the standard deviation
is always one by definition.

For the case of the sample mean, standard deviation and range, 
it is possible to derive their variances in analytic form, 
but, for the other case, it may be impossible.
In this case, 
the variances with \eqn{n=1,2,\ldots,100}{n = 1, 2, ... ,100} are obtained
using the extensive Monte Carlo simulation with 1E07 replicates.
For \eqn{n > 100}{n > 100}, the variances are approximated based on the method of Hayes (2014).

To obtain the relative efficiency value of the unbiased scale estimators, 
use \code{correction=TRUE} option.
Note that the location estimators 
(\code{"mean"}, \code{"median"}, \code{"HL1"}, \code{"HL2"}, \code{"HL3"}) are unbiased. 
}

\details{
\code{RE} calculates the relative efficiency value
of a location estimator (\code{"median"}, \code{"HL1"}, \code{"HL2"}, \code{"HL3"}) 
with respect to the sample mean and
that of a scale estimator (\code{"range"}, \code{"mad"}, \code{"shamos"}) 
with respect to the sample standard deviation.

Note that the relative efficiency (RE) of \eqn{\hat{\theta}_2}{\theta2_hat} with respect to 
\eqn{\hat{\theta}_1}{\theta1_hat} is defined as 
\deqn{\mathrm{RE}(\hat{\theta}_2 | \hat{\theta}_1) 
=\frac{\mathrm{Var}(\hat{\theta}_1)}{\mathrm{Var}(\hat{\theta}_2)}.}{RE(\theta2_hat|\theta1_hat) = Var(\theta1_hat) divided by Var(\theta2_hat).}
}

\value{It returns a numeric value.}

\seealso{
\code{n.times.eVar.of.HL1}\{rQCC\} 
for the empirical variance of the HL1 estimator (times \eqn{n}{n}).\cr
\code{n.times.eVar.of.HL2}\{rQCC\} 
for the empirical variance of the HL2 estimator (times \eqn{n}{n}).\cr
\code{n.times.eVar.of.HL3}\{rQCC\} 
for the empirical variance of the HL3 estimator (times \eqn{n}{n}).\cr
\code{n.times.eVar.of.mad}\{rQCC\} 
for the empirical variance of the MAD estimator (times \eqn{n}{n}).\cr
\code{n.times.eVar.of.median}\{rQCC\} 
for the empirical variance of the median estimator (times \eqn{n}{n}).\cr
\code{n.times.eVar.of.shamos}\{rQCC\} 
for the empirical variance of the Shamos estimator (times \eqn{n}{n}).
}

\references{
Park, C., H. Kim, and M. Wang (2020).
Investigation of finite-sample properties of robust location and scale estimators.
\emph{Communications in Statistics - Simulation and Computation}, To appear.\cr
\url{https://doi.org/10.1080/03610918.2019.1699114}

Hayes, K. (2014).
Finite-sample bias-correction factors for the median absolute deviation.
\emph{Communications in Statistics: Simulation and Computation},
\bold{43}, 2205--2212.
}

\examples{
# The RE of the Hodges-Lehmann (HL2) estimator 
#  with respect to the sample standard deviation under the normal distribution.
RE(n=5, method="HL2")

# RE of the Shamos estimator 
# with respect to the sample standard deviation under the normal distribution.
RE(n=5, method="shamos")

# RE of the unbiased Shamos estimator 
# with respect to the unbiased sample standard deviation under the normal distribution.
RE(n=5, method="shamos", correction=TRUE) 

# RE of the range (maximum minus minimum)
# with respect to the sample standard deviation under the normal distribution.
RE(n=6, method="range")

# RE of the unbiased range ( (maximum - minimum) / d2 )
# with respect to the unibased sample standard deviation under the normal distribution.
RE(n=6, method="range", correction=TRUE)
}

\keyword{variance}
\author{Chanseok Park and Min Wang}
