\encoding{utf-8}
\name{relative.efficiency}
\alias{RE}
\title{Relative efficiency (RE)}

\description{\code{RE} calculates the relative efficiency value
of a location estimator (\code{"median"}, \code{"HL1"}, \code{"HL2"}, \code{"HL3"})
with respect to the sample mean and
that of a scale estimator (\code{"range"}, \code{"mad"}, \code{"shamos"})
with respect to the sample standard deviation.}

\usage{
RE(n, estimator=c("mean", "median", "HL1", "HL2", "HL3", "sd", "range", "mad", "shamos"),
      poolType =c("A", "B", "C"), 
      baseEstimator, basePoolType, correction=TRUE, correctionBase)
}
\arguments{
 \item{n}{a vector of sample sizes (\eqn{n_i \ge 1}{ni \ge 1}).}
 \item{estimator}{a character string specifying the estimator, must be 
one of \code{"mean"} (default), \code{"median"}, \code{"HL1"}, \code{"HL2"}, \code{"HL3"}, 
 \code{"sd"}, \code{"range"}, \code{"mad"}, and \code{"shamos"}.}
 \item{poolType}{Type for how to pool estimators, must be 
one of \code{"A"} (default), \code{"B"}, and \code{"C"}.}
 \item{baseEstimator}{a character string specifying the baseline estimator 
on the numerator of the relative efficiency, must be
one of \code{"mean"} (default), \code{"median"}, \code{"HL1"}, \code{"HL2"}, \code{"HL3"},
 \code{"sd"}, \code{"range"}, \code{"mad"}, and \code{"shamos"}.}
 \item{basePoolType}{Type for how to pool baseline estimator, must be 
one of \code{"A"}, \code{"B"}, and \code{"C"}. 
If missing, \code{basePoolType <- poolType} is used.}
 \item{correction}{logical. A finite-sample bias correction for the estimator with a single sample.
\code{TRUE} (default) adjusts a finite-sample bias correction for a scale estimator using \code{c4.factor} function.}
 \item{correctionBase}{logical. 
A finite-sample bias correction for the baseline estimator with a single sample. 
If missing, \code{correctionBase <- correction} is used.}
}

\details{Under the assumption of the normal distribution, 
the function calculates the relative efficiency value of 
the mean, median and Hodges-Lehmann (HL1, HL2, HL3) estimators with respect to the selected baseline
estimator (default is the sample mean)
and that of the standard deviation, range, median absolute deviation (MAD) and 
Shamos estimators with respect to the selected baseline estimator 
(default is the sample standard deviation).

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

If \code{n} is a vector of multiple values (multiple samples), the \code{RE} function calculates 
the relative efficiency value of the pooled estimator. In this case, only unbiased pooled estimator
 and baseline estimator are considered. 
That is, we use \code{correction=TURE} and \code{correctionBase=TURE} for multiple samples. 

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
Park, C., H. Kim, and M. Wang (2022).
Investigation of finite-sample properties of robust location and scale estimators.
\emph{Communications in Statistics - Simulation and Computation}, 
\bold{51}, 2619-2645. \cr 
\doi{10.1080/03610918.2019.1699114}


Hayes, K. (2014).
Finite-sample bias-correction factors for the median absolute deviation.
\emph{Communications in Statistics: Simulation and Computation},
\bold{43}, 2205--2212.
}

\examples{
#################
# Single sample #
#################

# RE of the Hodges-Lehmann (HL2) estimator 
#  with respect to the sample standard deviation under the normal distribution.
RE(n=5, estimator="HL2")

# RE of the unbiased Shamos estimator 
# with respect to the unbiased sample standard deviation under the normal distribution.
RE(n=5, estimator="shamos")

# RE of the original Shamos estimator 
# with respect to the sample standard deviation under the normal distribution.
RE(n=5, estimator="shamos", correction=FALSE)

# RE of the unbiased range ( (maximum - minimum) / d2 )
# with respect to the unibased sample standard deviation under the normal distribution.
RE(n=6, estimator="range")

# RE of the original range (maximum minus minimum)
# with respect to the sample standard deviation under the normal distribution.
RE(n=6, estimator="range", correction=FALSE)


####################
# Multiple samples #
####################
# With multiple samples, only the unbiased pooled estimators are considered.

# RE of the pooled median (pooling type A) with respect to the mean (pooling type A)
RE( n=c(4,5), estimator="median" )

# RE of the pooled median (pooling type A) with respect to the median (pooling type C)
RE( n=c(4,5), estimator="median", baseEstimator="median", basePoolType="C")

# RE of the pooled mad (pooling type A) with respect to the standard deviation (pooling type A)
RE( n=c(4,5), estimator="mad")

# RE of the pooled mad (pooling type A) with respect to the standard deviation (pooling type C)
RE( n=c(4,5), estimator="mad", basePoolType="C")

# RE of the pooled standard deviation (pooling type A) with respect to the sd (pooling type C)
RE( n=c(4,5), estimator="sd", baseEstimator="sd", basePoolType="C" )
}

\author{Chanseok Park}
\keyword{variance}
