\encoding{utf-8}
\name{empirical.variance}
\alias{evar}
\title{Empirical variances of robust estimators}

\description{This function calculate or estimate the variances of the mean, median, Hodges-Lehmann (HL1, HL2, HL3), standard deviation, range,
median absolute deviation (MAD) and Shamos estimators.}

\usage{
evar (n, estimator=c("mean","median","HL1","HL2","HL3", "sd","range","mad","shamos"),
         poolType=c("A","B","C"), correction=TRUE )
}
\arguments{
 \item{n}{a vector of sample sizes.
          For \code{"HL1"}, \code{"sd"}, \code{"range"} and \code{"shamos"}, \eqn{n \ge 2}{n \ge 2}.
          For others, \eqn{n \ge 2}{n \ge 2}.}
 \item{estimator}{a character string specifying the estimator, must be 
one of \code{"mean"} (default), \code{"median"}, \code{"HL1"}, \code{"HL2"}, \code{"HL3"}, \code{"sd"}, \code{"range"}, \code{"mad"}, 
and \code{"shamos"}.}
 \item{poolType}{Type for how to pool estimators, must be
one of \code{"A"} (default), \code{"B"}, and \code{"C"}.}
 \item{correction}{logical. A finite-sample bias correction for the estimator with a single sample.
      \code{TRUE} (default) adjusts a finite-sample bias correction for a scale estimator
      using \code{c4.factor} function.}
}


\details{
This function calculates or estimates 
the variance of a specific estimator 
when a random sample is from the standard normal distribution.

For the mean, standard deviation (sd) and range, their exact variances
are calculated, but the others are empirically estimated through 
the extensive Monte Carlo simulation with 1E07 replicates
for  \eqn{n=1,2,\ldots,100}{n = 1, 2, ... ,100}. 
For the case of \eqn{n > 100}{n > 100}, 
the empirical variances are obtained using the method of Hayes (2014).
}

\value{It returns a numeric value.}

\seealso{
\code{RE}\{rQCC\} for the relative efficiency. \cr 
\code{n.times.eVar.of.HL1}\{rQCC\} for the empirical variance of the HL1 estimator (times n). \cr
\code{n.times.eVar.of.HL2}\{rQCC\} for the empirical variance of the HL2 estimator (times n). \cr
\code{n.times.eVar.of.HL3}\{rQCC\} for the empirical variance of the HL3 estimator (times n). \cr
\code{n.times.eVar.of.mad}\{rQCC\} for the empirical variance of the MAD estimator (times n). \cr
\code{n.times.eVar.of.median}\{rQCC\} for the empirical variance of the median estimator (times n). \cr
\code{n.times.eVar.of.shamos}\{rQCC\} for the empirical variance of the Shamos estimator (times n). 
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
# Empirical variance of the Hodges-Lehmann estimator (HL2) under the standard normal distribution.
evar (n=10, estimator="HL2")

# Multiple samples
evar (n=c(4,5), estimator="mad", poolType="C")
}

\author{Chanseok Park and Min Wang}
\keyword{variance}
\keyword{unbalanced}
