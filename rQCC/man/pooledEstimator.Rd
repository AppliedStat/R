\encoding{utf-8}
\name{pooled.Estimator}
\alias{pooledEstimator}

\title{Pooled Estimator}

\description{This function calculates the pooled estimator
based on the unbiased estimators such as the  mean, median, Hodges-Lehmann (HL1, HL2, HL3),
standard deviation, range,  
median absolute deviation (MAD) and Shamos estimators.}

\usage{
pooledEstimator(x, estimator = c("mean", "median", "HL1", "HL2", "HL3", 
                                 "sd", "range", "mad", "shamos"), 
                   poolType=c("A", "B", "C") )
}
\arguments{
 \item{x}{a numeric list of observations.}
 \item{estimator}{a character string specifying the estimator, must be
one of \code{"mean"} (default), \code{"median"}, \code{"HL1"}, \code{"HL2"}, \code{"HL3"},
 \code{"sd"}, \code{"range"}, \code{"mad"}, and \code{"shamos"}.}
 \item{poolType}{Type for how to pool estimators, must be
one of \code{"A"} (default), \code{"B"}, and \code{"C"}.}
}

\details{
This function calculates the pooled estimator based on 
one of \code{"mean"} (default), \code{"median"}, \code{"HL1"}, \code{"HL2"}, \code{"HL3"}, 
\code{"sd"}, \code{"mad"}, and \code{"shamos"}, which are all unbiased. 
There are three different methods of pooling the estimators, denoted by 
\code{"A"} (default), \code{"B"}, and \code{"C"}.
For more details on how to pool them, refer to vignette.
}

\value{They return a numeric value.}

\references{
Park, C. and M. Wang (2020).
A study on the X-bar and S control charts with unequal sample sizes.
\emph{Mathematics}, \bold{8}(5), 698. \cr
\url{https://doi.org/10.3390/math8050698}

Park, C., H. Kim, and M. Wang (2022).
Investigation of finite-sample properties of robust location and scale estimators.
\emph{Communications in Statistics - Simulation and Computation}, To appear.\cr
\url{https://doi.org/10.1080/03610918.2019.1699114}
}

\examples{
x1 = c(1,2,4)
x2 = c(2,5,9)
x  = list(x1,x2)

# Pooled sample mean (default) by type "A" pooling
pooledEstimator(x) 

# Pooled sample sd by type "B" pooling
pooledEstimator(x, estimator="sd", pool="B")
}

\author{Chanseok Park}
\keyword{univar}
\keyword{unbalanced}
