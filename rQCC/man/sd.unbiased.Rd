\encoding{utf-8}
\name{sd.unbiased}
\alias{sd.unbiased}
\title{Unbiased standard deviation}

\description{This function computes the unbiased standard deviation of the values in \code{x}.
  If \code{na.rm} is \code{TRUE}, then missing values are removed before computation proceeds.
}

\usage{
sd.unbiased(x, na.rm = FALSE)
}

\arguments{
  \item{x}{a numeric vector or an \R object which is coercible to one
    by \code{as.double(x)}.}
  \item{na.rm}{logical.  If \code{TRUE}, then missing values are removed before computation proceeds.}
}

\details{
  \code{\link{sd}} is not unbiased. 
  This function computes \code{\link{sd}} divided by \eqn{c_4(n)}{c4(n)}, where \eqn{c_4(n)}{c4(n)}
   is obtained by \code{c4.factor(n, estimator="sd")}\{rQCC\}.
}

\seealso{
  \code{\link{sd}} for the square root of \code{\link{var}}, but it is biased. \cr
   Refer to \code{mad.unbiased}\{rQCC\} for finite-sample unbiased median absolute deviation 
    (MAD) estimator, the most robust alternative. \cr
}

\examples{
sd.unbiased(1:2)
}

\author{Chanseok Park}
\keyword{univar}
