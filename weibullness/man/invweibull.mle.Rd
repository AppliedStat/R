\name{invweibull.mle}
\alias{invweibull.mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximum likelihood estimates of the two-parameter inverse Weibull distribution}
\description{Calculates the maximum likelihood estimates of 
the two-parameter Weibull distribution.
}
\usage{
invweibull.mle(x, interval, tol = .Machine$double.eps^0.25, maxiter = 1000, trace = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector of observations.}
  \item{interval}{a vector containing the end-points of the interval to be estimated
for the shape parameter.}
  \item{tol}{the desired accuracy (convergence tolerance).}
  \item{maxiter}{the maximum number of iterations.}
  \item{trace}{integer number; if positive, tracing information is produced. Higher values giving more details.}
}
\details{
The two-parameter inverse Weibull distribution has the cumulative distribution function
\deqn{F(X)=\exp(-(\theta/x)^\beta)}{F(x) = exp[-(s/x)^b]}
where \eqn{x>0}{x>0}, \eqn{\beta>0}{b>0} and \eqn{\theta>0}{s>0}.

The shape (\eqn{\beta}{b}) and scale (\eqn{\theta}{s}) parameters are estimated 
by calling \code{\link{weibull.mle}} on the reciprocally transformed data.
The maximum likelihood estimation on the the reciprocally transformed data is performed using the method by 
Farnum and Booth (1997). 
If \code{interval} is missing, the interval is given by the method in 
Farnum and Booth (1997).

Convergence is declared either if \code{f(x) == 0} 
or the change in \code{x} for one step of the algorithm is less than 
\code{tol}  (see also \code{\link{uniroot}}).

If the algorithm does not converge in \code{maxiter} steps, 
a warning is printed and the current approximation is returned
(see also \code{\link{uniroot}}).
}
\value{An object of class \code{"weibull.estimate"}, a list with 
two parameter estimates.
}
\references{
Farnum, N. R. and P. Booth (1997). 
Uniqueness of Maximum Likelihood Estimators of the 2-Parameter Weibull Distribution.
\emph{IEEE Transactions on Reliability}, \bold{46}, 523-525.
}
\author{Chanseok Park}
%% ~Make other sections like Warning with \section{Warning }{\ldots} ~
\examples{
# Three-parameter Weibull
data = c(355,725,884,462,1092,190,166,172,188,224,267,298,355,471,
        154,101,76,811,80,249,752,305,301,386,667,212,186,127,
        121,214,242,237,355,210,253,400,401,514,211,285)
invweibull.mle(data)
}

\keyword{maximum likelihood estimator}
\keyword{parameter estimate}
