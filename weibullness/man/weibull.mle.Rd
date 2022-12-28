\name{weibull.mle}
\alias{weibull.mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximum likelihood estimates of three-parameter Weibull distribution}
\description{Calculates the maximum likelihood estimates of 
three-parameter Weibull distribution.
}
\usage{
weibull.mle(x, threshold, interval, interval.threshold, extendInt="downX", 
            a, tol = .Machine$double.eps^0.25, maxiter = 1000, trace = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector of observations.}
  \item{threshold}{the threshold parameter value.}
  \item{interval}{a vector containing the end-points of the interval to be estimated
for the shape parameter.}
  \item{interval.threshold}{a vector containing the end-points of the interval 
      to be estimated for the threshold parameter.}
  \item{extendInt}{character string specifying if the interval c(left,right) should be extended or directly produce an error when f() has no differing signs at the endpoints. The default, "downX", keep lowering the the left end of the interval so that f() has different signs. See \code{\link{uniroot}}.}
  \item{a}{the offset fraction to be used; typically in (0,1).}
  \item{tol}{the desired accuracy (convergence tolerance).}
  \item{maxiter}{the maximum number of iterations.}
  \item{trace}{integer number; if positive, tracing information is produced. Higher values giving more details.}
}
\details{
The three-parameter Weibull distribution has the cumulative distribution function
\deqn{F(x) = 1 - \exp\Big[ - \Big( \frac{x-\theta}{\beta} \Big)^{\alpha} \Big],}
where  \eqn{x>\theta}.
The shape (\eqn{\alpha}) and scale (\eqn{\beta}) parameters are estimated 
using the maximum likelihood.
The maximum likelihood estimation is performed using the method by 
Farnum and Booth (1997). 
If the \code{threshold} (\eqn{\theta}) is missing, it is estimated by 
\code{\link{weibull.threshold}}. 
If \code{threshold=0}, then \code{weibull.mle} calculates the maximum likelihood 
estimates of the two-parameter Weibull distribution. 

If \code{interval} is missing, the interval is given by the method in 
Farnum and Booth (1997).

If \code{interval.threshold} is missing, the interval is initally given
by \code{(min(x)-sd(x), min(x))}. If this interval does not include 
the estimate, its lower bound is extended (see also \code{\link{uniroot}}). 

The choice of \code{a} follows \code{\link{ppoints}} function.

Convergence is declared either if \code{f(x) == 0} 
or the change in \code{x} for one step of the algorithm is less than 
\code{tol}  (see also \code{\link{uniroot}}).

If the algorithm does not converge in \code{maxiter} steps, 
a warning is printed and the current approximation is returned
(see also \code{\link{uniroot}}).
}
\value{An object of class \code{"weibull.estimate"}, a list with 
two parameter estimates (if \code{threshold} is given) or three parameter estimates.
}
\references{
Park, C. (2018). 
A Note on the Existence of the Location Parameter Estimate of the Three-Parameter Weibull Model Using the Weibull Plot.
\emph{Mathematical Problems in Engineering}, \bold{2018}, 6056975. \cr
\doi{10.1155/2018/6056975}

Park, C. (2017).
Weibullness test and parameter estimation of the three-parameter
Weibull model using the sample correlation coefficient. 
\emph{International Journal of Industrial Engineering - Theory,
Applications and Practice},  
\bold{24}(4), 376-391. \cr 
\url{http://journals.sfu.ca/ijietap/index.php/ijie/article/view/2848}

Farnum, N. R. and P. Booth (1997). 
Uniqueness of Maximum Likelihood Estimators of the 2-Parameter Weibull Distribution.
\emph{IEEE Transactions on Reliability},
\bold{46}, 523-525.
}
\author{Chanseok Park}
%% ~Make other sections like Warning with \section{Warning }{\ldots} ~
\seealso{
\code{\link{weibull.wp}} for the parameter estimation using the Weibull plot.

\code{\link{weibull.threshold}} for the estimate of the threshold parameter.

\code{\link[MASS:fitdistr]{fitdistr}} for maximum-likelihood fitting of univariate distributions in package \CRANpkg{MASS}. 
}
\examples{
library(weibullness)

# Three-parameter Weibull
data = c(355,725,884,462,1092,190,166,172,188,224,267,298,355,471,
        154,101,76,811,80,249,752,305,301,386,667,212,186,127,
        121,214,242,237,355,210,253,400,401,514,211,285)
weibull.mle(data)

# Two-parameter Weibull
weibull.mle(data, threshold=0)
}

\keyword{Maximum likelihood estimator (MLE)}
\keyword{parameter}
\keyword{estimate}
