\name{weibull.threshold}
\alias{weibull.threshold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate of threshold parameter of three-parameter Weibull distribution} 
\description{Calculates the estimate of the threshold parameter.}
\usage{
weibull.threshold(x, a, interval.threshold, extendInt="downX")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector of observations.}
  \item{a}{the offset fraction to be used; typically in (0,1).}
  \item{interval.threshold}{a vector containing the end-points of the interval to be estimated for the threshold parameter.}
  \item{extendInt}{character string specifying if the interval c(left,right) should be extended or directly produce an error when f() has no differing signs at the endpoints. The default, "downX", keep lowering the the left end of the interval so that f() has different signs. See \code{\link{uniroot}}.}
}
\details{
The three-parameter Weibull distribution has the cumulative distribution function
\deqn{F(x) = 1 -\exp\Big[-\Big( \frac{x-\theta}{\beta}\Big)^{\alpha}\Big],}
where \eqn{x>\theta}.
The threshold parameter (\eqn{\theta}) is estimated 
by maximizing the correlation function from the Weibull plot.

The choice of \code{a} follows \code{\link{ppoints}} function.

If \code{interval.threshold} is missing, the interval is initially given
by \code{(min(x)-sd(x), min(x))}. If this interval does not include
the estimate, its lower bound is extended (see also \code{\link{uniroot}}).
}
\value{\code{weibull.threshold} returns a numeric value.}

%% ~put references to the literature/web site here ~
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
\doi{10.23055/ijietap.2017.24.4.2848}
}
\author{Chanseok Park}

\seealso{
\code{\link{weibull.mle}} for the maximum likelihood estimate.

\code{\link{weibull.wp}} for the parameter estimation using the Weibull plot.
} 
\examples{
# Data
data = c(355,725,884,462,1092,190,166,172,188,224,267,298,355,471,
        154,101,76,811,80,249,752,305,301,386,667,212,186,127,
        121,214,242,237,355,210,253,400,401,514,211,285)
weibull.threshold(data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Weibull plot}
\keyword{parameter estimate}
