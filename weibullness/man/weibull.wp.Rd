\name{weibull.wp}
\alias{weibull.wp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate of shape and scale parameters of Weibull using the Weibull plot
%%  ~~function to do ... ~~
}
\description{Calculates the estimates of the shape and scale
parameters.}
\usage{
weibull.wp(x, n, a=0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{a numeric vector of observations.}
 \item{n}{The number of observations is needed if there is right-censoring.}
 \item{a}{the offset fraction to be used; typically in (0,1).}
}
\details{ 
\code{weibull.wp} obtains the estimates of the shape and scale
parameters using the intercept and slope estimates from the 
Weibull plot.
}
\value{An object of class \code{"weibull.estimate"}, a list with
two parameter estimates}
\references{
Park, C. (2018).
A Note on the Existence of the Location Parameter Estimate of the Three-Parameter Weibull Model Using the Weibull Plot.
\emph{Mathematical Problems in Engineering},
\bold{2018}, 1-6. \cr
\url{https://doi.org/10.1155/2018/6056975}

Park, C. (2017).
Weibullness test and parameter estimation of the three-parameter
Weibull model using the sample correlation coefficient.
\emph{International Journal of Industrial Engineering - Theory,
Applications and Practice},
\bold{24}(4), 376-391. \cr
\url{http://journals.sfu.ca/ijietap/index.php/ijie/article/view/2848}
}
\author{Chanseok Park}
\seealso{
\code{\link{weibull.mle}} for the parameter estimation 
using the maximum likelihood method.

\code{\link[MASS:fitdistr]{fitdistr}} for maximum-likelihood fitting of univariate distributions in package \CRANpkg{MASS}.
}
\examples{
library(weibullness)

data = c(355,725,884,462,1092,190,166,172,188,224,267,298,355,471,
        154,101,76,811,80,249,752,305,301,386,667,212,186,127,
        121,214,242,237,355,210,253,400,401,514,211,285)
weibull.wp(data)
}
\keyword{Weibull plot}% use one of  RShowDoc("KEYWORDS")
\keyword{parameter}
\keyword{estimate}
