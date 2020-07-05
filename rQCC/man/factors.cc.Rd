\encoding{utf-8}
\name{factors.for.control.chart}
\alias{factors.cc}

\title{Factors for constructing control charts}

\description{
Factors for constructing control charts.
}

\usage{
factors.cc(n, factor=c( "A","A1","A2","A3", "B1","B2","B3","B4","B5","B6", "c2","c4",
                        "d2","d3","D1","D2","D3","D4", "E1", "E2", "E3"), sigma.factor=3)
}

\arguments{
 \item{n}{sample size (\eqn{n \ge 1}{n \ge 2}).}
 \item{factor}{a character string specifying the foctor. \cr
       \code{"c2"}, \code{"c4"}, \code{"d2"}, \code{"d3"} for control chart lines;
       \code{"A"}, \code{"A1"}, \code{"A2"}, \code{"A3"} 
                 for averages for computing control limits; 
       \code{"B1"}, \code{"B2"}, \code{"B3"}, \code{"B4"}, \code{"B5"}, \code{"B6"}
                 for standard devations; 
       \code{"D1"}, \code{"D2"}, \code{"D3"}, \code{"D4"} for ranges;
       and 
       \code{"E1"}, \code{"E2"}, \code{"E3"} for individuals.}
 \item{sigma.factor}{a factor for the standard deviation (\eqn{\sigma}{\sigma}). \cr 
For example, the American Standard uses "\emph{3*sigma}" limits (0.27\% false alarm rate),
while the British Standard uses "\emph{3.09*sigma}" limits (0.20\% false alarm rate).}
}


\details{
The values of the factors are used for constructing 
various control charts.

For example, the conventional \eqn{\bar{X}}{X-bar} chart 
with the sample standard deviation is given by
\deqn{\bar{\bar{X}} \pm A_3 \bar{S}.}{"X.bar.bar + A3*S.bar" and "X.bar.bar - A3*S.bar."}

For more details, refer to \code{vignette("factors.cc", package="rQCC")}.
}

\value{It returns a numeric value.}

\references{
ASTM (2010).
\emph{Manual on Presentation of Data and Control Chart Analysis} (STP 15-D),
8th edition. 
American Society for Testing and Materials, West Conshohocken, PA.

ASTM (1951).
\emph{Manual on Quality Control of Materials} (STP 15-C), 
American Society for Testing and Materials, Philadelphia, PA.
}

\seealso{
\code{c4.factor}\{rQCC\} for \eqn{c_4}{c4} factor
for the finite-sample unbiasing factor to estimate 
the standard deviation (\eqn{\sigma}{\sigma}) under the normal distribution
using various estimators such as the sample standard deviation,
the sample range, the median absolute deviation (MAD), and 
the Shamos estimate.
}


\examples{
## A3 is used for constructing the conventional X-bar chart 
#  with the sample standard deviation.
factors.cc(n=10, factor="A3")


## Unbiasing factor for the standard deviation 
#           using the sample standard deviation.
factors.cc(n=10, factor="c4")
#  The above is the same as below:
c4.factor(n=10, method="sd")


## Unbiasing factor for the standard deviation 
#           using the sample range. 
factors.cc(n=10, factor="d2")
#  The above is the same as below:
c4.factor(n=10, method="range") 


## Table B2 in Supplement B of ASTM (1951).
char = c("A","A1","A2","c2", "B1","B2","B3","B4", "d2","d3","D1","D2","D3","D4")
nn = 2L:25L

res=NULL
for(n in nn){tmp=NULL;for(ch in char) tmp=c(tmp,factors.cc(n,ch));res=rbind(res,tmp)}
rownames(res) = paste0("n=",nn)
round(res,4)


## Table 49 in Chapter 3 of ASTM (2010).
char = c("A","A2","A3","c4", "B3","B4","B5","B6", "d2","d3","D1","D2","D3","D4")
nn = 2L:25L

res=NULL
for(n in nn){tmp=NULL;for(ch in char) tmp=c(tmp,factors.cc(n,ch));res=rbind(res,tmp)}
rownames(res) = paste0("n=",nn)
round(res,4)


## Table 50 in Chapter 3 of ASTM (2010).
char = c("E2", "E3")
nn = 2L:25L

res=NULL
for(n in nn){tmp=NULL;for(ch in char) tmp=c(tmp,factors.cc(n,ch));res=rbind(res,tmp)}
rownames(res) = paste0("n=",nn)
round(res,3)


## vignette ##
if (interactive()) vignette("factors.cc", package="rQCC")
}


\keyword{factor}
\keyword{control chart}
\author{Chanseok Park}

