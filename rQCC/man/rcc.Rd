\name{robust.control.chart}
\alias{rcc}
\title{Robust statistical control chart}
\usage{
rcc (x, location = c("median","HL1","HL2","HL3","mean"), scale = c("mad","shamos","sd"), 
     bias.correction = TRUE, sigma.factor = 3) 
}
\arguments{
 \item{x}{a numeric vector of observations.}
 \item{location}{a character string specifying the location estimator, must be 
one of "median" (default), "HL1", "HL2" and "HL3."}
 \item{scale}{a character string specifying the scale estimator, must be 
one of "mad" (default), "shamos" and "sd."}
 \item{bias.correction}{A finite-sample bias correction.
\code{TRUE} adjusts a finite-sample bias correction using \code{A3.factor} function.} 
 \item{sigma.factor}{a factor for the standard deviation (\eqn{\sigma}{\sigma}).  
The American Standard uses "\emph{3*sigma}" limits (0.27\% false alarm rate),
while the British Standard uses "\emph{3.09*sigma}" limits (0.20\% false alarm
 rate).}
}

\description{Constructs robust statistical control chart using 
the median and Hodges-Lehmann (HL1, HL2, HL3) estimators, 
and the median absolute deviation (MAD) and Shamos estimators.
Note that a conventional X-bar chart is also constructed by 
\code{rcc(x, location="mean", scale="sd")}.
}

\value{\code{rcc} provides the center line and the lower and upper control limits.}

\references{
Park, C., H. Kim, and M. Wang (2019).
Finite-sample properties of robust location and scale estimators. 
\emph{arXiv:1908.00462}.
}


\examples{
# Key in data (10 subgroups with 5 observations each).
x1 = c(4.13, 5.97, 5.54, 5.26, 5.73)
x2 = c(4.11, 6.27, 5.84, 5.17, 5.24)
x3 = c(4.51, 4.83, 5.67, 5.99, 5.85)
x4 = c(3.77, 3.78, 5.69, 4.31, 4.42)
x5 = c(4.83, 5.15, 3.83, 3.48, 2.52)
x6 = c(4.25, 3.66, 5.23, 5.96, 4.75)
x7 = c(4.93, 5.17, 6.44, 4.19, 4.42)
x8 = c(5.28, 4.83, 7.94, 4.89, 3.51)
x9 = c(4.58, 5.60, 5.84, 6.22, 4.54)
x10= c(5.14, 4.31, 3.23, 6.95, 5.35)
data = rbind(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)

# Print LCL, CL, UCL.
result = rcc(data)
print(result)

# Plot a control chart
xbar = rowMeans(data)
plot(xbar, type="b", ylim=range(result,xbar) ) 
abline(h=result, lty=2, col="red")
}

\keyword{control chart}
\keyword{robust}
\author{Chanseok Park and Min Wang}
