\encoding{utf-8}
\name{robust.Quality.Control.Chart}
\alias{rcc}
\alias{rqcc}
\alias{rQCC}
\alias{print.rcc}
\alias{summary.rcc}
\title{Robust quality control chart}

\description{Constructs a robust quality control chart.}

\usage{
rcc (x, location = c("mean", "median","HL1","HL2","HL3"), 
     scale = c("sd", "range", "mad","shamos"), 
     correction = TRUE, sigma.factor = 3) 
}
\arguments{
 \item{x}{a numeric vector of observations.}
 \item{location}{a character string specifying the location estimator, must be 
one of \code{"mean"} (default), \code{"median"}, \code{"HL1"}, \code{"HL2"} and \code{"HL3"}.}
 \item{scale}{a character string specifying the scale estimator, must be 
one of \code{"sd"} (default), \code{"range"}, \code{"mad"} and \code{"shamos"}.}
 \item{correction}{A finite-sample bias correction.
\code{TRUE} adjusts a finite-sample bias correction using \code{A3.factor} function.} 
 \item{sigma.factor}{a factor for the standard deviation (\eqn{\sigma}{\sigma}).  
The American Standard uses "\emph{3*sigma}" limits (0.27\% false alarm rate),
while the British Standard uses "\emph{3.09*sigma}" limits (0.20\% false alarm
 rate).}
}

\details{
\code{rcc} constructs a Xbar-S-type control chart.
For the Xbar (location) and S (scale) parts, 
various location and scale estimates can be used to obtain the robustness.
Also, one can construct the conventional \eqn{\bar{X}}{Xbar}-\eqn{S}{S} and 
 \eqn{\bar{X}}{Xbar}-\eqn{R}{R} charts.
}

\value{\code{rcc} returns an object of \link{class} "rcc".
The function \code{summary} is used to obtain and print a summary of the results
and the function \code{plot} is used to plot the control chart.
}

\references{
Park, C., H. Kim, and M. Wang (2020).
Investigation of finite-sample properties of robust location and scale estimators. 
\emph{Communications in Statistics - Simulation and Computation}, To appear.\cr
\url{https://doi.org/10.1080/03610918.2019.1699114}

ASTM (2010).
Manual on Presentation of Data and Control Chart Analysis,
8th edition. American Society for Testing and Materials. West Conshohocken, PA.
}


\examples{
# =========
# Example A (the conventional Xbar-S chart).
# ---------
# See Example 3 in Section 3.31 of ASTM (2010).

# The data below are from Table 29 in Section 3.31 of ASTM (2010).
x1 = c(0.5005, 0.5000, 0.5008, 0.5000, 0.5005, 0.5000)
x2 = c(0.4998, 0.4997, 0.4998, 0.4994, 0.4999, 0.4998)
x3 = c(0.4995, 0.4995, 0.4995, 0.4995, 0.4995, 0.4996)
x4 = c(0.4998, 0.5005, 0.5005, 0.5002, 0.5003, 0.5004)
x5 = c(0.5000, 0.5005, 0.5008, 0.5007, 0.5008, 0.5010)
x6 = c(0.5008, 0.5009, 0.5010, 0.5005, 0.5006, 0.5009)
x7 = c(0.5000, 0.5001, 0.5002, 0.4995, 0.4996, 0.4997)
x8 = c(0.4993, 0.4994, 0.4999, 0.4996, 0.4996, 0.4997)
x9 = c(0.4995, 0.4995, 0.4997, 0.4992, 0.4995, 0.4992)
x10= c(0.4994, 0.4998, 0.5000, 0.4990, 0.5000, 0.5000)
data1 = rbind(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)

# Print LCL, CL, UCL.
result = rcc(data1, loc="mean", scale="sd") 
print(result)

# Note: Xbar-S chart is a default so that the below is the same as the above.
rcc(data1) 

# Summary of a control chart
summary(result)

# Plot a control chart
plot(result, cex.text=0.8, x.text=4.2)
abline(v=5.5, lty=1, lwd=2, col="gold")
text( c(3,8), c(0.5005, 0.5005), labels=c("Group 1", "Group 2") )


# =========
# Example B (the conventional Xbar-R chart).
# ---------
# See Example 5 in Section 3.31 of ASTM (2010).
# The data are the same as in Example A.

# Print LCL, CL, UCL.
result = rcc(data1, loc="mean", scale="range")
print(result)

# Summary of a control chart
# Note: the RE is calculated based on the unbiased estimators.
summary(result)
RE(n=6, method="range", correction=TRUE)


# Plot a control chart
plot(result, cex.text=0.8, x.text=4.2)
abline(v=5.5, lty=1, lwd=2, col="gold")
text( c(3,8), c(0.5005, 0.5005), labels=c("Group 1", "Group 2") )


# ========= 
# Example C (median-MAD chart)
# --------- 
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
data2 = rbind(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)

# Print LCL, CL, UCL.
result = rcc(data2, loc="median", scale="mad")
print(result)

# Summary of a control chart
summary(result)

# Plot a control chart
plot(result)


# =========
# Example D (HL2-Shamos chart)
# ---------
# The data are the same as in Example C.

# Print LCL, CL, UCL.
result = rcc(data2, loc="HL2", scale="shamos")
print(result)

# Summary of a control chart
summary(result)

# Plot a control chart
plot(result)
}

\keyword{control chart}
\keyword{robust}
\author{Chanseok Park}
