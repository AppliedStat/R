\encoding{utf-8}
\name{robust.control.chart}
\alias{rcc}
\alias{rQCC}
\alias{print.rcc}
\alias{summary.rcc}
\title{Robust quality control chart}

\description{Constructs a robust control chart.}

\usage{
rcc (x, location = c("mean","median","HL1","HL2","HL3"), 
     scale = c("sd", "range","mad","shamos"), type = c("Xbar","S","R"), sigma.factor = 3) 
}
\arguments{
 \item{x}{a numeric vector of observations.}
 \item{location}{a character string specifying the location estimator, must be 
one of \code{"mean"} (default), \code{"median"}, \code{"HL1"}, \code{"HL2"} and \code{"HL3"}.}
 \item{scale}{a character string specifying the scale estimator, must be 
one of \code{"sd"} (default), \code{"range"}, \code{"mad"} and \code{"shamos"}.}
 \item{type}{a character string specifying the type of control chart.}
 \item{sigma.factor}{a factor for the standard deviation (\eqn{\sigma}{\sigma}).  
For example, the American Standard uses "\emph{3*sigma}" limits (0.27\% false alarm rate),
while the British Standard uses "\emph{3.09*sigma}" limits (0.20\% false alarm
 rate).}
}

\details{
\code{rcc} constructs a robust X-bar control chart.
Using various robust location and scale estimates, one can construct a robust X-bar chart.
The location and scale estimates used in this fuction are all unbiased.
In addition, one can also construct the conventional \eqn{\bar{X}}{X-bar} chart with \eqn{S}{S} and 
 \eqn{\bar{X}}{X-bar} chart with \eqn{R}{R}.
For more details, see the vignette.
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
\emph{Manual on Presentation of Data and Control Chart Analysis} (STP 15-D),
8th edition. 
American Society for Testing and Materials, West Conshohocken, PA.

Ryan (2000). 
\emph{Statistical Methods For Quality Improvement}, 2nd edition.
John Wiley & Sons, New York, NY.
}


\examples{
##############
# X-bar chart #
##############

# ========== 
# Example 1a 
# ---------- 
# The conventional X-bar chart with the standard deviation. 
# Refer to Example 3 in Section 3.31 of ASTM (2010). 

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
# The mean and standard deviation are used.
result = rcc(data1, loc="mean", scale="sd", type="Xbar") 
print(result)

# Note: X-bar chart is a default with the mean and sd
#       so the below is the same as the above.
rcc(data1) 

# Summary of a control chart
summary(result)
RE(n=6, method="sd", correction=TRUE)

# The above limits are also calculated as 
A3 = factors.cc(n=6, "A3")
xbarbar = mean(data1)
sbar = mean( apply(data1, 1, sd) )
c(xbarbar-A3*sbar, xbarbar, xbarbar+A3*sbar)

# Plot a control chart
plot(result, cex.text=0.8, x.text=4.1)
abline(v=5.5, lty=1, lwd=2, col="gold")
text( c(3,8), c(0.5005, 0.5005), labels=c("Group 1", "Group 2") )

# ==========
# Example 1b 
# ----------
# The conventional X-bar chart with the range.
# Refer to Example 5 in Section 3.31 of ASTM (2010).
# The data are the same as in Example 1a.

# Print LCL, CL, UCL.
# The range is used for the scale estimate.
result = rcc(data1, loc="mean", scale="range")
print(result)

# Summary of a control chart
# Note: the RE is calculated based on the unbiased estimators.
summary(result)
RE(n=6, method="range", correction=TRUE)

# The above limits are also calculated as 
A2 = factors.cc(n=6, "A2")
xbarbar = mean(data1)
Rbar = mean( apply(data1, 1, function(x) {diff(range(x))}) )
c(xbarbar-A2*Rbar, xbarbar, xbarbar+A2*Rbar)

# Plot a control chart
plot(result, cex.text=0.8, x.text=4.2)
abline(v=5.5, lty=1, lwd=2, col="gold")
text( c(3,8), c(0.5005, 0.5005), labels=c("Group 1", "Group 2") )

# ==========
# Example 1c 
# ----------
# The median-MAD chart.
# Refer to Table 4.2 in Section 4.7 of Ryan (2000).
# Data: 20 subgroups with 4 observations each.
tmp = c(
72, 84, 79, 49, 56, 87, 33, 42, 55, 73, 22, 60, 44, 80, 54, 74,
97, 26, 48, 58, 83, 89, 91, 62, 47, 66, 53, 58, 88, 50, 84, 69,
57, 47, 41, 46, 13, 10, 30, 32, 26, 39, 52, 48, 46, 27, 63, 34,
49, 62, 78, 87, 71, 63, 82, 55, 71, 58, 69, 70, 67, 69, 70, 94,
55, 63, 72, 49, 49, 51, 55, 76, 72, 80, 61, 59, 61, 74, 62, 57 )
data2 = matrix(tmp, ncol=4, byrow=TRUE)

# Print LCL, CL, UCL.
# The median (location) and MAD (scale) are used.
result = rcc(data2, loc="median", scale="mad")
print(result)

# Summary of a control chart
summary(result)
# Note: the RE is calculated based on the unbiased estimators.
RE(n=4, method="median", correction=TRUE)
RE(n=4, method="mad", correction=TRUE)

# Plot a control chart
plot(result)

# ==========
# Example 1d 
# ----------
# The HL2-Shamos chart.
# The data are the same as in Example 1c.

# Print LCL, CL, UCL.
# The HL2 (location) and Shamos (scale) are used.
result = rcc(data2, loc="HL2", scale="shamos")
print(result)

# Summary of a control chart
summary(result)
# Note: the RE is calculated based on the unbiased estimators.
RE(n=4, method="HL2", correction=TRUE)
RE(n=4, method="shamos", correction=TRUE)

# Plot a control chart
plot(result)


############
# S chart  #
############

# ==========
# Example 2a 
# ----------
# The conventional S chart with the standard deviation.
# Refer to Example 3 in Section 3.31 of ASTM (2010). 
# The data are the same as in Example 1a.

# Print LCL, CL, UCL.
# The standard deviaion (default) is used for the scale estimate.
result = rcc(data1, type="S")
print(result)

# Summary of a control chart
# Note: the RE is calculated based on the unbiased estimators.
summary(result)

# The above limits are also calculated as 
B3 = factors.cc(n=6, "B3")
B4 = factors.cc(n=6, "B4")
sbar = mean( apply(data1, 1, sd) )
c(B3*sbar, sbar, B4*sbar)

# Plot a control chart
plot(result, cex.text=0.8, x.text=4.1)
abline(v=5.5, lty=1, lwd=2, col="gold")
text( c(3,8), c(0.0005, 0.0005), labels=c("Group 1", "Group 2") )


# ==========
# Example 2b
# ----------
# The S-type chart with the MAD.
# The data are the same as in Example 2a.

# Print LCL, CL, UCL.
# The mad (scale) are used.
result = rcc(data1, scale="mad", type="S")
print(result)

# Summary of a control chart
# Note: the RE is calculated based on the unbiased estimators.
summary(result)

# Plot a control chart
plot(result, cex.text=0.8, x.text=4.1)
abline(v=5.5, lty=1, lwd=2, col="gold")
text( c(3,8), c(0.00045, 0.00045), labels=c("Group 1", "Group 2") )


############
# R chart  #
############

# ==========
# Example 3a 
# ----------
# The conventional R chart with the range.
# Refer to Example 5 in Section 3.31 of ASTM (2010). 
# The data are the same as in Example 1a.

# Print LCL, CL, UCL.
# The range is used for the scale estimate.
# Unlike the S chart, scale="range" is not a default. 
# Thus, for the conventional R chart, use the option (scale="range") as below.
result = rcc(data1, scale="range", type="R")
print(result)

# Summary of a control chart
# Note: the RE is calculated based on the unbiased estimators.
summary(result)

# The above limits are also calculated as 
D3 = factors.cc(n=6, "D3")
D4 = factors.cc(n=6, "D4")
Rbar = mean( apply(data1, 1, function(x) {diff(range(x))}) )
c(D3*Rbar, Rbar, D4*Rbar)

# Plot a control chart
plot(result, cex.text=0.8, x.text=4.1)
abline(v=5.5, lty=1, lwd=2, col="gold")
text( c(3,8), c(0.0014, 0.0014), labels=c("Group 1", "Group 2") )


# ==========
# Example 3b
# ----------
# The R-type chart with the Shamos.
# Refer to Example 5 in Section 3.31 of ASTM (2010). 
# The data are the same as in Example 3a.

# Print LCL, CL, UCL.
# The mad (scale) are used.
result = rcc(data1, scale="shamos", type="R")
print(result)

# Summary of a control chart
# Note: the RE is calculated based on the unbiased estimators.
summary(result)

# Plot a control chart
plot(result, cex.text=0.8, x.text=4.1)
abline(v=5.5, lty=1, lwd=2, col="gold")
text( c(3,8), c(0.0014, 0.0014), labels=c("Group 1", "Group 2") )


############
# vignette #
############

vignette("rcc", package="rQCC")

}

\keyword{control chart}
\keyword{robust}
\author{Chanseok Park}
