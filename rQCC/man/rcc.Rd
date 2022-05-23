\encoding{utf-8}
\name{robust.control.chart.unbalanced}
\alias{rQCC}
\alias{rcc}
\alias{print.rcc}
\alias{summary.rcc}

\title{Robust quality control chart with balanced/unbalanced samples}

\description{Constructs the robust control charts with balanced/unbalanced samples.}

\usage{
rcc (x, location = c("mean", "median", "HL1", "HL2", "HL3"), 
     scale = c("sd", "range", "mad", "shamos"), type = c("Xbar", "S", "R"), 
     poolLoc = c("A", "B", "C"), poolScale = c("A", "B", "C"), sigmaFactor=3, nk)
}

\arguments{
 \item{x}{a numeric matrix or list of vectors. Each row or vector contains a sub-sample.}
 \item{location}{a character string specifying the location estimator, must be 
one of \code{"mean"} (default), \code{"median"}, \code{"HL1"}, \code{"HL2"} and \code{"HL3"}.}
 \item{scale}{a character string specifying the scale estimator, must be 
one of \code{"sd"} (default), \code{"range"}, \code{"mad"} and \code{"shamos"}.}
 \item{type}{a character string specifying the type of control chart.}
 \item{poolLoc}{pooling type for location.}
 \item{poolScale}{pooling type for scale.}
 \item{sigmaFactor}{a factor for the standard deviation (\eqn{\sigma}{\sigma}).  
For example, the American Standard uses "\emph{3*sigma}" limits (0.27\% false alarm rate),
while the British Standard uses "\emph{3.09*sigma}" limits (0.20\% false alarm rate).}
 \item{nk}{sample size for Phase-II. If \code{nk} is missing, the average sample size is used.}
}

\details{
\code{rcc} constructs the robust X-bar, S and R control charts.
Using various robust location and scale estimators, one can construct the robust control charts.
For more details on the control charts, refer to \code{vignette("rcc", package="rQCC")}.

Note that the location and scale estimators used in this fuction are all unbiased.
For more details on how to pool the location and scale estimators, refer to 
\code{vignette("pooledEstimator", package="rQCC")}.

In addition, one can also construct the conventional \eqn{\bar{X}}{X-bar} chart with \eqn{S}{S} and 
 \eqn{\bar{X}}{X-bar} chart with \eqn{R}{R}.
}

\value{\code{rcc} returns an object of \link{class} "rcc".
The function \code{summary} is used to obtain and print a summary of the results
and the function \code{plot} is used to plot the control chart.
}

\references{
Park, C., H. Kim, and M. Wang (2022).
Investigation of finite-sample properties of robust location and scale estimators. 
\emph{Communications in Statistics - Simulation and Computation}, 
\bold{51}, 2619-2645. \cr 
\doi{10.1080/03610918.2019.1699114}

ASTM (2010).
\emph{Manual on Presentation of Data and Control Chart Analysis} (STP 15-D),
8th edition. 
American Society for Testing and Materials, West Conshohocken, PA.

Ryan (2000). 
\emph{Statistical Methods For Quality Improvement}, 2nd edition.
John Wiley & Sons, New York, NY.

Montgomery (2013). 
\emph{Statistical Quality Control: An Modern Introduction}, 7th edition. 
John Wiley & Sons, New York, NY.
}


\examples{
###############
# X-bar chart #
###############

# ========== 
# Example 1a 
# ---------- 
# The conventional X-bar chart with the standard deviation. 
# Refer to Example 3 in Section 3.31 of ASTM (2010). 

# The data below are from Table 29 in Section 3.31 of ASTM (2010). 
# Each subgroup has a sample of size n=6. There are m=10 subgroups.
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
RE(n=6, estimator="sd", correction=TRUE)

# The above limits are also calculated as 
A3 = factors.cc(n=6, "A3")

xbarbar = mean(data1)
 # xbarbar = mean(unlist(data1)) # for list

sbar = mean( apply(data1, 1, sd) )
c(xbarbar-A3*sbar, xbarbar, xbarbar+A3*sbar)

# Plot a control chart
plot(result, cex.text=0.8)
abline(v=5.5, lty=1, lwd=2, col="gold")
text( c(3,8), c(0.5005, 0.5005), labels=c("Group 1", "Group 2") )

# ==========
# Example 1b 
# ----------
# The conventional X-bar chart with the range.
# Refer to Example 5 in Section 3.31 of ASTM (2010).
# The data are the same as in Example 1a.

# Print LCL, CL, UCL.
# The range is used for the scale estimator.
result = rcc(data1, loc="mean", scale="range")
print(result)

# Summary of a control chart
# Note: the RE is calculated based on the unbiased estimators.
summary(result)
RE(n=6, estimator="range", correction=TRUE)

# The above limits are also calculated as 
A2 = factors.cc(n=6, "A2")

xbarbar = mean(data1)
 # xbarbar = mean(unlist(data1)) # for list

Rbar = mean( apply(data1, 1, function(x) {diff(range(x))}) )
 # Rbar = mean( apply(sapply(data1,range),2,diff) ) # for list

c(xbarbar-A2*Rbar, xbarbar, xbarbar+A2*Rbar)

# Plot a control chart
plot(result, cex.text=0.8)
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
rcc(data2, loc="median", scale="mad")

# Note: the RE is calculated based on the unbiased estimators.
RE(n=4, estimator="median", correction=TRUE)

# ==========
# Example 1d 
# ----------
# The HL2-Shamos chart.
# The data are the same as in Example 1c.

# Print LCL, CL, UCL.
# The HL2 (location) and Shamos (scale) are used.
rcc(data2, loc="HL2", scale="shamos")

# Note: the RE is calculated based on the unbiased estimators.
RE(n=4, estimator="HL2", correction=TRUE)


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
# The standard deviaion (default) is used for the scale estimator.
result = rcc(data1, type="S")
print(result)

# The above limits are also calculated as 
B3 = factors.cc(n=6, "B3")
B4 = factors.cc(n=6, "B4")
sbar = mean( apply(data1, 1, sd) )
c(B3*sbar, sbar, B4*sbar)

# Plot a control chart
plot(result, cex.text=0.8)
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

# Plot a control chart
plot(result, cex.text=0.8)
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
# The range is used for the scale estimator.
# Unlike the S chart, scale="range" is not a default. 
# Thus, for the conventional R chart, use the option (scale="range") as below.
result = rcc(data1, scale="range", type="R")
print(result)

# The above limits are also calculated as 
D3 = factors.cc(n=6, "D3")
D4 = factors.cc(n=6, "D4")

Rbar = mean( apply(data1, 1, function(x) {diff(range(x))}) )
 # Rbar = mean( apply(sapply(data1,range),2,diff) ) # for list

c(D3*Rbar, Rbar, D4*Rbar)

# Plot a control chart
plot(result, cex.text=0.8)
abline(v=5.5, lty=1, lwd=2, col="gold")
text( c(3,8), c(0.00135, 0.00135), labels=c("Group 1", "Group 2") )


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

# Plot a control chart
plot(result, cex.text=0.8)
abline(v=5.5, lty=1, lwd=2, col="gold")
text( c(3,8), c(0.00135, 0.00135), labels=c("Group 1", "Group 2") )


###################
# Unbalanced Data #
###################

# ==========
# Example 4
# ----------
# Refer to Example 4 in Section 3.31 of ASTM (2010).
# Data set is from Table 30 in Section 3.31 of ASTM (2010).
 x1 = c( 73, 73, 73, 75, 75)
 x2 = c( 70, 71, 71, 71, 72)
 x3 = c( 74, 74, 74, 74, 75)
 x4 = c( 70, 70, 70, 72, 73)
 x5 = c( 70, 70, 70, 70, 70)
 x6 = c( 65, 65, 66, 69, 70)
 x7 = c( 72, 72, 74, 76)
 x8 = c( 69, 70, 71, 73, 73)
 x9 = c( 71, 71, 71, 71, 72)
x10 = c( 71, 71, 71, 71, 72)
x11 = c( 71, 71, 72, 72, 72)
x12 = c( 70, 71, 71, 72, 72)
x13 = c( 73, 74, 74, 75, 75)
x14 = c( 74, 74, 75, 75, 75)
x15 = c( 72, 72, 72, 73, 73)
x16 = c( 75, 75, 75, 76)
x17 = c( 68, 69, 69, 69, 70)
x18 = c( 71, 71, 72, 72, 73)
x19 = c( 72, 73, 73, 73, 73)
x20 = c( 68, 69, 70, 71, 71)
x21 = c( 69, 69, 69, 69, 69)

# For unbalanced data set, use list.
data = list(x1, x2,  x3, x4,  x5,  x6,  x7,  x8,  x9, x10,
            x11,x12,x13,x14, x15, x16, x17, x18, x19, x20, x21)

# Xbar chart (witn nk=5)
rcc(data, nk=5)

# S chart (witn nk=4)
rcc(data, type="S", nk=4)

# ==========
# Example 5a
# ----------
# Data set is from Example 6.4 of Montgomery (2013)
#    Statistical Quality Control (7th ed), Wiley.
# Data set for Phase I
 x1 = c(74.030, 74.002, 74.019, 73.992, 74.008)
 x2 = c(73.995, 73.992, 74.001)
 x3 = c(73.988, 74.024, 74.021, 74.005, 74.002)
 x4 = c(74.002, 73.996, 73.993, 74.015, 74.009)
 x5 = c(73.992, 74.007, 74.015, 73.989, 74.014)
 x6 = c(74.009, 73.994, 73.997, 73.985)
 x7 = c(73.995, 74.006, 73.994, 74.000)
 x8 = c(73.985, 74.003, 73.993, 74.015, 73.988)
 x9 = c(74.008, 73.995, 74.009, 74.005)
x10 = c(73.998, 74.000, 73.990, 74.007, 73.995)
x11 = c(73.994, 73.998, 73.994, 73.995, 73.990)
x12 = c(74.004, 74.000, 74.007, 74.000, 73.996)
x13 = c(73.983, 74.002, 73.998) 
x14 = c(74.006, 73.967, 73.994, 74.000, 73.984)
x15 = c(74.012, 74.014, 73.998) 
x16 = c(74.000, 73.984, 74.005, 73.998, 73.996)
x17 = c(73.994, 74.012, 73.986, 74.005) 
x18 = c(74.006, 74.010, 74.018, 74.003, 74.000)
x19 = c(73.984, 74.002, 74.003, 74.005, 73.997)
x20 = c(74.000, 74.010, 74.013) 
x21 = c(73.982, 74.001, 74.015, 74.005, 73.996)
x22 = c(74.004, 73.999, 73.990, 74.006, 74.009)
x23 = c(74.010, 73.989, 73.990, 74.009, 74.014)
x24 = c(74.015, 74.008, 73.993, 74.000, 74.010)
x25 = c(73.982, 73.984, 73.995, 74.017, 74.013) 

# For unbalanced data set, use list.
data = list(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, 
           x16, x17, x18, x19, x20, x21, x22, x23, x24, x25)

# Xbar chart (witn nk=5)
rcc(data, nk=5)

# S chart (witn nk=5)
rcc(data, type="S", nk=5)


# ==========
# Example 5b
# ----------
# With contaminated data set. 
# Two contaminated observations are added 
#   in the first subgroup (70.5, 77.0) of Example 5a.
datan = data
datan[[1]] = c(data[[1]], 70.5, 77.0)

# Xbar chart with non-robust estimators
rcc(datan, nk=5)

# Xbar chart with robust estimators (median and mad)
rcc(datan, loc="median", sc="mad", nk=5)
}

\author{Chanseok Park}
\keyword{control chart}
\keyword{robust}
\keyword{unbalanced}
