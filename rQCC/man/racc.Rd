\encoding{utf-8}
\name{robust.attributes.chart.unbalanced}
\alias{racc}
\alias{print.racc}
\alias{summary.racc}

\title{Robust attributes control charts with balanced/unbalanced samples}

\description{Constructs the robust g and h attributes control charts with balanced/unbalanced samples.}

\usage{
racc (x, gamma, type=c("g","h","t"), parameter, gEstimator=c("cdf", "MM"), 
      tModel=c("E", "W"), location.shift = 0, sigmaFactor=3, nk)
}

\arguments{
 \item{x}{a numeric vector of the number of nonconforming units.}
 \item{gamma}{a numeric value for a inlier proportion.
gamma should be between 0 and 1 (smaller value means more trimming). }
 \item{type}{a character string specifying the type of control chart.}
 \item{parameter}{a known Bernoulli parameter value for the \eqn{g}{g} and \eqn{h}{h} charts.
If not known, it is estimated. For more details, refer to \code{vignette("racc", package="rQCC")}.}
 \item{gEstimator}{a method for estimating the Bernoulli parameter for \eqn{g}{g} and \eqn{h}{h} charts.
"cdf" is based on the memoryless property and "MM" is based on the truncated geometric distribution.
}
 \item{tModel}{Probability model for \eqn{t}{t} chart. "E" for Exponential and "W" for Weibull.}
 \item{location.shift}{a known location shift parameter value for  \eqn{g}{g} and \eqn{h}{h} charts.}
 \item{sigmaFactor}{a factor for the standard deviation (\eqn{\sigma}{\sigma}).  
For example, the American Standard uses "\emph{3*sigma}" limits (0.27\% false alarm rate),
while the British Standard uses "\emph{3.09*sigma}" limits (0.20\% false alarm rate).}
 \item{nk}{sample size for Phase-II. If \code{nk} is missing, the average of the subsample sizes is used.}
}

\details{\code{racc} constructs the attributes control charts 
for nonconforming units (\eqn{p}{p} and \eqn{np}{np} charts)
and for nonconformities per unit (\eqn{c}{c} and \eqn{u}{u} charts).}

\value{\code{racc} returns an object of \link{class} "racc".
The function \code{summary} is used to obtain and print a summary of the results
and the function \code{plot} is used to plot the control chart.
}

\references{
Park, C., L. Ouyang, and M. Wang (2021).
Robust g-type quality control charts for monitoring nonconformities. 
\emph{Computers and Industrial Engineering}, \bold{162}, 107765.

Kaminsky, F. C., J. C. Benneyan and R. D. Davis (1992).
Statistical Control Charts Based on a Geometric Distribution.
\emph{Journal of Quality Technology}, \bold{24}, 63-69.
}

\examples{
# ===============================
# Example 1: g and h charts
# -------------------------------
# Refer to Kaminsky et al. (1992) and Table 2 of Park, et al. (2021).
tmp = c(
11,  2,  8,  2, 4,   1,  1, 11,  2, 1,   1,  7,  1,  1, 9, 
 5,  1,  3,  6, 5,  13,  2,  3,  3, 4,   3,  2,  6,  1, 5,  
 2,  2,  8,  3, 1,   1,  3,  4,  6, 5,   2,  8,  1,  1, 4,  
13, 10, 15,  5, 2,   3,  6,  1,  5, 8,   9,  1, 18,  3, 1,  
 3,  7, 14,  3, 1,   7,  7,  1,  8, 1,   4,  1,  6,  1, 1, 
 1, 14,  2,  3, 7,  19,  9,  7,  1, 8,   5,  1,  1,  6, 1,  
 9,  5,  6,  2, 2,   8, 15,  2,  3, 3,   4,  7, 11,  4, 6,  
 7,  5,  1, 14, 8,   3,  3,  5, 21,10,  11,  1,  6,  1, 2,  
 4,  1,  2, 11, 5,   3,  5,  4, 10, 3,   1,  4,  7,  3, 2, 
 3,  5,  4,  2, 3,   5,  1,  4, 11,17,   1, 13, 13,  2, 1)  
data = matrix(tmp, byrow=TRUE, ncol=5)

# g chart with cdf (trimming) method.
# Print LCL, CL, UCL.
result = racc(data, gamma=0.9, type="g", location=1)
print(result)

# Summary of a control chart
summary(result)

plot(result, cex.text=0.8)

# h chart with MM (truncated geometric) method.
racc(data, gamma=0.9, type="h", location=1, gEstimator="MM")


# ===============================
# Example 2: g and h charts (unbalanced data)
# -------------------------------
x1 = c(11, 2,  8,  2, 4)
x2 = c(1,  1, 11,  2, 1)
x3 = c(1,  7,  1)
x4 = c(5,  1,  3,  6, 5)
x5 = c(13, 2,  3,  3)
x6 = c(3,  2,  6,  1, 5)
x7 = c(2,  2,  8,  3, 1)
x8 = c(1,  3,  4,  6, 5)
x9 = c(2,  8,  1,  1, 4)
data = list(x1, x2, x3, x4, x5, x6, x7, x8, x9)

result = racc(data, gamma=0.9, type="g", location=1, gEstimator="cdf", nk=5)
summary(result)
plot(result)


# ===============================
# Example 3: t charts 
# -------------------------------
x = c(0.35, 0.92, 0.59, 4.28, 0.21, 0.79, 1.75, 0.07, 3.3,
1.7, 0.33, 0.97, 0.96, 2.23, 0.88, 0.37, 1.3, 0.4, 0.19, 1.59)

# Exponential t chart
result = racc(x, type="t", tModel="E")
summary(result)

plot(result, cex.text=0.8)
text(10, 6, labels="Robust exponential t chart" )


# Weibull t chart
result = racc(x, type="t", tModel="W")
summary(result)

plot(result, cex.text=0.8)
text(10, 5.5, labels="Robust Weibull t chart" )

}

\author{Chanseok Park}
\keyword{control chart}
\keyword{robust}
\keyword{unbalanced}
