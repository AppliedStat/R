\encoding{utf-8}
\name{attributes.chart.unbalanced}
\alias{acc}
\alias{print.acc}
\alias{summary.acc}

\title{Attributes control chart with balanced/unbalanced samples}

\description{Constructs the attributes control charts with balanced/unbalanced samples.}

\usage{
acc (x, n, type=c("p","np","c","u","g","h","t"), parameter, pEstimator=c("Wald","Wilson"),
     gEstimator = c("ML", "MVU"), tModel=c("E","W"), location.shift = 0,
     sigmaFactor=3, nk)
}

\arguments{
 \item{x}{a numeric vector of the number of nonconforming units.}
 \item{n}{a numeric vector of sample sizes.}
 \item{type}{a character string specifying the type of control chart.}
 \item{parameter}{a known Bernoulli parameter value for the corresponding chart.
If not known, it is estimated. For more details, refer to \code{vignette("acc", package="rQCC")}.}
 \item{pEstimator}{a method for estimating the Bernoulli parameter for \eqn{p}{p} and \eqn{np}{np} charts.}
 \item{gEstimator}{a method for estimating the Bernoulli parameter for \eqn{g}{g} and \eqn{h}{h} charts.}
 \item{tModel}{Probability model for \eqn{t}{t} chart. "E" for Exponential and "W" for Weibull.}
 \item{location.shift}{a known location shift parameter value for  \eqn{g}{g} and \eqn{h}{h} charts.}
 \item{sigmaFactor}{a factor for the standard deviation (\eqn{\sigma}{\sigma}).  
For example, the American Standard uses "\emph{3*sigma}" limits (0.27\% false alarm rate),
while the British Standard uses "\emph{3.09*sigma}" limits (0.20\% false alarm rate).}
 \item{nk}{sample size for Phase-II. If \code{nk} is missing, the average of the subsample sizes is used.}
}

\details{\code{acc} constructs the attributes control charts 
for nonconforming units (\eqn{p}{p} and \eqn{np}{np} charts)
and for nonconformities per unit (\eqn{c}{c} and \eqn{u}{u} charts).}

\value{\code{acc} returns an object of \link{class} "acc".
The function \code{summary} is used to obtain and print a summary of the results
and the function \code{plot} is used to plot the control chart.
}

\references{
ASTM (2010).
\emph{Manual on Presentation of Data and Control Chart Analysis} (STP 15-D), 8th edition. 
American Society for Testing and Materials, West Conshohocken, PA.

Kaminsky, F. C., J. C. Benneyan and R. D. Davis (1992).
Statistical Control Charts Based on a Geometric Distribution.
\emph{Journal of Quality Technology}, \bold{24}, 63-69.

Montgomery, D. C. (2013). 
\emph{Statistical Quality Control: An Modern Introduction}, 7th edition.
John Wiley & Sons, New York, NY.

Park, C. (2013).
An Improved \emph{p} Chart Based on the Wilson Interval.
\emph{Journal of Statistics and Management Systems}, \bold{16}(2-03), 201-221.\cr
\doi{10.1080/09720510.2013.777576}

Park, C. and M. Wang (2022). 
A study on the g and h control charts. 
\emph{Communications in Statistics - Theory and Methods}, \bold{To Appear}.\cr
\doi{10.1080/03610926.2022.2044492}
}

\examples{
# ==============================
# Example 1a: p chart (balanced)
# ------------------------------
# Refer to Table 31 of ASTM (2010).
x = c(1, 3, 0, 7, 2, 0, 1, 0, 8, 5, 2, 0, 1, 0, 3)
n = 400

# The conventional p chart with the balanced samples.
# Print LCL, CL, UCL.
result = acc(x, n)
print(result)

# Summary of a control chart.
summary(result)

# Plot of a control chart.
plot(result, cex.text=0.8)
text(15, 0.5, labels="p chart (with balanced sample)" )

# The p chart based on the Wilson confidence interval.
acc(x, n, pEstimator="Wilson")


# ===============================
# Example 1b: np chart (balanced)
# -------------------------------
# The data are the same as in Example 1a.
# The conventional np chart with the balanced samples.
# Print LCL, CL, UCL.
result = acc(x, n, type="np")
print(result)
summary(result)
plot(result, cex.text=0.8)
text(15, 25, labels="np chart" )

# The np chart based on the Wilson confidence interval.
acc(x, n, type="np", pEstimator="Wilson")


# ================================
# Example 2a: p chart (unbalanced)
# --------------------------------
# Refer to Table 32 of ASTM (2010).
x = c( 9, 7, 3, 9,13,14,14,10,12,14, 6,12, 7,11, 5,
       4, 2, 4, 7, 9, 7,12, 8, 5,18, 7, 8, 8,15, 3, 5)
n = c( 580, 550, 580, 640, 880, 880, 640, 550, 580, 880,
       800, 800, 580, 580, 550, 330, 330, 640, 580, 550,
       510, 640, 300, 330, 880, 880, 800, 580, 880, 880, 330)

# The conventional p chart with the unbalanced samples.
# Print LCL, CL, UCL.
result = acc(x, n, nk=880)
print(result)

# Summary of a control chart.
summary(result)

# Plot of a control chart.
plot(result, cex.text=0.8)
text(15, 0.2, labels="p chart (with unbalanced sample)" )


# ================================
# Example 2b: p chart (unbalanced)
# --------------------------------
# Refer to Table 7.4 of Montgomery (2013).
x = c(12, 8, 6, 9, 10, 12, 11, 16, 10, 6, 20, 15, 9, 8, 6, 8, 10, 7, 5, 8, 5, 8, 10, 6, 9)
n = c(100,80,80,100,110,110,100,100,90,90,110,120,120,120,110,80,80,80,90,100,100,100,100,90,90)

# The conventional p chart with the unbalanced samples.
# Print LCL, CL, UCL.
# If nk is missing, the average sample size is used.
result = acc(x, n)
print(result)

# Summary of a control chart.
summary(result)

# Plot of a control chart.
# Refer to Figure 7.8 of Montgomery (2013).
plot(result, cex.text=0.8)
text(15, 0.2, labels="p chart (with unbalanced sample)" )


# ================================
# Example 2c: p chart (unbalanced)
#             p is known 
# --------------------------------
# Refer to Table 41 of ASTM (2010).
x = c(2, 2, 1, 1, 5, 2, 0, 3, 0, 15, 7, 2, 5, 2, 0, 3, 0, 4, 8, 4)
n = c(600,1300,2000,2500,1550,2000,1550,780,260,2000,1550,950,950,950,35,330,200,600,1300,780)

# The fraction nonconforming is known as 0.0020
# The control limits at the size nk=600.
# If nk (sample size for Phase II) is unknown, the average of subsample sizes is used.
result = acc(x, n, parameter=0.002, nk=600)
summary(result)


# ===============================
# Example 3a: u chart (balanced)
# -------------------------------
# Refer to Table 33 of ASTM (2010).
x = c(17, 14, 6, 23,  5,  7, 10, 19, 29, 18, 25, 5,
       8, 11, 18, 13, 22, 6, 23, 22,  9, 15, 20, 6, 24)
n = 10

# The u chart with the balanced samples.
# Print LCL, CL, UCL.
result = acc(x, n, type="u")
print(result)

# Summary of a control chart
summary(result)

# Plot of a control chart
plot(result, cex.text=0.8)
text(13, 3, labels="u chart" )


# ================================
# Example 3b: u chart (unbalanced)
# --------------------------------
# Refer to Table 34 of ASTM (2010).
x = c(72, 38, 76, 35, 62, 81, 97, 78, 103, 56,
      47, 55, 49, 62, 71, 47, 41, 52, 128, 84)
n = c(20, 20, 40, 25, 25, 25, 40, 40, 40, 40,
      25, 25, 25, 25, 25, 20, 20, 20, 40, 40)

# The u chart with the unbalanced samples.
# Print LCL, CL, UCL.
result = acc(x, n, type="u", nk=20)
print(result)

# Summary of a control chart
summary(result)

# Plot of a control chart
plot(result, cex.text=0.8)
text(12, 3.5, labels="u chart (with unbalanced sample)" )


# ===============================
# Example 4: c chart 
# -------------------------------
# Refer to Table 35 of ASTM (2010).
x = c(0, 1, 1, 0, 2, 1, 3, 4, 5, 3, 0, 1, 1, 1, 2, 4, 0, 1, 1, 0,
      6, 4, 3, 2, 0, 0, 9,10, 8, 8, 6,14, 0, 1, 2, 4, 5, 7, 1, 3,
      3, 2, 0, 1, 5, 3, 4, 3, 5, 4, 2, 0, 1, 2, 5, 9, 4, 2, 5, 3)

# Print LCL, CL, UCL.
result = acc(x, type="c")
print(result)

# Summary of a control chart
summary(result)

# Plot of a control chart
plot(result, cex.text=0.8)
text(40, 14, labels="c chart" )


# ===============================
# Example 5: g and h charts
# -------------------------------
# Refer to Kaminsky et al. (1992).
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

# g chart with ML method.
# Print LCL, CL, UCL.
result = acc(data, type="g", location=1)
print(result)

# Summary of a control chart
summary(result)

plot(result, cex.text=0.8)

# h chart with MVU method.
acc(data, type="h", location=1, gEstimator="MVU")


# ===============================
# Example 6: g and h charts (unbalanced data)
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

result = acc(data, type="g", location=1, gEstimator="MVU", nk=5)
summary(result)
plot(result)


# ===============================
# Example 7: t charts 
# -------------------------------
x = c(0.35, 0.92, 0.59, 4.28, 0.21, 0.79, 1.75, 0.07, 3.3, 
1.7, 0.33, 0.97, 0.96, 2.23, 0.88, 0.37, 1.3, 0.4, 0.19, 1.59)

# Exponential t chart
result = acc(x, type="t", tModel="E")
summary(result)

plot(result, cex.text=0.8)
text(10, 6, labels="Exponential t chart" )


# Weibull t chart
result = acc(x, type="t", tModel="W")
summary(result)

plot(result, cex.text=0.8)
text(10, 6, labels="Weibull t chart" )

}

\author{Chanseok Park}
\keyword{control chart}
\keyword{unbalanced}
