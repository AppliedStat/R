\name{lp.plot}
\alias{lp.plot}
\title{Gumbel Probability Plot}
\description{\code{lp.plot} produces a Gumbel probability plot.}
\usage{
lp.plot(x, plot.it=TRUE, a, col.line="black", lty.line=1,
        xlim=NULL, ylim=NULL, main=NULL, sub=NULL, xlab=NULL, ylab="Probability", ...)
}

\arguments{
 \item{x}{a numeric vector of data values. Missing values are allowed.}
 \item{plot.it}{logical. Should the result be plotted?}
 \item{a}{the offset fraction to be used; typically in (0,1). See \code{\link{ppoints}}.}
 \item{col.line}{the color of the straight line.}
 \item{lty.line}{the line type of the straight line.}
 \item{xlim}{the x limits of the plot.}
 \item{ylim}{the y limits of the plot.}
 \item{main}{a main title for the plot, see also \code{\link{title}}.}
 \item{sub}{a sub title for the plot.}
 \item{xlab}{a label for the x axis, defaults to a description of x.}
 \item{ylab}{a label for the y axis, defaults to "Probability".}
 \item{...}{graphical parameters.}
}

\details{The Gumbel probability plot is based on taking the logarithm of the Laplace cumulative distribution function.}

\value{A list with the following components:
    \item{x}{The sorted data}
    \item{y}{p5=ppoints(length(x))-0.5; y=-sign(p5)*log(1-2*abs(p5))}

}
\references{
Bain, L. J., and Englehardt, M. (1973). 
Interval Estimation for the Two-Parameter Double Exponential Distribution. 
\emph{Technometrics}, \bold{15}, 875-887.
}
\author{Chanseok Park}

\seealso{
\code{\link{plot}}, \code{\link{qqnorm}}, \code{\link{qqplot}}. \cr
}

\examples{
# See Bain and Englehardt (1973).
x = x = c( 1.96, 1.97, 3.60, 3.80, 4.79, 5.66, 5.76, 5.78, 6.27, 6.30, 6.76,
 7.65, 7.84, 7.99, 8.51, 9.18,10.13,10.24,10.25,10.43,11.45,11.48,
11.75,11.81,12.34,12.78,13.06,13.29,13.98,14.18,14.40,16.22,17.06 )

# With cosmetic lines 
lp.plot(x, main="Gumbel Probability Plot", col.line="red",
        xlab="Lifetimes of bearings", lty.line=1, pch=3)
p5 = c((1:5)/100, (1:9)/10, 0.95, 0.99) - 0.5 
hline = -sign(p5)*log(1-2*abs(p5))
abline( h=hline, col=gray(0.1), lty=3, lwd=0.5 )
abline( v=seq(0, 20 ,by=1), col=gray(0.1), lty=3, lwd=0.5 )
}


