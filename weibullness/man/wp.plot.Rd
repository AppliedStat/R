\name{wp.plot}
\alias{wp.plot}
\title{Weibull Probability Plot}
\description{\code{wp.plot} produces a Weibull probability plot.}
\usage{
wp.plot(x, plot.it=TRUE, a, col.line="black", lty.line=1, 
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

\details{The Weibull probability plot is based on taking the logarithm of the Weibull cumulative distribution function twice}

\value{A list with the following components:
    \item{x}{The sorted data}
    \item{y}{log(-log(1-ppoints(n,a=a)))}

}

\author{Chanseok Park}

\seealso{
\code{\link{plot}}, \code{\link{qqnorm}}, \code{\link{qqplot}}.
}

\examples{ 
attach(Wdata)

wp.plot( bearings )

# With cosmetic lines 
wp.plot(bearings, main="Weibull Probability Plot", col.line="red", 
        xlab="Lifetimes of bearings", lty.line=1, pch=3)
hline = log(-log(1- c( (1:5)/100, (1:9)/10)  ))
abline( h=hline, col=gray(0.1), lty=3, lwd=0.5 )
abline( v=  seq(15, 200,by=5), col=gray(0.1), lty=3, lwd=0.5 )
}

