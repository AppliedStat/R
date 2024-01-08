\name{iwp.plot}
\alias{iwp.plot}
\title{Inverse Weibull Probability Plot}
\description{\code{iwp.plot} produces an inverse Weibull probability plot.}
\usage{
iwp.plot(x, plot.it=TRUE, a, col.line="black", lty.line=1, 
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

\details{The inverse Weibull probability plot is based on taking the logarithm of the inverse Weibull cumulative distribution function twice.}

\value{A list with the following components:
    \item{x}{The sorted data}
    \item{y}{-log(-log(ppoints(n,a=a)))}

}

\author{Chanseok Park}

\seealso{
\code{\link{plot}}, \code{\link{qqnorm}}, \code{\link{qqplot}}, \code{\link{wp.plot}}.
}

\examples{ 
x = c(0.38, 0.41, 1.2, 0.52, 0.69, 0.89, 0.67, 1.59, 0.55, 0.59)
iwp.plot( x )

# With cosmetic lines 
iwp.plot(x, main="Inverse Weibull Probability Plot", col.line="red", 
        xlab="Lifetimes of bearings", lty.line=1, pch=3)
hline = -log(-log( c( (1:5)/100, (1:9)/10)  ))
abline( h=hline, col=gray(0.1), lty=3, lwd=0.5 )
abline( v=  seq(0, 2,by=0.1), col=gray(0.1), lty=3, lwd=0.5 )
}

