\name{gp.plot}
\alias{gp.plot}
\title{Gumbel Probability Plot}
\description{\code{gp.plot} produces a Gumbel probability plot.}
\usage{
gp.plot(x, plot.it=TRUE, a, col.line="black", lty.line=1, 
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

\details{The Gumbel probability plot is based on taking the logarithm of the Gumbel cumulative distribution function twice.}

\value{A list with the following components:
    \item{x}{The sorted data}
    \item{y}{-log(-log(ppoints(n,a=a)))}

}

\author{Chanseok Park}

\seealso{
\code{\link{plot}}, \code{\link{qqnorm}}, \code{\link{qqplot}}, \code{\link{wp.plot}}, \code{\link{iwp.plot}}, \code{\link{ep.plot}}. \cr 
\code{\link[bsgof:bs.plot]{bs.plot}} for the Birnbaum-Saunders probability plot in package \CRANpkg{bsgof}.
}

\examples{ 
x = c(-3.16, -3.07, -2.24, -1.8, -1.48, -0.92, -0.87, -0.41, -0.06, 1.15)
# With cosmetic lines 
gp.plot(x, main="Gumbel Probability Plot", col.line="red", 
        xlab="Lifetimes of bearings", lty.line=1, pch=3)
hline = -log(-log(c( (1:5)/100, (1:9)/10)))
abline( h=hline, col=gray(0.1), lty=3, lwd=0.5 )
abline( v=seq(-4, 2,by=0.5), col=gray(0.1), lty=3, lwd=0.5 )
}

