\name{bs.plot}
\alias{bs.plot}
\title{Birnbaum-Saunders Probability Plot}
\description{\code{bs.plot} produces a Birnbaum-Saunders probability plot.}
\usage{
bs.plot(x, plot.it=TRUE, a, col.line="black", lty.line=1,  
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

\details{The Birnbaum-Saunders probability plot is based on the linearization proposed by Chang and Tang (1994).}

\value{A list with the following components:
    \item{x}{The sorted data}
    \item{w}{sqrt(x)*qnorm(p)}
}

\references{Chang, D. S and Tang, L. C. (1994). Graphical analysis for Birnbaum-Saunders distribution.
Microelectronics Reliability 34: 17-22. \cr 

Birnbaum, Z. W. and Saunders, S. C. (1969). Estimation for a Family of Life Distributions with Applications to Fatigue.  J. Appl. Probab. 6(2): 328-347.
}

\author{Chanseok Park}

\seealso{
\code{\link{qqnorm}}, \code{\link{qqplot}}.

\code{\link[weibullness:wp.plot]{wp.plot}} for the Weibull probability plot in package \CRANpkg{weibullness}.
}

\examples{ 
# Data set from Birnbaum and Saunders (1969).
attach(BSdata)
data = psi21k
bs.plot(data)

# Adding cosmetic lines
bs.plot(data, main="BS probability plot", lty.line=2, pch=3, col.line="red")

ticklabels=c(0.01, seq(0.1,0.9,by=0.1), seq(0.91,0.99,by=0.01) )
qn = quantile(data, probs=ticklabels) 
ticksat= qnorm(ticklabels)* sqrt( qn )
hline =  qnorm( ticklabels  ) *  sqrt( qn )

abline( h=hline,    col=gray(0.5), lty=3, lwd=0.6 )
abline( v=  seq(0, 2500, by=100), col=gray(0.5), lty=3, lwd=0.5 )

abline( h= qnorm(0.5)*sqrt(median(data)), col=gray(0.1), lty=1, lwd=0.6 )
abline( v= median(data), col=gray(0.1), lty=1, lwd=0.6 )
}

