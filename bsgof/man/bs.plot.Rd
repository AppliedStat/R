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
}

\examples{ 
# Data set from Birnbaum and Saunders (1969).
data = c(0.37, 0.706, 0.716, 0.746, 0.785, 0.797, 0.844, 0.855, 0.858,
0.886, 0.886, 0.93, 0.96, 0.988, 0.99, 1, 1.01, 1.016, 1.018,
1.02, 1.055, 1.085, 1.102, 1.102, 1.108, 1.115, 1.12, 1.134,
1.14, 1.199, 1.2, 1.2, 1.203, 1.222, 1.235, 1.238, 1.252, 1.258,
1.262, 1.269, 1.27, 1.29, 1.293, 1.3, 1.31, 1.313, 1.315, 1.33,
1.355, 1.39, 1.416, 1.419, 1.42, 1.42, 1.45, 1.452, 1.475, 1.478,
1.481, 1.485, 1.502, 1.505, 1.513, 1.522, 1.522, 1.53, 1.54,
1.56, 1.567, 1.578, 1.594, 1.602, 1.604, 1.608, 1.63, 1.642,
1.674, 1.73, 1.75, 1.75, 1.763, 1.768, 1.781, 1.782, 1.792, 1.82,
1.868, 1.881, 1.89, 1.893, 1.895, 1.91, 1.923, 1.94, 1.945, 2.023,
2.1, 2.13, 2.215, 2.268, 2.44)

bs.plot(data)
bs.plot(data, main="BS probability plot", lty.line=2, pch=3, col.line="red")

ticklabels=c(0.01, seq(0.1,0.9,by=0.1), seq(0.91,0.99,by=0.01) )
qq = quantile(data, probs=ticklabels) 
ticksat= qnorm(ticklabels)* sqrt( qq )
hline =  qnorm( ticklabels  ) *  sqrt( qq )

abline( h=hline,    col=gray(0.5), lty=3, lwd=0.6 )
abline( v=  seq(0, 2.5, by=0.1), col=gray(0.5), lty=3, lwd=0.5 )

abline( h= qnorm(0.5)*sqrt(median(data)), col=gray(0.1), lty=1, lwd=0.6 )
abline( v= median(data), col=gray(0.1), lty=1, lwd=0.6 )
}


