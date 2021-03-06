\encoding{utf-8}
\name{plot.rcc}
\alias{plot.rcc}

\title{Plot robust quality control chart}

\description{Plot method for an object of class \code{"rcc"}.}

\usage{
\method{plot}{rcc}(x, digits=getOption("digits")-2, cex.text=0.7, x.text=1, \dots)
}

\arguments{
 \item{x}{an object of \link{class} \code{"rcc"}.}
 \item{digits}{number of significant digits to use, see \link{print}.}
 \item{cex.text}{magnification to be used for the text labels (LCL, CL, UCL).}
 \item{x.text}{x-coordinate where the text labels (LCL, CL, UCL) should be written.}
 \item{\dots}{additional parameters to plot.}
}

\seealso{
\code{rcc}\{rQCC\}, \link{plot}
}

\author{Chanseok Park}
