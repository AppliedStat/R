\encoding{utf-8}
\name{plot}
\alias{plot.rcc}
\alias{plot.acc}

\title{Plot quality control chart}

\description{Plot method for an object of class \code{"rcc"} or \code{"acc"}.}

\usage{
\method{plot}{acc}(x, digits= getOption("digits")-2, col.line="cyan4", 
                   col.background="ivory", col.boundary="ivory4",
                   cex.text=0.7, x.text, text.offset=0, LCL, CL, UCL, \dots)

\method{plot}{rcc}(x, digits= getOption("digits")-2, col.line="cyan4", 
                   col.background="ivory", col.boundary="ivory4", 
                   cex.text=0.7, x.text, text.offset=0, LCL, CL, UCL, \dots)
}

\arguments{
 \item{x}{an object of \link{class} \code{"rcc"}.}
 \item{digits}{number of significant digits to use, see \link{print}.}
 \item{col.line}{line color connecting subgroup estimates.}
 \item{col.background}{background color for in-control region.}
 \item{col.boundary}{boundary color for in-control region.}
 \item{cex.text}{magnification to be used for the text labels (LCL, CL, UCL).}
 \item{x.text}{a vector of x-coordinates where the text labels (LCL, CL, UCL) should be written.
          If missing, the text labels will be written on the right side.}
 \item{text.offset}{an offset value for the text labes (LCL, CL, UCL).}
 \item{LCL}{lower control limit. If missing, use the value from an object of class \code{"rcc"} or \code{"acc"}.} 
 \item{CL}{center line.  If missing, use the value from an object of class \code{"rcc"} or \code{"acc"}.}
 \item{UCL}{upper control limit. If missing, use the value from an object of class \code{"rcc"} or \code{"acc"}.}
 \item{\dots}{additional parameters to plot.}
}

\seealso{
\code{rcc}\{rQCC\}, \code{rcc}\{rQCC\},  \code{acc}\{rQCC\},  \link{plot}
}

\author{Chanseok Park}
