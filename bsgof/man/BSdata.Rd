\name{BSdata}
\alias{BSdata}
\alias{psi21k}
\alias{psi26k}
\alias{psi31k}
\alias{bearings}
\alias{fatigue}
\alias{repair}

\title{Data set}

\description{
\code{psi21k}, \code{psi26k}, and \code{psi31k} are from  Birnbaum and Saunders (1969). 
The fatigue lifetimes of aluminum specimens exposed to a maximum stress of 21,000 psi, 26,000 psi, 31,000 psi,
respectively. 

\code{bearings} is from McCool (1974). The fatigue lifetimes (in hours) of ten bearings.

\code{fatigue} is from Brown and Miller (1978). The fatigue lifetimes of cylindrical specimens 
subjected to combined torsional and axial loads over constant-amplitude cycles until failure.

\code{repair} is from Hsieh (1990). This is a maintenance data set 
on active repair times (in hours) for an airborne communications transceiver.
}

\usage{
BSdata
}

\references{
Birnbaum, Z. W. and Saunders, S. C. (1969). A new family of life distributions. J. Appl. Probab. 6(2): 637-652. \cr 

McCool, J. I. (1974). Inferential techniques for Weibull populations. Aerospace Research Laboratories Report ARL T
R74-0180, Wright-Patterson Air Force Base, Dayton, OH.  \cr 

Rieck, J. R. and Nedelman, J. (1991). 
A Log-Linear Model for the Birnbaum-Saunders Distribution. Technometrics. 33, 51-60.  \cr 

Brown, M. W. and Miller, K. J. (1978). 
Biaxial Fatigue Data. Report CEMR1/78. University of Sheffield, Dept. of Mechanical Engineering.  \cr 

Hsieh, H. K. (1990). 
Estimating the Critical Time of Inverse Gaussian Hazard Rate. IEEE Transactions on Reliability, 39(10): 342-345.
}


\examples{
# Attach data sets
attach(BSdata)
psi21k
psi26k
psi31k
bearings
fatigue
repair
}


\keyword{datasets}
