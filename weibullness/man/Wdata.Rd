\name{Wdata}
\alias{Wdata}
\alias{bearings}
\alias{glassfiber1.5}
\alias{glassfiber15}
\alias{urinary}
\alias{radiotherapy}
\alias{radio.chemotherapy}

\title{Data set}

\description{
\code{bearings}: It is from Lieblein and Zelen (1956). 
These data are deep-groove ball bearings failure times (number of millions of revolutions) in eudurance tests.

\code{glassfiber1.5} and \code{glassfiber15}: They are from Smith and Naylor (1987).
These data sets are from experimental data for the strength of glass fiber of length 1.5cm and 15cm, respectively.

\code{urinary}: It is from Santiago and Smith (2013).
It is about the days in between discharge of males in nosocomial urinary tract infections in patients.

\code{radiotherapy} and \code{radio.chemotherapy}: They are from Finkelstein, D. M. (1986) 
and Lindsey, J. C. and L. M. Ryan (1998). 
These data are interval-censored observations from a study of patients with breast cancer. 
The measurement is the time to cosmetic deterioration of the breast for women who received radiotherapy 
and women who received radio-chemotherapy.
}

\usage{
Wdata
}

\references{
Lieblein, J. and M. Zelen (1956). 
Statistical Investigation of the Fatigue Life of Deep-Groove Ball Bearings.
Journal of Research of the National Bureau of Standards, 57 (5), 273-316. 

Smith, R. L. and J. C. Naylor (1987). 
A comparison of maximum likelihood and Bayesian estimators for the three-parameter Weibull distribution. 
Journal of the Royal Statistical Society. Series C (Applied Statistics), 36 (3), 358-369.  

Santiago, E. and J. Smith (2013). Control Charts Based on the Exponential Distribution: 
Adapting Runs Rules for the t Chart, Quality Engineering, 25 (2), 85-96.  

Finkelstein, D. M. (1986), A proportional hazards model for interval-censored failure time data. Biometrics, 42, 845-865.

Lindsey, J. C. and L. M. Ryan (1998). Tutorial in biostatistics: Methods for interval-censored data. Stat. Med., 17, 219-238.
}


\examples{
# Attach data sets
attach(Wdata)
bearings
glassfiber1.5
glassfiber15
urinary
radiotherapy
radio.chemotherapy
}


\keyword{datasets}
