\name{news}
\title{Package news for rQCC package}


\section{version 2.22.5}{
  \itemize{
    \item Some references are updated.
    \item new \code{racc} function for the robust attributes control charts is added.
    \item References are updated. 
}}

\section{version 2.22.1}{
  \itemize{
    \item The argument '\code{method}' is changed to '\code{estimator}.'
    \item The new \code{pooledEstimator} function is added. 
          This function calculates the pooled estimator based on the unbiased estimators such as the mean, median, Hodges-Lehmann (HL1, HL2, HL3), standard deviation, range, median absolute deviation (MAD) and Shamos estimators..
    \item The new \code{acc} function for the attributes control charts is added. 
    \item The new \code{sd.unbiased} function is added.
    \item The updated \code{rcc}  function can handle unbalanced samples.
    \item The updated \code{evar} function can handle multiple samples.
    \item The updated \code{RE}   function can handle multiple samples.
    \item The \code{factors.cc} is improved for the case of '\code{factor="d3"}.'
}}

\section{version 1.20.7}{
  \itemize{
    \item The S and R charts are added in the \code{rcc} function.
    \item The \code{vignettes} for the \code{factors.cc} and \code{rcc} functions are added.
    \item The \code{A3.factor} function is replaced with the updated \code{factors.cc} function.
    \item Typos are corrected. For example, some typos in the help page of shamos function. 
}}

\section{version 1.20.3}{
  \itemize{
    \item The graphics parts were added.
    \item The \code{summary} function for the output from \code{rcc} was added.
    \item The \code{print}   function for the output from \code{rcc} was added.
    \item The \code{plot}    function for the output from \code{rcc} was added.
}}

\section{version 0.19.8.2}{
  \itemize{
    \item The first \code{rQCC} R package was developed.
    \item Unbiased robust estimates were implemented.
    \item Incorporated into control charting applications.
}}

