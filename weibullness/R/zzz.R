
#.onAttach <- function(...){ cat("\n weibullness Package is installed. \n\n") }
.onAttach <- function(...){
  packageStartupMessage("\n weibullness Package is installed. \n\n")
}
