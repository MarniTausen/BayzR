#' Print function for bayz class
#'
#' Prints the bayz class
#'
#' @param x    bayz object
#' @param ...  additional parameters
#'
#' @return Nothing, only prints
#' @import stats
#' @export
print.bayz <- function(x, ...){
    cat("Bayz object - Bayesian mixed model", "\n\n")
    cat("    model:", deparse(x$terms), "\n\n")
    if(x$nError>0){
        cat("Bayz encountered errors while running:\n")
        for(errormsg in x$Errors){
            cat("\t",errormsg,"\n")
        }
    }
}
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
