#' Print function for bayz output class
#'
#' This method is usually not seen or used directly by the user, but allows to see a compact version of the output
#' in the console without overflowing the screen, for instance when running bayz() without catching the output, or
#' to quickly see error messages without the need to use summary().
#'
#' @param x    bayz object
#' @param ...  additional parameters
#'
#' @return Nothing, only prints
#' @import stats
#' @export
print.bayz <- function(x, ...){
    cat("Bayz output object - Bayesian mixed and shrinkage models", "\n\n")
#    cat("    model:", deparse(x$terms), "\n\n")
    
    if(x$nError>0){
        cat("Bayz encountered errors while running:\n")
        for(errormsg in x$Errors){
            cat("\t",errormsg,"\n")
        }
    }
}
