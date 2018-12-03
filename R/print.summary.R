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
print.summarybayz <- function(x, ...){
    object <- x
    cat("Summary statistics of Bayz object\n\n")
    cat("   model formula:", deparse(object$terms), "\n\n")
    if(object$nError>0){
        cat("Bayz encountered errors while running:\n")
        for(errormsg in object$Errors){
            cat("  ",errormsg,"\n")
        }
    } else {
        int_vars <- object[["Random_variables"]]

        cat("Random variable estimates:\n")
        print(object$Estimates[int_vars,])

        cat("\n")

        cat("Variance explained (Heritability):\n")
        print(object$Heritability)
    }
}
