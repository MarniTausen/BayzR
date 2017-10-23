#' Print function for bayz class
#'
#' Prints the bayz class
#'
#' @param x    bayz object
#'
#' @return Nothing, only prints
#' @import stats
#' @export
print.bayz <- function(x){
    cat(x$modelname, "\n\n")
    cat("    model:", deparse(x$terms), "\n\n")
    cat("Model function:", x$modelfunction,"\n")
}
