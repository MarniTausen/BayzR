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
summary.bayz <- function(x, ...){
    cat(x$modelname, "\n\n")
    cat("    model formula:", deparse(x$terms), "\n\n")
    cat("Model function:", x$modelfunction,"\n\n")

    cat("Estimates:\n")
    print(x$Estimates)

    cat("\n")

    cat("Heritability:\n")
    print(x$Estimates)
}
