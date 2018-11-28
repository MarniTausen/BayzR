#' Print function for bayz class
#'
#' Prints the bayz class
#'
#' @param object    bayz object
#' @param ...  additional parameters
#'
#' @return Nothing, only prints
#' @import stats
#' @export
summary.bayz <- function(object, ...){
    cat(object$modelname, "\n\n")
    cat("    model formula:", deparse(object$terms), "\n\n")
    cat("Model function:", object$modelfunction,"\n\n")

    cat("Estimates:\n")
    print(object$Estimates)

    cat("\n")

    cat("Heritability:\n")
    print(object$Estimates)
}
