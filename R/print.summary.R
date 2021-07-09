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
    cat("Summary of bayz model fit\n\n")
    cat("model formula:", deparse(object$terms), "\n\n")
    if(object$nError>0){
        cat("Bayz encountered errors while running:\n")
        for(errormsg in object$Errors){
            cat("  ",errormsg,"\n")
        }
    } else {
#       No errors

        if (length(x$Estimates_notShown)==0) cat("Estimates for model coefficients:\n")
        else cat("Estimates for model coefficients(*):\n")
        print(object$Estimates)
        if (length(x$Estimates_notShown)>0) {
            cat("* Some estimates are not shown because they have more than maxLevel levels:\n")
            cat(paste(" ",x$Estimates_notShown))
            cat("\n")
        }
        cat("\n")

        cat("Estimates and HPD intervals for hyper-parameters and other 'logged' parameters:\n")
        if(x$ConvergenceStatus == 1) {
            cat("*** This table is not printed because there are fewer than 10 output samples ***\n")
        }
        else if (x$ConvergenceStatus == 0) {
            print(object$Convergence[,c("postMean","postSD","HPDleft","HPDright")])
        }
        cat("\n")

        cat("Convergence diagnostics on 'logged' parameters:\n")
        if (x$ConvergenceStatus == 2) {
            cat("*** This table is not printed because the coda package is not installed ***\n")
        }
        else if (x$ConvergenceStatus == 0) {
            print(object$Convergence[,c("effSize","GewekeZ","MCSE","MCCV%")])
        }
        cat("\n")

        cat("Estimates for explained variances (variances as proportions of total):\n")
        if (nrow(object$variance_table) > 0) print(object$variance_table)
        else cat("*** This table is not printed because there is only one variance in the model\n")
    }
}
