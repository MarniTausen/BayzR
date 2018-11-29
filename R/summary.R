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
    cat("Summary statistics of Bayz object\n\n")
    cat("    model formula:", deparse(object$terms), "\n\n")
    if(object$nError>0){
        cat("Bayz encountered errors while running:\n")
        for(errormsg in object$Errors){
            cat("\t",errormsg,"\n")
        }
    } else {
        # Collect parameter information
        par <- object$Parameters
        # Get variables of interrest
        int_vars <- rownames(par)[par$Hyper]

        cat("Random variable estimates:\n")
        print(object$Estimates[int_vars,])

        cat("\n")

        cat("Broad sense heritability:\n")
        total_variance <- sum(object$Estimates[int_vars,"postMean"])
        for(var in int_vars){
            cat(paste(rev(unlist(strsplit(var, split=".", fixed=TRUE)))[1], ":", sep=""),
                "\t\t",(object$Estimates[var,"postMean"]/total_variance)*100,"%\n")
        }
    }
}
