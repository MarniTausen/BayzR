#' Summary statistics of a bayz object
#'
#' Produced summary statistics of the bayz object and returns a summarybayz object.
#'
#' @param object    bayz object
#' @param ...       additional parameters
#'
#' @return summarybayz object
#' @import stats
#' @export
summary.bayz <- function(object, ...){

    new_object <- list()
    class(new_object) <- "summarybayz"
    new_object[['terms']] <- object[['terms']]
    new_object[['nError']] <- object[['nError']]
    new_object[['Errors']] <- object[['Errors']]
    if(new_object$nError>0){
        return(new_object)
    }

    new_object[['Estimates']] <- object[['Estimates']]
    new_object[['Samples']] <- object[['Samples']]
    # Collect parameter information
    par <- object$Parameters
    # Get variables of interrest
    int_vars <- rownames(par)[par$Hyper]

    variables <- c()
    heritabilities <- c()

    total_variance <- sum(object$Estimates[int_vars,"postMean"])
    for(var in int_vars){
        variables <- c(variables, rev(unlist(strsplit(var, split=".", fixed=TRUE)))[1])
        heritabilities <- c(heritabilities, (object$Estimates[var,"postMean"]/total_variance))
    }

    heritability_table <- data.frame("Variable"=variables, "Heritability"=heritabilities)
    new_object[["Heritability"]] <- heritability_table
    new_object[["Random_variables"]] <- int_vars

    return(new_object)

}
