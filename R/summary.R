#' Summary bayz model fit
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
    narrow_sense <- c()

    ## Groups the variables by model type and variable name
    items <- .group.variable.names(int_vars)

    ## Calculate Heritability based on each model type and variable name
    ## If there are multiple estimate a broadsense heritability
    ## Narrowsense is only calculated if ranf has a V or if ran2f has V1 and/or V2
    total_variance <- sum(object$Estimates[int_vars,"postMean"])
    for(item in items){
        if(length(item$var.names)>1){
            broadsense <- sum(object$Estimates[item$var.index,"postMean"])/total_variance
            heritabilities <- c(heritabilities, broadsense)
            if(item$model.type=="ranf"){
                variables <- c(variables, unlist(strsplit(item$var.names[1], split=".", fixed=TRUE))[1])
                n_elements <- sapply(item$var.names, function(x) length(unlist(strsplit(x, split=".", fixed=TRUE))))
                narrow_sense <- c(narrow_sense, object$Estimates[item$var.index[which(n_elements==2)], "postMean"]/total_variance)
            } else if(item$model.type=="ran2f"){
                current_name <- paste(unlist(strsplit(item$var.names[1], split=".", fixed=TRUE))[1:2], collapse=" x ")
                variables <- c(variables, current_name)
                n_elements <- sapply(item$var.names, function(x) length(unlist(strsplit(x, split=".", fixed=TRUE))))
                narrow_sense <- c(narrow_sense, object$Estimates[item$var.index[which(n_elements>2)], "postMean"]/total_variance)
            } else {
                variables <- c(variables, unlist(strsplit(item$var.names[1], split=".", fixed=TRUE))[1])
                narrow_sense <- c(narrow_sense, NA)
            }
        } else {
            variables <- c(variables, item$var.names)
            heritabilities <- c(heritabilities, (object$Estimates[item$var.index,"postMean"])/total_variance)
            narrow_sense <- c(narrow_sense, NA)
        }
    }

    if(all(is.na(narrow_sense))){
        heritability_table <- data.frame("Variable"=variables, "Heritability"=heritabilities, stringAsFactors=FALSE)
        heritability_table$stringAsFactors <- NULL
        new_object[["Heritability"]] <- heritability_table
        new_object[["Random_variables"]] <- int_vars
    } else {
        heritability_table <- data.frame("Variable"=variables, "Broad.sense.Heritability"=heritabilities,
                                         "Narrow.sense.Heritability"=narrow_sense)
        new_object[["Heritability"]] <- heritability_table
        new_object[["Random_variables"]] <- int_vars
    }

    return(new_object)

}

.group.variable.names <- function(int_vars){
    targets <- c("ranf", "ran2f")
    n <- length(int_vars)
    items <- list()
    ## ASSUME resid is always first
    items[["resid"]] <- list(var.names="Residual (resid)", var.index=int_vars[1],
                             model.type="resid")
    for(var in int_vars[2:n]){
        labels <- unlist(strsplit(var, split=".", fixed=TRUE))
        m <- length(labels)
        element <- paste(labels[2:3], collapse=".")
        name <- paste(labels[3:m], collapse=".")
        if(element %in% names(items)){
            items[[element]]$var.names <- c(items[[element]]$var.names, name)
            items[[element]]$var.index <- c(items[[element]]$var.index, var)
        } else {
            items[[element]] <- list(var.names=name, var.index=var,
                                     model.type=labels[2])
        }
    }
    items
}
