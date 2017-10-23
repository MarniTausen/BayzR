#' Bayz function
#'
#' Fits a model, given as a formula, optionally with data provided through the "..." parameter.
#'
#' @param model   A formula describing the model.
#' @param data    Data frame to collect data from
#' @param fct     Model function
#' @param ...     Additional parameters passed onto the Model function (for example )
#'
#' @return A fitted bayz model
#' @import stats
#' @export
bayz <- function(model, data=NULL, fct=NULL, ...){
    model_data <- model.frame(model, data=data)
    if(is.function(fct)){
        result <- fct(data=model_data, model=model, ...)
    } else {
        stop("No model function provided")
    }
    class(result) <- "bayz"
    result[['modelname']] <- fct()
    result[['modelfunction']] <- deparse(substitute(fct))
    result[['terms']] <- terms(model)
    return(result)
}

fixed <- function(x) x
random <- function(x) x
