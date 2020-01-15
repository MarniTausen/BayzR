#' Bayz function
#'
#' Fits a model, given as a formula, optionally with data provided through the "..." parameter.
#'
#' @param model   A formula describing the model.
#' @param data    Data frame to collect data from
#' @param fct     Model function
#' @param chain   Vector describing the number of iterations to be run.
#' @param ...     Additional parameters passed onto the Model function.
#'
#' @return A fitted bayz model
#' @import stats
#' @export
#'
#' @useDynLib BayzR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
bayz <- function(model, data=NULL, fct=NULL, chain=NULL, ...){
    bayz_check_modelterms(model)
    model_data <- model.frame(model, data=data)
    if (is.null(chain)){
        chain=c(1100,100,10)
        warning("running the default chain of 1100 cycles, this may be too short for many analyses\n")
    }
    chain <- as.integer(chain)
    if(is.function(fct)){
        result <- fct(data=model_data, model=model, chain=chain, ...)
    } else {
        result <- rbayz_cpp(model_data, chain)
    }
    class(result) <- "bayz"
    #result[['modelname']] <- fct()
    #result[['modelfunction']] <- deparse(substitute(fct))
    result[['terms']] <- terms(model)
    return(result)
}
