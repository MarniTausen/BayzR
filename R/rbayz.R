#' Bayz test function
#'
#' Fits a model, given as a formula, optionally with data provided through the "..." parameter.
#'
#' @param model   A formula describing the model.
#' @param data    Data frame to collect data from
#' @param chain   Number of iteration/sample to generate.
#'
#' @return A fitted bayz model
#' @import stats
#' @export
#'
#' @useDynLib BayzR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
rbayz <- function(model, data=NULL, chain=NULL){

    modeldata <- model.frame(model, data=data)
    if (is.null(chain)){
        chain=c(1100,100,10)
        cat("Warning: running the default chain of 1100 cycles, this may be too short for many analyses\n")
    }
    chain <- as.integer(chain)

    return(rbayz_cpp(modeldata, chain))
}
