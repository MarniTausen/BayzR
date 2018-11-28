#' Bayesian mixed model
#'
#' Fits a model, given as a formula, optionally with data provided through the "..." parameter.
#'
#' @param data    Data frame to collect data from
#' @param ...     Additional parameters passed onto the Model function (for example )
#'
#' @return Fitted mixed model
#' @import stats
#' @export
mm <- function(data=NULL, chain=NULL, ...){
    if(is.null(data)) return("Bayesian mixed model")
    return(rbayz_cpp(data, chain))
}
