#' mfactor function to build multi-column/matrix factor from list of single factors
#'
#' @param ... list of factors to be merge into one multi-column factor
#'
#' @return object of type mfactor
#' @export

mfactor <- function(...) {
   fact_list <- list(...)
   return(mfactor_cpp(fact_list))
}
