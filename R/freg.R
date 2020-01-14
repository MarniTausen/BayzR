#' Fixed-effect regression on continuous covariable
#'
#' @param x Vector of covariates (will be forced to be numeric) to fit in the model as a fixed regression.
#'
#' @return Numeric version of the vector
#' @export

freg <- function(x) as.numeric(x)
