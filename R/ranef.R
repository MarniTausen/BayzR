#' Extract posterior means and SDs for `random´ effects
#'
#' Extract from the Bayz output the coefficients (posterior mean and SD) for rn() model terms.
#'
#' @param object        A bayz model output
#' @param splitLabels   Whether labels that contain % (in interactions) such as a%b should be split
#'                      in multiple columns (default TRUE). 
#' @param ...           Additional parameters passed onto the Model function.
#'
#' @return a list with one member (a data frame) for each random effect
#' @importFrom lme4 ranef
#' @export
ranef.bayz <- function(object, splitLabels=TRUE, ...){
    par = object$Parameters
    par_select = ( par$ModelTerm=="rn" )
    return(coef.bayz(object, par_select, splitLabels))
}
