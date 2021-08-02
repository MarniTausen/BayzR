#' Extract posterior means and SDs for `fixed´ effects
#'
#' Extract from the Bayz output the coefficients (posterior mean and SD) for terms fitted with fx() and rg().
#' With the default splitLabels=TRUE, labels for interactions such as a%b will be split in multiple columns.
#'
#' @param object        A bayz model output
#' @param splitLabels   Whether labels that contain % (in interactions) such as a%b should be split
#'                      in multiple columns (default TRUE). 
#' @param ...           Additional parameters passed onto the Model function.
#'
#' @return a list with one member (a data frame) for each fixed effect
##' @method fixef bayz
##' @import nlme
#' @export fixef.bayz
#' @export
fixef.bayz <- function(object, splitLabels=TRUE, ...){
    par = object$Parameters
    par_select = ( par$ModelTerm=="fx" | par$ModelTerm=="rg" )
    return(coef.bayz(object, par_select, splitLabels))
}
