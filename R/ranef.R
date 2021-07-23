#' Extract posterior means and SDs for `random´ effects
#'
#' Extract from the Bayz output the estimates (posterior mean and SD) for those terms that would be
#' called random in mixed models,
#' i.e. those that have shrinkage priors such as a Normal prior, but also other shrinkage priors.
#'
#' @param object  A bayz model output
#' @splitLabels   Whether labels that contain % (in interactions) such as a%b should be split in multiple columns (default TRUE). 
#' @param ...     Additional parameters passed onto the Model function.
#'
#' @return 
#' @export
ranef.bayz <- function(object, splitLabels=TRUE, ...){
}
