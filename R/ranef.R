#' Extract posterior means and SDs for `random´ effects
#'
#' Extract from the Bayz output the estimates (posterior mean and SD) for those terms that would be
#' called random in mixed models,
#' i.e. those that have shrinkage priors such as a Normal prior, but also other shrinkage priors.
#'
#' @param object        A bayz model output
#' @param splitLabels   Whether labels that contain % (in interactions) such as a%b should be split
#'                      in multiple columns (default TRUE). 
#' @param ...           Additional parameters passed onto the Model function.
#'
#' @return 
#' @export
ranef.bayz <- function(object, splitLabels=TRUE, ...){
    new_object = list()
    par = object$Parameters
    est = object$Estimates
    par_select = ( par$ModelTerm=="fx" | par$ModelTerm=="rn" )
    par = par[par_select,]
    for(i in 1:nrow(par)) {
        estim_select = est[par$EstStart[i]:par$EstEnd[i],]
        if(splitLabels) {
            splitlabels = t(as.data.frame(strsplit(rownames(estim_select),"%")))
            colnm = splitlabels[1,1]
            splitlabels = splitlabels[,-1]
            estim_table = data.frame(splitlabels, estim_select)
            colnames(estim_table) = c(unlist(strsplit(colnm,":")),colnames(estim_select))
            rownames(estim_table) = rownames(estim_select)
        }
        else {
            estim_table = estim_select
        }
        new_object[[rownames(par)[i]]] = estim_table
    }
    return(new_object)
}
