#' Extract posterior means and SDs for all or a selection of coefficients
#'
#' Extract from the Bayz output all or a selection of coefficients (posterior mean and SD), i.e. estimates from the
#' fx(), rg() and rn() terms.
#'
#' @param object        A bayz model output
#' @param which         A vector TRUE/FALSE to indicate which model coefficients to extract (should match the rows
#'                      in the $Parameters table in a bayz output object). It can be omitted to extract all coefficients.
#' @param splitLabels   Whether labels that contain % (in interactions) such as a%b should be split
#'                      in multiple columns (default TRUE). 
#' @param ...           Additional parameters.
#'
#' @return a list with one member (a data frame) for each coefficient
#' @import stats
#' @export
coef.bayz <- function(object, which=NULL, splitLabels=TRUE, ...){
    par = object$Parameters
    est = object$Estimates
    if(is.null(which)) {
       which = ( par$ModelTerm=="fx" | par$ModelTerm=="rg" | par$ModelTerm=="rn" | par$ModelTerm=="rr" )
    }
    else {
       if(length(which) != nrow(par)) stop("In coef.R which= is not of right length")
       # can also check that which vector is proper TRUE/FALSE
    }
    par = par[which,]
    new_object = list()
    for(i in 1:nrow(par)) {
        estim_select = est[par$EstStart[i]:par$EstEnd[i],]
        if (!grepl("%",rownames(estim_select)[1],fixed=TRUE)) splitLabels = FALSE
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
