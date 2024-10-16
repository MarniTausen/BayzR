#' Summary of bayz model fit
#'
#' Produced summary statistics of the bayz object and returns a summarybayz object.
#'
#' @param object    bayz output object
#' @param maxLevel  controls the inclusion of estimates in the "Main estimates" table by selecting
#'                  only the model-terms with less <= maxlevel coefficients
#' @param HPDprob   the probability for the Highest Posterior Density intervals reported in the
#'                  table "Convergence diagnostics"
#' @param ...       additional parameters
#'
#' @return summarybayz object
#' @import stats coda
#' @export
summary.bayz <- function(object, maxLevel=10, HPDprob=0.95, ...){

    new_object <- list()
    class(new_object) <- "summarybayz"
    new_object[['terms']] <- object[['terms']]
    new_object[['nError']] <- object[['nError']]
    new_object[['Errors']] <- object[['Errors']]
    new_object[['Chain']] <- object[['Chain']]
    if(new_object$nError>0){
        return(new_object)
    }

    # Main coefficient estimates
    par <- object$Parameters
    par_select = ( par$Hyper==FALSE & par$ModelTerm!="rp" & par$Size <= maxLevel & substr(rownames(par),0,3)!="var")
    par_not_shown = ( par$Hyper==FALSE & par$ModelTerm!="rp" & par$Size > maxLevel & substr(rownames(par),0,3)!="var")
    par_select_names = object$Parameters[par_select, "Param"]
    coeff_select = c()
#    for(i in 1:nrow(par)) {
#        if(par_select[i]) estim_select = c(estim_select,par$EstStart[i]:par$EstEnd[i])
#    }
#    new_object[['Estimates']] = object$Estimates[estim_select,]
#    new_object[['Estimates_notShown']] = rownames(par)[par_not_shown]


    # 'logged' parameters analysis with convergence diagnostics and HPD intervals
    output_cycles = as.numeric(rownames(object$Samples))
    if(length(output_cycles)<10) convergence_status = 1       # fail because too few output
    convergence_table = data.frame()
    samp = coda::mcmc(object$Samples, start=output_cycles[1], end=output_cycles[length(output_cycles)],
                 thin=output_cycles[2]-output_cycles[1])
    effSizes = coda::effectiveSize(samp)
    postMeans =  apply(object$Samples,2,mean)
    postSDs  = apply(object$Samples,2,sd)
    MCSEs = postSDs / sqrt(effSizes)
    MCCVpct = 100 * MCSEs / abs(postMeans)
    HPDbounds = rep("none",ncol(object$Samples))
    HPDbounds[substr(colnames(object$Samples),0,3)=="var"] = "var"
    HPDs = HPDbayz(object$Samples,prob=HPDprob, bound=HPDbounds)
    GewekeZ = abs(coda::geweke.diag(object$Samples)$z)
    convergence_table = data.frame(postMeans,postSDs,effSizes,GewekeZ,MCSEs,MCCVpct,HPDs)
    colnames(convergence_table) = c("postMean","postSD","effSize","GewekeZ","MCSE","MCCV%","HPDleft","HPDright")
    rownames(convergence_table) = colnames(object$Samples)
    convergence_status = 0                                # success
    new_object[['Convergence']] = convergence_table
    new_object[['ConvergenceStatus']] = convergence_status

    new_object[['HPDprob']] = HPDprob

    return(new_object)

}
