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
#' @import stats
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

    # Main estimates table
    par <- object$Parameters
    par_select = ( par$Hyper==FALSE & par$ModelTerm!="rp" & par$Size <= maxLevel & substr(rownames(par),0,3)!="var")
    par_not_shown = ( par$Hyper==FALSE & par$ModelTerm!="rp" & par$Size > maxLevel & substr(rownames(par),0,3)!="var")
    estim_select = c()
    for(i in 1:nrow(par)) {
        if(par_select[i]) estim_select = c(estim_select,par$EstStart[i]:par$EstEnd[i])
    }
    new_object[['Estimates']] = object$Estimates[estim_select,]
    new_object[['Estimates_notShown']] = rownames(par)[par_not_shown]


    # 'logged' parameters analysis with convergence diagnostics and HPD intervals
    output_cycles = as.numeric(rownames(object$Samples))
    if(length(output_cycles<10)) convergence_status = 1       # fail because too few output
    convergence_table = data.frame()
    if(require("coda")) {
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
    }
    else {   # no coda package available
        convergence_status = 2                                # fail because of missing coda package
    }
    new_object[['Convergence']] = convergence_table
    new_object[['ConvergenceStatus']] = convergence_status

    # Proportions of variances.
    varnames = ( substr(colnames(object$Samples),0,3)=="var" )
    variance_table = data.frame()
    if (sum(varnames)>1) {   # if <= 1 variances the variance_table remains an empty data frame
       var_samples = as.matrix(object$Samples[,varnames])
       var_proportions = t(apply(var_samples,1,function(x){x/sum(x)}))
       postMeans =  apply(var_proportions,2,mean)
       postSDs  = apply(var_proportions,2,sd)
       HPDs = HPDbayz(var_proportions,prob=HPDprob, bound="prob")
       variance_table = data.frame(postMeans,postSDs,HPDs)
       colnames(variance_table) = c("postMean","postSD","HPDleft","HPDright")
    }
    new_object[['variance_table']] = variance_table

    new_object[['HPDprob']] = HPDprob

    return(new_object)

}
