#' Bayesian mixed model with genotypes
#'
#' Bayesian mixed model with table of genotypes for random regression (rrBLUP/SNP-BLUP) fit
#'
#' @param data    Data frame to collect data from
#' @param geno    Table with genotypes
#' @param ...     Additional parameters passed onto the Model function (for example )
#'
#' @return Fitted mixed model
#' @import stats
#' @export
mmg <- function(data=NULL, geno=NULL, ...){
    if(is.null(data)) return("Bayesian mixed model with genotypes")
    if(is.null(geno)){
        stop("No genotype table has been provided! \n  Please use geno=..., to provide a genotype matrix")
    }
    data
}
