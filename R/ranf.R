#' Specify fit of a random factor, optionally with a kernel/similarity matrix
#'
#' @param x     Variable (forced to be factor) to be used as random effect in the model.
#' @param V     Kernel/similarity/var-covar matrix to be included so that x~N(0,Vsigma^2). Uses eigen decomposition.
#' @param rrankpct Reduced Rank threshold for including eigenvectors, given in percent (default 99)
#' @param prior Prior information to be added.
#'
#' @return Object storing the factor information.
#' @export

ranf <- function(x, V=NULL, rrankpct=99, prior=NULL) {
    x <- as.factor(x)
    if(!is.null(V)) {
        if(is.null(rownames(V))) {
            stop("There are no rownames on V")
        }
        EVdecomp <- eigen(V)
        colnames(EVdecomp$vectors) <- paste("evec",1:ncol(EVdecomp$vectors),sep="")
        rownames(EVdecomp$vectors) <- rownames(V)
        attr(x, "evalues") <- EVdecomp$values
        attr(x, "evectors") <- EVdecomp$vectors
    }
    if (!is.null(prior)) attr(x,"prior") <- prior
    attr(x,"rrankpct") <- rrankpct
    return(x)
}
