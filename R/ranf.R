#' Set Random variable
#'
#'  Sets a variable to a random variable in a linear mixed model.
#'
#' @param x     Variable vector which should be set to random in a model.
#' @param V     Similarity/var-covar matrix to be included so that x~N(0,Vsigma^2). Uses eigen decomposition
#' @param rrankpct Reduced Rank threshold for including eigenvectors in percent (default 99\%)
#' @param prior Prior information to be added.
#'
#' @return Random variable information.
#' @export

ranf <- function(x, V=NULL, rrankpct=99, prior=NULL) {
    x <- as.factor(x)
    if(!is.null(V)) {
        EVdecomp <- eigen(V)
        colnames(EVdecomp$vectors) <- paste("evec",1:ncol(EVdecomp$vectors),sep="")
        attr(x, "evalues") <- EVdecomp$values
        attr(x, "evectors") <- EVdecomp$vectors
    }
    if (!is.null(prior)) attr(x,"prior") <- prior
    attr(x,"rrankpct") <- rrankpct
    return(x)
}
