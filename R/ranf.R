#' Set Random variable
#'
#'  Sets a variable to a random variable in a linear mixed model.
#'
#' @param x     Variable vector which should be set to random in a model.
#' @param V   Correlation matrix to be included. Default uses eigen decomposition
#' @param prior Prior information to be added.
#'
#' @return Random variable information.
#' @export

ranf <- function(x, V=NULL, prior=NULL) {
    x <- as.factor(x)
    if(!is.null(V)) {
        EVdecomp <- eigen(V)
        colnames(EVdecomp$vectors) <- paste("evec",1:ncol(EVdecomp$vectors),sep="")
        attr(x, "evalues") <- EVdecomp$values
        attr(x, "evectors") <- EVdecomp$vectors
    }
    if (!is.null(prior)) attr(x,"prior") <- prior
    return(x)
}
