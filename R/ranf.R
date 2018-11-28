#' Set Random variable
#'
#'  Sets a variable to a random variable in a linear mixed model.
#'
#' @param x     Variable vector which should be set to random in a model.
#' @param cor   Correlation matrix to be included. Default uses eigen decomposition
#' @param prior Prior information to be added.
#'
#' @return Random variable information.
#' @export

ranf <- function(x, cor=NULL, prior=NULL) {
    x <- as.factor(x)
    if(!is.null(cor)) {
        EVdecomp <- eigen(cor)
        attr(x, "evalues") <- EVdecomp$values
        attr(x, "evectors") <- EVdecome$vectors
    }
    if (!is.null(prior)) attr(x,"prior") <- prior
    return(x)
}
