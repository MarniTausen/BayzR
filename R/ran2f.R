#' Fit two random factors with interactions, and optional attach covariance structures to first or both variables.
#'
#' With two covariance structures specified it will build the interaction-relationship matrix.
#'
#' @param x1        First variable vector (a factor) to build two-factor interaction
#' @param x2        Second variable vector (a factor) to build two-factor interaction
#' @param V1        Similarity/var-covar matrix to be included on first variable. Uses/stores the eigen decomposition of V1.
#' @param V2        Similarity/var-covar matrix to be included on second variable. Uses/stores the eigen decomposition of V2.
#' @param rrankpct  Reduced Rank threshold for including eigenvectors in percent (default 99\%)
#' @param prior     Prior information for the variance component.
#'
#' @return The first factor data is returned, the second factor and V1 and V2 eigenvectors (when available) are stored as attributes.
#' @export

ran2f <- function(x1, x2, V1=NULL, V2=NULL, rrankpct=99, prior=NULL) {
    x1 <- as.factor(x1)
    x2 <- as.factor(x2)
    if(!is.null(V1)) {
        if(is.null(rownames(V1))) {
            stop("There are no rownames on V1")
        }
        EVdecomp <- eigen(V1)
        colnames(EVdecomp$vectors) <- paste("evec",1:ncol(EVdecomp$vectors),sep="")
        rownames(EVdecomp$vectors) <- rownames(V1)
        attr(x1, "evalues") <- EVdecomp$values
        attr(x1, "evectors") <- EVdecomp$vectors
    }
    if(!is.null(V2)) {
        if(is.null(rownames(V2))) {
            stop("There are no rownames on V2")
        }
        EVdecomp <- eigen(V2)
        colnames(EVdecomp$vectors) <- paste("evec",1:ncol(EVdecomp$vectors),sep="")
        rownames(EVdecomp$vectors) <- rownames(V2)
        attr(x2, "evalues") <- EVdecomp$values
        attr(x2, "evectors") <- EVdecomp$vectors
    }
    attr(x1, "factor2") <- x2
    if (!is.null(prior)) attr(x1,"prior") <- prior
    attr(x1,"rrankpct") <- rrankpct
    return(x1)
}
