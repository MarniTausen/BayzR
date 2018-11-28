#' prior function
#'
#' Sets a variable to a fixed variable in a linear mixed model.
#'
#' @param scale  Scale
#' @param df     Degrees of freedom
#'
#' @return Factorized version of the vector
#' @export
#'
ichi <- function(scale=0, df=-2){
    prior <- list("scale"=scale, "df"=df, update=FALSE)
    class(prior) <- c("bayzPrior", "list")
    attr(prior, "dist") <- "ichi"
    return(prior)
}
