#' Predict new response variables from the bayz object
#'
#' Predict new data points from a fitted bayz model, if no new data is given returns the fitted response variables. If new data is given, it gives predicted response variables.
#'
#' @param object  A formula describing the model.
#' @param ...     Additional parameters passed onto the Model function.
#'
#' @return fitted
#' @export
predict.bayz <- function(object, ...){
    parameters <- list(...)

    findDF <- Vectorize(function(x) class(x)=="data.frame")

    data_included <- FALSE

    modeldata <- parameters$data
    if(is.null(modeldata)){
        ## Find any data.frame included in the predict function
        ## Makes it possible to provide a data.frame without using the name data=df
        if(sum(unlist(findDF(parameters)))==1) {
            modeldata <- parameters[[findDF(parameters)]]
            data_included
        }
    } else {
        data_included <- TRUE
    }

    if(data_included){
        ## Data was provided

        ## check for matching variable names between

        ### compute the new response based on the fitted values! ###

        predictions <- 0

        return(predictions) ## Return new predicted values

    } else {
        ## No Data has been provided

        ## Return fitted values

        fitted_values <- 0

        return(fitted_values) ## Return
    }

}
