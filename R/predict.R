#' Obtain / retrieve predictictions from a bayz model for NA response values in the original data input 
#'
#' Prediction in bayz currently only supports prediction of the NA responses that were in the
#' original data frame. Predict is simply a wrapper to retrieve them from the output. 
#' To predict data points that were not in the original data frame for analysis, the training model
#' needs to be re-run with these data points added in the training data.
#' 
#' @param object  An output object of a bayz model run of class 'bayz'.
#' @param id      An optional ID that can be attached to the predicted values for reference.
#'                This should be a vector with length and order matching the original input data frame
#'                - typically it is a column of the original data frame such as id=mydata$myID.
#'                Omit for not attaching any reference ID.
#' @param ...     Additional parameters passed onto the Model function.
#'
#' @return fitted
#' @export
predict.bayz <- function(object, id=NULL, ...){
    resid = object$Residuals
    if(!is.null(id)) {
        if(length(id) != nrow(resid)) {
           stop("Error: the length of the supplied id-vector does not match the original input data")
        }
        resid =  cbind(resid,id)
    }
    NAresiduals = is.na(resid[,1])   # residuals are set to NA where response was missing
    if (sum(NAresiduals)==0) {
        cat("There are no predicted values because there seem to have been no NA in the original data\n")
        return(NULL)
    }
    predicted = resid[NAresiduals,-1]
    return(as.data.frame(predicted))
}
