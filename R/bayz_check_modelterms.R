#' Bayz check model terms
#'
#' Checks if a model is prepared to run through bayz.
#' Detects if the support functions, fixf(), ranf() ... etc is being used.
#' Lists all of the
#'
#' @param model   A formula describing the model.
#'
#' @return A fitted bayz model
#' @export
bayz_check_modelterms <- function(model){
    model_terms <- terms(model)
    accepted_functions_terms <- c("fixf(", "ranf(", "freg(", "ran2f(")
    function_names <- "fixf ranf freg ran2f"
    for(label in attr(model_terms, "term.labels")){
        passed <- FALSE
        for(afunc in accepted_functions_terms){
            passed <- grepl(afunc, label, fixed=TRUE)
            if(passed) break
        }
        if(!passed){
            stop(paste0("Term type is not defined: ", label,
            "\n  Please use one of these model function definitions: ",
            function_names))
        }
    }
}
