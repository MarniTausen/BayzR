#' Print function for bayz class
#'
#' Prints the bayz class
#'
#' @param x    bayz object
#' @param ...  additional parameters
#'
#' @return Nothing, only prints
#' @import stats
#' @export
print.summarybayz <- function(x, ...){
    object <- x
    cat("Summary statistics of Bayz object\n\n")
    cat("   model formula:", deparse(object$terms), "\n\n")
    if(object$nError>0){
        cat("Bayz encountered errors while running:\n")
        for(errormsg in object$Errors){
            cat("  ",errormsg,"\n")
        }
    } else {
        int_vars <- object[["Random_variables"]]

        cat("Random variable estimates:\n")
        print(object$Estimates[int_vars,])

        cat("\n")

        cat("Variance explained (Heritability):\n")
        Heritability_table <- object$Heritability
        percent <- function(x, digits = 2, format = "f", ...) {
            if(is.na(x)) return(x)
            paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
        }
        percent <- Vectorize(percent)
        sum.ignore.na <- function(x){
            total <- 0
            for(i in x){
                if(is.na(i)) next
                total <- total + i
            }
            total
        }
        if(is.null(Heritability_table$Heritability)){
            last_row <- data.frame(Variable="TOTAL",
                        Broad.sense.Heritability=sum(Heritability_table$Broad.sense.Heritability),
                        Narrow.sense.Heritability=sum.ignore.na(Heritability_table$Narrow.sense.Heritability),
                        stringAsFactors=FALSE)
            last_row$stringAsFactors <- NULL
            Heritability_table <- rbind(Heritability_table, last_row)
            Heritability_table$Broad.sense.Heritability <- percent(Heritability_table$Broad.sense.Heritability)
            Heritability_table$Narrow.sense.Heritability <- percent(Heritability_table$Narrow.sense.Heritability)
        } else {
            last_row <- data.frame(Variable="TOTAL",
                                   Heritability=sum(Heritability_table$Heritability),
                                   stringAsFactors=FALSE)
            last_row$stringAsFactors <- NULL
            Heritability_table <- rbind(Heritability_table, last_row)
            Heritability_table$Heritability <- percent(Heritability_table$Heritability)
        }
        print(Heritability_table)
    }
}
