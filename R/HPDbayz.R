#' HPD intervals accounting for boundaries on parameter space
#'
#' @param x    samples from a (posterior) distribution; if the input is a matrix or dataframe with multiple columns
#'             an HPD interval is computed for every column
#' @param prob the target probability (0,1) to collect in the HPD interval; as the underlying density estimate
#'             is discrete, HPDbayz will obtain the interval that includes just more than the given 'prob'
#' @param bound a single character string or vector of character strings giving the types of bounds on the parameters,
#'              which must be any of "none", "var" (parameter is restricted to (0,INF)), or "prob" (parameter is
#'              restricted to (0,1)). 
#'
#' @return A matrix with 2 columns for the lower and upper values of the HPD interval, and a row for every input column.
#' @export
HPDbayz <- function(x, prob=0.95, bound="none") {
    x = as.matrix(x)
    if( prob < 0 | prob > 1) {
      stop("In HPDBayz the HPD probability is not from 0-1")
    }
    if ( !( (ncol(x) == length(bound)) | (ncol(x)>1 & length(bound)==1) ) ) {
       stop("In HPDBayz the length of 'bound' does not match the number of columns in the input data")
    }
    if ( ncol(x)>1 & length(bound)==1 ) bound = rep(bound,ncol(x))
    result = matrix(0,ncol(x),2)
    for(col in 1:ncol(x)) {
        standarddens = density(x[,col])  # the density without any reflections for bound "none"
        if(bound[col]=="var") {
           xreflected = c(-x[,col],x[,col])
           densestim = density(xreflected,bw=standarddens$bw,from=0)
        }
        else if (bound[col]=="prob") {
           xreflected = c(-x[,col],x[,col],2-x[,col])
           densestim = density(xreflected,bw=standarddens$bw,from=0, to=1)
        }
        else if (bound[col]=="none") {
           densestim = standarddens
        }
        else {
           stop("In HPDbayz bound-type is not any of 'none','var' or 'prob'")
        }
        n = length(densestim$x)
        densestim$y = densestim$y/(sum(densestim$y)-0.5*densestim$y[1]-0.5*densestim$y[n])
        hpdleft = 1
        hpdright = n
        tailprob = 0
        leftbinprob = 0.5*(densestim$y[hpdleft]+densestim$y[hpdleft+1])
        rightbinprob = 0.5*(densestim$y[hpdright]+densestim$y[hpdright-1])
        while ( (tailprob - min(leftbinprob,rightbinprob) < (1-prob) ) ) {
            if (leftbinprob < rightbinprob) {
                tailprob = tailprob + leftbinprob
                hpdleft = hpdleft + 1
                leftbinprob = 0.5*(densestim$y[hpdleft]+densestim$y[hpdleft+1])
            }
            else {
                tailprob = tailprob + rightbinprob
                hpdright = hpdright -1
                rightbinprob = 0.5*(densestim$y[hpdright]+densestim$y[hpdright-1])
            }
        }
        result[col,] = c(densestim$x[hpdleft],densestim$x[hpdright])
   }
   colnames(result) = c("HPDleft","HPDright")
   rownames(result) = colnames(x)
   return(result)
}
