#' Plotting a bayz object
#'
#' Produced summary statistics of the bayz object and returns a summarybayz object.
#'
#' @param x    bayz object
#' @param ...  additional parameters
#'
#' @return nothing, plots the convergence of the parameters
#' @import graphics
#' @export
plot.bayz <- function(x, ...){
    npar <- ncol(x$Samples)

    ncolnrow <- function(n){
        if(n==1) return(c(1,1))
        grid <- ncolnrow(n-1)
        if(prod(grid)>=n){
            return(grid)
        } else {
            argmin = which(grid==min(grid))
            grid[argmin[1]] <- grid[argmin[1]]+1
            return(grid)
        }
    }

    par(mfrow=ncolnrow(npar),
        mar=c(3.5,2.5,3,2), mgp=c(1.5,0.5,0))
    for(j in 1:npar){
        variable_name <- rev(unlist(strsplit(colnames(x$Samples)[j],
                                             split=".", fixed=TRUE)))[1]
        plot(1:length(x$Samples[,j]), x$Samples[,j], main=variable_name,
             xlab = "Number of cycles", ylab = colnames(x$Samples)[j], pch=16)
    }
}
