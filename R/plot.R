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
    cyclenr = as.integer(rownames(x$Samples))
    for(j in 1:npar){
        coldata = x$Samples[,j];
        plot(cyclenr, coldata, main=colnames(x$Samples)[j],
             xlab = "Cycle Number", ylab = colnames(x$Samples)[j], pch=16, cex=0.7)
        lines(cyclenr, cumsum(coldata)/seq(1,length(coldata)), col=2, lwd=2)
    }
}
