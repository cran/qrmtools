### Plot of fitted GEV shape as a function of the block size ###################

##' @title Fitted GEV Shape as a Function of the Block Size
##' @param x numeric vector of data
##' @param blocksize numeric vector of block sizes for which to fit a GEV to
##'        the block maxima
##' @param estimate.cov logical indicating whether confidence intervals are
##'        computed
##' @param conf.level confidence level of the confidence intervals
##' @param lines.args list of arguments passed to underlying lines() for
##'        drawing the confidence intervals
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param xlab2 label of the secondary x-axis
##' @param plot logical indicating whether a plot is done
##' @param ... additional arguments passed to the underlying plot()
##' @return invisibly returns the block sizes, the list of corresponding block
##'         maxima and the list of fitted GEV distributions.
##' @author Marius Hofert
GEV_shape_plot <- function(x, blocksize = tail(pretty(seq_len(length(x)/20), n = 64), -1),
                           estimate.cov = TRUE, conf.level = 0.95,
                           lines.args = list(lty = 2), xlab = "Block size",  ylab = NULL,
                           xlab2 = "Number of blocks", plot = TRUE, ...)
{
    ## Checks
    stopifnot(length(blocksize) >= 2, is.logical(estimate.cov), 0 <= conf.level, conf.level <= 1,
              is.logical(plot))
    ## Compute block maxima for block size n (last one not filled)
    x <- as.numeric(x)
    block_maxima <- function(n) sapply(split(x, ceiling(seq_along(x)/n)), max)
    ## Block maxima data for each considered block size
    blockmax <- lapply(blocksize, function(n) block_maxima(n))
    ## Fit GEV models to the given block sizes
    fits <- lapply(blockmax, fit_GEV_MLE, estimate.cov = estimate.cov) # for all block sizes, fit a GEV to the block maxima
    ## Extract the fitted shape parameters and compute CIs
    xi <- sapply(fits, function(f) f$par[["shape"]])
    if(estimate.cov) {
        xi.SE <- sapply(fits, function(f) f$SE[["shape"]])
        q <- qnorm(1 - (1 - conf.level)/2)
        xi.CI.low <- xi - xi.SE * q
        xi.CI.up  <- xi + xi.SE * q
        ylim <- range(xi, xi.CI.low, xi.CI.up)
    } else ylim <- range(xi)
    ## Plot
    if(plot) {
        if(is.null(ylab)) ylab <- paste0("Estimated GEV shape parameter",
                                         if(estimate.cov) paste0(" with ", 100*conf.level,
                                                                 "% confidence intervals") else "")
        plot(blocksize, xi, type = "l", ylim = ylim,
             xlab = xlab, ylab = ylab, ...)
        if(estimate.cov) {
            do.call(lines, args = c(list(x = blocksize, y = xi.CI.low), lines.args))
            do.call(lines, args = c(list(x = blocksize, y = xi.CI.up),  lines.args))
        }
        pb <- pretty(blocksize) # where actual x labels are (even if those block sizes are not considered)
        axis(3, at = pb, labels = sapply(pb, function(n) ceiling(length(x)/n)))
        mtext(xlab2, side = 3, line = 3)
    }
    ## Return
    invisible(list(block.sizes = blocksize, block.maxima = blockmax, GEV.fits = fits))
}
