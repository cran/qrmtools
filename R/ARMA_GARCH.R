### ARMA-GARCH #################################################################

##' @title Fitting ARMA-GARCH processes
##' @param x A matrix-like data structure
##' @param ugarchspec.list An object of class uGARCHspec (as returned by
##'        ugarchspec()) or a list of such
##' @param verbose A logical indicating whether verbose output is given
##' @param ... Additional arguments passed to the underlying ugarchfit() or HoltWinters()
##' @return A list of length equal to the number of columns of x containing
##'         the fitted objects
##' @author Marius Hofert
fit_ARMA_GARCH <- function(x, ugarchspec.list = ugarchspec(), verbose = TRUE, ...)
{
    ## Checking and expanding ugarchspec.list to a list of length d
    if(!is.matrix(x)) x <- cbind(x) # is.matrix() is also true for 'xts' objects
    d <- ncol(x)
    isGARCHspec <- function(spec) is(spec, class2 = "uGARCHspec")
    if(isGARCHspec(ugarchspec.list)) ugarchspec.list <- rep(list(ugarchspec.list), d) # repeat
    if(is.list(ugarchspec.list))
        stopifnot(length(ugarchspec.list) == d && all(sapply(ugarchspec.list, isGARCHspec)))
    stopifnot(is.logical(verbose))

    ## Iterate over all time series
    fit  <- vector("list", length = d)
    warn <- vector("list", length = d)
    err  <- vector("list", length = d)
    if(verbose) {
        pb <- txtProgressBar(max = d, style = if(isatty(stdout())) 3 else 1)
        on.exit(close(pb)) # on exit, close progress bar
    }
    for(j in seq_len(d)) {
        res <- catch(ugarchfit(ugarchspec.list[[j]], data = x[,j], ...)) # fitting
        if(!is.null(res$value)) fit[[j]] <- res$value
        if(!is.null(res$warning)) warn[[j]] <- res$warning
        if(!is.null(res$error)) err[[j]]  <- res$error
        if(verbose) setTxtProgressBar(pb, j) # update progress bar
    }

    ## Return
    list(fit = fit, warning = warn, error = err)
}
