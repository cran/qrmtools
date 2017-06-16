### ARMA-GARCH #################################################################

##' @title Fitting ARMA-GARCH processes
##' @param x matrix-like data structure
##' @param ugarchspec.list object of class uGARCHspec (as returned by
##'        ugarchspec()) or a list of such
##' @param solver string indicating the solver used; see ?ugarchfit
##' @param verbose logical indicating whether verbose output is given
##' @param ... additional arguments passed to the underlying ugarchfit() or HoltWinters()
##' @return list of length equal to the number of columns of x containing
##'         the fitted objects
##' @author Marius Hofert
##' @note The default ugarchspec.list fits an ARMA(1,1)-GARCH(1,1) with N(0,1)
##'       standardized residuals
fit_ARMA_GARCH <- function(x, ugarchspec.list = ugarchspec(), solver = "hybrid",
                           verbose = TRUE, ...)
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
        res <- catch(ugarchfit(ugarchspec.list[[j]], data = x[,j], solver = solver,
                               ...)) # fitting
        if(!is.null(res$value)) fit[[j]] <- res$value
        if(!is.null(res$warning)) warn[[j]] <- res$warning
        if(!is.null(res$error)) err[[j]]  <- res$error
        if(verbose) setTxtProgressBar(pb, j) # update progress bar
    }

    ## Return
    list(fit = fit, warning = warn, error = err)
}
