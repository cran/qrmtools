### General graphical tools ####################################################

##' @title Image Indicating NAs in a Data Set
##' @param x matrix (ideally an xts object)
##' @param col colors for NA and non-NA, respectively
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param text see mtext()
##' @param side see mtext()
##' @param line see mtext()
##' @param adj see mtext()
##' @param ... additional arguments passed to image()
##' @return invisible()
##' @author Marius Hofert
NA_plot <- function(x, col = c("black", "white"), xlab = "Time", ylab = "Component",
                    text = "Black: NA; White: Available data", side = 4, line = 1, adj = 0,
                    ...)
{
    stopifnot(is.matrix(x))
    x. <- if(inherits(x, "xts")) {
        index(x) # use the time points
    } else {
        rn <- rownames(x)
        if(is.null(rn)) seq_len(nrow(x)) else rn # if available, use row names, otherwise numbers
    }
    image(x = x., y = seq_len(ncol(x)), z = is.na(x),
          col = rev(col), xlab = xlab, ylab = ylab, ...)
    if(inherits(text, "call") || nchar(text) > 0)
        mtext(text, side = side, line = line, adj = adj)
    invisible()
}

##' @title Plot of a Matrix
##' @param x matrix
##' @param ylim y-axis limits in reverse order (for the rows to appear 'top down')
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param scales see levelplot(); if NULL, labels and ticks are omitted
##' @param at see levelplot()
##' @param colorkey see levelplot()
##' @param col default colors for the color key
##' @param col.regions see levelplot()
##' @param ... additional arguments passed to levelplot()
##' @return the plot, a Trellis object
##' @author Marius Hofert
##' @note - Another option would be:
##'         corrplot::corrplot(err, method = "color", col = grey(seq(0.4, 1, length.out = 200)),
##'                            tl.col = "black", is.corr = FALSE)
##'       - Check via example on ?matrix_plot
matrix_plot <- function(x, ylim = rev(c(0.5, nrow(x) + 0.5)),
                        xlab = "Column", ylab = "Row",
                        scales = list(alternating = c(1,1), tck = c(1,0),
                                      x = list(at = pretty(1:ncol(x)), rot = 90),
                                      y = list(at = pretty(1:nrow(x)))),
                        at = NULL, colorkey = NULL, col = c("royalblue3", "white", "maroon3"),
                        col.regions = NULL, ...)
{
    stopifnot(is.matrix(x), (nr <- nrow(x)) >= 1,
              is.null(scales) || is.list(scales), is.numeric(at) || is.null(at),
              is.list(colorkey) || is.null(colorkey))
    ## Determine colors for the color key
    ran <- range(x, na.rm = TRUE)
    if(all(ran >= 0)) { # all non-NA entries >= 0
        if(is.null(at)) at <- seq(0, ran[2], length.out = 200)
        if(is.null(col.regions))
           col.regions <- colorRampPalette(c(col[2], col[3]), space = "Lab")(200)
    } else if(all(ran <= 0)) { # all non-NA entries <= 0
        lcol <- length(col)
        stopifnot(2 <= lcol, lcol <= 3)
        if(is.null(at)) at <- seq(ran[1], 0, length.out = 200)
        if(is.null(col.regions))
           col.regions <- colorRampPalette(c(col[1], col[2]), space = "Lab")(200)
    } else { # entries < 0 && entries > 0
        lcol <- length(col)
        stopifnot(2 <= lcol, lcol <= 3)
        if(lcol == 2) col <- c(0, col)
        ## Scale so that 0 gets the 'middle' color col[2]
        frac1 <- -ran[1]/diff(ran)
        frac2 <- ran[2]/diff(ran) # => frac1 + frac2 = 1 (and both in [0,1])
        if(is.null(at)) at <- seq(ran[1], ran[2], length.out = 200)
        if(is.null(col.regions))
            col.regions <- c(colorRampPalette(col[1:2], space = "Lab")(floor(200*frac1)),
                             colorRampPalette(col[2:3], space = "Lab")(ceiling(200*frac2)))
    }
    if(min(x, na.rm = TRUE) < at[1] || max(x, na.rm = TRUE) > at[length(at)])
        stop("'x' values outside the range spanned by 'at'. Choose 'at' appropriately.")
    levelplot(x, ylim = ylim, # correct (x-axis = column; y-axis = row (and decreasing))
              xlab = xlab, ylab = ylab, col.regions = col.regions,
              scales = if(is.null(scales)) list(alternating = c(0,0), tck = c(0,0)) else scales,
              at = at, colorkey = if(is.null(colorkey)) list(at = at) else colorkey, ...)
}

##' @title Density Plot of the Values from a Lower Triangular Matrix
##' @param x matrix
##' @param xlab x-axis label; see plot()
##' @param main title; see plot()
##' @param text see mtext()
##' @param side see mtext()
##' @param line see mtext()
##' @param adj see mtext()
##' @param ... additional arguments passed to plot()
##' @return invisible()
##' @author Marius Hofert
matrix_density_plot <- function(x, xlab = "Entries in the lower triangular matrix",
                                main = "", text = NULL, side = 4, line = 1, adj = 0, ...)
{
    if(!is.matrix(x)) x <- rbind(x, deparse.level = 0L)
    d <- ncol(x)
    x.vec <- x[lower.tri(x)] # grab out values from the lower triangular matrix
    dens.x <- density(x.vec) # density
    plot(dens.x, main = main, xlab = xlab, ...) # plot
    rug(x.vec, col = "black") # rugs
    if(is.null(text))
        text <- substitute("Dimension: "*d.*"; Sample size: "*n.*"; Bandwidth: "*b.,
                           list(d. = d, n. = dens.x$n, b. = format(dens.x$bw, digits = 3, scientific = TRUE)))
    if(inherits(text, "call") || nchar(text) > 0)
        mtext(text, side = side, line = line, adj = adj)
    invisible()
}

##' @title P-P Plot
##' @param x data (a vector of convertible to such)
##' @param FUN hypothesized *distribution* function
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param ... additional arguments passed to the underlying plot()
##' @return invisible()
##' @author Marius Hofert
##' @note as.vector() required since sort(x) == x for an 'xts' object!
pp_plot <-  function(x, FUN = pnorm,
                     xlab = "Theoretical probabilities", ylab = "Sample probabilities",
                     ...)
{
    p <- ppoints(length(x)) # theoretical probabilities of sorted data
    y <- FUN(sort(as.vector(x))) # hypothesized quantiles of sorted data
    plot(p, y, xlab = xlab, ylab = ylab, ...)
    abline(a = 0, b = 1)
    invisible()
}

##' @title Q-Q Plot
##' @param x data (a vector or convertible to such)
##' @param FUN hypothesized *quantile* function
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param do.qqline logical indicating whether a Q-Q line is plotted
##' @param method method used to construct the Q-Q line:
##'        "theoretical": theoretically true line which helps
##'                       deciding whether x comes from exactly
##'                       the distribution specified by FUN
##'        "empirical": qqline() which interpolates
##'                     (F^-(0.25), F_n^-(0.25)) and (F^-(0.75)), F_n^-(0.75))
##'                     and which helps deciding whether x comes from
##'                     a location-scale transformed distribution
##'                     specified by FUN.
##'        Example:
##'        z <- rnorm(1000)
##'        qq_plot(z, qnorm) # fine
##'        mu <- 3
##'        sig <- 2
##'        z. <- mu+sig*z
##'        qq_plot(z., qnorm) # not fine
##'        qq_plot((z.-mean(z.))/sd(z.), qnorm) # fine again
##'        qq_plot(z., qnorm, method = "empirical") # fine again
##' @param qqline.args list containing additional arguments passed to the
##'        underlying abline() functions
##' @param ... additional arguments passed to the underlying plot()
##' @return invisible()
##' @author Marius Hofert
##' @note - as.vector() required since sort(x) == x for an 'xts' object!
##'       - This is a convenience function which does the same as
##'         qqplot(FUN(ppoints(length(x))), x)
##'         qqline(x, distribution = function(p) FUN(p))
##'         ... but has nicer default labels, does the Q-Q line by default
##'         and uses the *theoretical* Q-Q line by default, not an
##'         empirically estimated one via
##'         x <- FUN(probs)
##'         y <- quantile(x, probs = c(0.25, 0.75))
##'         slope <- diff(y) / diff(x)
##'         int <- y[1] - slope * x[1]
##'         abline(a = int, b = slope)
qq_plot <-  function(x, FUN = qnorm, xlab = "Theoretical quantiles", ylab = "Sample quantiles",
                     do.qqline = TRUE, method = c("theoretical", "empirical"),
                     qqline.args = NULL, ...)
{
    x. <- x[!is.na(x)] # grab out non-NA data
    q <- FUN(ppoints(length(x.))) # theoretical quantiles of sorted data
    y <- sort(as.vector(x.)) # compute order statistics (sample quantiles)
    plot(q, y, xlab = xlab, ylab = ylab, ...)
    if(do.qqline) {
        method <- match.arg(method)
        switch(method,
        "theoretical" = { # intercept 0, slope 1
            if(is.null(qqline.args)) qqline.args <- list(a = 0, b = 1) # defaults
            do.call(abline, qqline.args)
        },
        "empirical" = { # estimates intercept and slope
            if(is.null(qqline.args)) qqline.args <- list()
            arg.lst <- c(list(x., distribution = FUN), qqline.args)
            do.call(qqline, arg.lst)
        },
        stop("Wrong 'method'"))
    }
    invisible()
}

##' @title Plot of an Empirical Distribution Function
##' @param x data of which the empirical distribution function is to be plotted
##' @param do.points logical indicating whether points are to be plotted
##' @param log either "" or "x", indicating whether a logarithmic x-axis is used
##' @param xlim x-axis limits; default avoids possible failure if log = "x"
##'        and data points are all positive (plot.stepfun() extends the range,
##'        possibly below 0)
##' @param main title
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param ... additional arguments passed to the underlying plot method of stepfun();
##'        see plot.stepfun()
##' @return see plot.stepfun()
##' @author Marius Hofert
##' @note MWE:
##'       x <- 1:100/100
##'       plot(stepfun(x = x, y = c(0, ecdf(x)(x))), log = "x")
##'       plot(stepfun(x = x, y = c(0, ecdf(x)(x))), log = "x", xlim = range(x))
##'       plot.stepfun # => extends the range (below 0)
##'       Note: Manually extending the range a little bit to the left
##'             does not make sense (because of log-scale
##'             => artificial extension to the left). Best to leave it like this.
edf_plot <- function(x, do.points = length(x) <= 100, log = "",
                     xlim = range(x, na.rm = TRUE),
                     main = "", xlab = "Value", ylab = "Probability", ...)
{
    x <- sort(as.numeric(x)) # required by ecdf()
    y <- c(0, ecdf(x)(x)) # plot.stepfun() requires 'y' to be one longer than 'x' (y = values *between* x's)
    ## => log = "y" does not make sense anymore at this point
    if(grepl("y", log)) stop('log = "y" not available.')
    sf <- stepfun(x = x, y = y) # ok, does not extend x-range
    plot(sf, do.points = do.points, log = log, xlim = xlim,
         main = main, xlab = xlab, ylab = ylab, ...) # see plot.stepfun()
}
