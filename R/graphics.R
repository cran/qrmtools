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
##' @param ran range over which to plot the colors (can be used to force (-1, 1),
##'        for example)
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
##'       - Check via examples on ?matrix_plot
matrix_plot <- function(x, ran = range(x, na.rm = TRUE), ylim = rev(c(0.5, nrow(x) + 0.5)),
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
    ## Plot
    ## Note: t() does not affect the above as we only determine ranges/max/min there
    levelplot(t(x), ylim = ylim, # correct (x-axis = column; y-axis = row (and decreasing))
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
##' @param FUN hypothesized *quantile* function (vectorized)
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

##' @title Plot of (a) Step Function(s)
##' @param x 2-column matrix (with x and y values) or list of such (in which case)
##'        each element corresponds to one discrete distribution function.
##' @param yleft (vector of) limit(s) to the left of the plotted functions
##' @param do.points logical indicating whether points are to be plotted
##' @param log either "" or "x", indicating whether a logarithmic x-axis is used
##' @param xlim x-axis limits
##' @param ylim y-axis limits
##' @param col color
##' @param main title
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param ... additional arguments passed to the underlying method plot.stepfun()
##' @return see plot.stepfun()
##' @author Marius Hofert
stepfun_plot <- function(x, yleft, do.points = NA, log = "",
                         xlim = NULL, ylim = NULL, col = NULL, main = "",
                         xlab = "x", ylab = "Function value at x", ...)
{
    ## Basics
    if(!is.list(x)) x <- list(x)
    len <- length(x)
    stopifnot(len >= 1)
    if(length(do.points) == 1) do.points <- rep(do.points, len)
    if(length(yleft) == 1) yleft <- rep(yleft, len)
    if(is.null(col)) col <- "black"
    if(length(col) == 1) col <- rep(col, len)
    stopifnot(length(yleft) == len, length(do.points) == len, is.logical(do.points),
              length(col) == len)

    ## x- and y-values
    sf <- vector("list", len)
    for(l in 1:len) {
        x. <- x[[l]]
        if(!is.matrix(x.))
            stop(paste0("Element ",l," of 'x' is not a matrix"))
        if(ncol(x.) != 2)
            stop(paste0("Element ",l," of 'x' is not a 2-column matrix"))
        ## x-values
        ord <- order(as.numeric(x.[,1])) # get order of x
        x.. <- x.[,1][ord] # ordered x (as required by stepfun())
        ## y-values
        y.. <- as.numeric(x.[,2])[ord] # order as x (to not change the function)
        y.. <- c(yleft[l], y..) # y needs to be one element longer than x
        ## Determine if points shall be used
        if(is.na(do.points[l])) do.points[l] <- length(y..) <= 100
        ## Keep result object
        sf[[l]] <- list(FUN = stepfun(x = x.., y = y..), # stepfun() object
                        xlim = range(x.., na.rm = TRUE), # x-axis limits
                        ylim = range(y.., na.rm = TRUE), # y-axis limits
                        do.points = do.points[l], # whether we shall show points for jumps
                        col = col[l]) # color
    }

    ## Axis limits
    if(is.null(xlim)) {
        xlims <- t(sapply(1:len, function(l) sf[[l]]$xlim)) # (len, 2)-matrix
        xlim <- c(min(xlims[,1]), max(xlims[,2]))
    }
    x0 <- xlim[1] - 0.04 * diff(xlim) # x-value at left endpoint for the line extension y = yleft
    if(is.null(ylim)) {
        ylims <- t(sapply(1:len, function(l) sf[[l]]$ylim))
        ylim <- c(min(ylims[,1]), max(ylims[,2]))
    }

    ## Plot
    for(l in 1:len) {
        obj <- sf[[l]]
        if(l == 1) {
            plot(obj$FUN, log = log, xlim = xlim, ylim = ylim,
                 do.points = obj$do.points, col = obj$col,
                 main = main, xlab = xlab, ylab = ylab, ...) # ... can pass 'verticals = FALSE' (exists for plot.stepfun())
        } else {
            plot(obj$FUN, add = TRUE, # cannot have 'log' here then
                 do.points = obj$do.points, col = obj$col, ...)
        }
        ## Extend the plot(s) to the left (line segment with y = yleft).
        ## This is especially required if log == "x" and d >= 1 or log == "" and d > 1.
        segments(x0 = x0, y0 = yleft,
                 x1 = obj$xlim[1], y1 = yleft, col = obj$col) # doesn't pass '...' though (could remove those of plot.stepfun())
        ## This overplots possibly existing line segments for y = yleft but that's okay
        ## (and in fact the easiest solution because if log == "x" and d > 1, some line
        ## segments are already plotted anyways)
    }
}

##' @title Plotting (an) Empirical Distribution Function(s)
##' @param x vector or list of such; if a list, each element corresponds
##'        to one empirical distribution function
##' @param do.points see discrete_df_plot()
##' @param log see discrete_df_plot()
##' @param xlim see discrete_df_plot()
##' @param ylim see discrete_df_plot()
##' @param col see discrete_df_plot()
##' @param main see discrete_df_plot()
##' @param xlab see discrete_df_plot()
##' @param ylab see discrete_df_plot()
##' @param ... see discrete_df_plot()
##' @return see discrete_df_plot()
##' @author Marius Hofert
edf_plot <- function(x, yleft = 0, do.points = NA, log = "", xlim = NULL, ylim = NULL,
                     col = NULL, main = "", xlab = "x", ylab = "Distribution function at x", ...)
{
    if(!is.list(x)) x <- list(x)
    len <- length(x)
    stopifnot(len >= 1)
    ## Now add the y-values required for discrete_df_plot()
    x <- lapply(1:len, function(l) {
        x. <- sort(as.numeric(x[[l]]))
        cbind(x = x., y = ecdf(x.)(x.))
    })
    ## Plot
    stepfun_plot(x, yleft = yleft, do.points = do.points, log = log,
                 xlim = xlim, ylim = ylim, col = col,
                 main = main, xlab = xlab, ylab = ylab, ...)
}
