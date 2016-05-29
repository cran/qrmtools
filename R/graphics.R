### Graphical tools ############################################################

##' @title Image Indicating NAs in a Data Set
##' @param x A matrix (ideally an xts object)
##' @param col The colors for NA and non-NA, respectively
##' @param xlab The x-axis label
##' @param ylab The y-axis label
##' @param text See mtext()
##' @param side See mtext()
##' @param line See mtext()
##' @param adj See mtext()
##' @param ... Additional arguments passed to image()
##' @return invisible()
##' @author Marius Hofert
plot_NA <- function(x, col = c("black", "white"), xlab = "Time", ylab = "Component",
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
##' @param x A matrix
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param scales See levelplot(); if NULL, labels and ticks are omitted
##' @param at See levelplot()
##' @param colorkey See levelplot()
##' @param col The default colors for the color key
##' @param col.regions See levelplot()
##' @param ... Additional arguments passed to levelplot()
##' @return The level plot
##' @author Marius Hofert
##' @note Another option would be:
##'       corrplot::corrplot(err, method="color", col=grey(seq(0.4, 1, length.out=200)),
##'                          tl.col="black", is.corr=FALSE)
plot_matrix <- function(x, xlab = "Column", ylab = "Row",
                        scales = list(alternating = c(1,1), tck = c(1,0), x = list(rot = 90)),
                        at = NULL, colorkey = NULL, col = c("royalblue3", "white", "maroon3"),
                        col.regions = NULL, ...)
{
    stopifnot(is.matrix(x), (nr <- nrow(x)) >= 1,
              is.null(scales) || is.list(scales), is.numeric(at) || is.null(at),
              is.list(colorkey) || is.null(colorkey), length(col) == 3)
    ran <- range(x, na.rm=TRUE)
    if(all(ran >= 0)) {
        if(is.null(at)) at <- seq(0, ran[2], length.out = 200)
        if(is.null(col.regions))
           col.regions <- colorRampPalette(c(col[2], col[3]), space = "Lab")(200)
    } else if(all(ran <= 0)) {
        if(is.null(at)) at <- seq(ran[1], 0, length.out = 200)
        if(is.null(col.regions))
           col.regions <- colorRampPalette(c(col[1], col[2]), space = "Lab")(200)
    } else { # ran[1] < 0 && ran[2] > 0
        ## Scale so that 0 gets the 'middle' color col[2]
        frac1 <- -ran[1]/diff(ran)
        frac2 <- ran[2]/diff(ran) # => frac1 + frac2 = 1 (and both in [0,1])
        if(is.null(at)) at <- seq(ran[1], ran[2], length.out = 200)
        if(is.null(col.regions))
            col.regions <- c(colorRampPalette(c(col[1], col[2]), space = "Lab")(floor(200*frac1)),
                             colorRampPalette(c(col[2], col[3]), space = "Lab")(ceiling(200*frac2)))
    }
    if(min(x, na.rm=TRUE) < at[1] || max(x, na.rm=TRUE) > at[length(at)])
        stop("'x' values outside the range spanned by 'at'. Choose 'at' appropriately.")
    levelplot(t(x)[,nr:1], xlab=xlab, ylab=ylab,
              col.regions=col.regions,
              scales=if(is.null(scales)) list(alternating=c(0,0), tck=c(0,0)) else scales,
              at=at, colorkey=if(is.null(colorkey)) list(at=at) else colorkey, ...)
}

##' @title Density Plot of the Values from a Lower Triangular Matrix
##' @param x A matrix
##' @param xlab The x-axis label
##' @param main The title
##' @param text See mtext()
##' @param side See mtext()
##' @param line See mtext()
##' @param adj See mtext()
##' @param ... Additional arguments passed to plot()
##' @return invisible()
##' @author Marius Hofert
density_plot_matrix <- function(x, xlab = "Entries in the lower triangular matrix",
                                main = "", text = NULL, side = 4, line = 1, adj = 0, ...)
{
    if(!is.matrix(x)) x <- rbind(x, deparse.level=0L)
    d <- ncol(x)
    x.vec <- x[lower.tri(x)] # grab out values from the lower triangular matrix
    dens.x <- density(x.vec) # density
    plot(dens.x, main = main, xlab = xlab, ...) # plot
    rug(x.vec, col = "black") # rugs
    if(is.null(text))
        text <- substitute("Dimension: "*d.*"; Sample size: "*n.*"; Bandwidth: "*b.,
                           list(d. = d, n. = dens.x$n, b. = dens.x$bw))
    if(inherits(text, "call") || nchar(text) > 0)
        mtext(text, side = side, line = line, adj = adj)
    invisible()
}
