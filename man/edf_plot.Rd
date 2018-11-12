\name{edf_plot}
\alias{edf_plot}
\title{Plot of an Empirical Distribution Function}
\description{
  Plotting an empirical distribution function.
}
\usage{
edf_plot(x, do.points = length(x) <= 100, log = "",
         xlim = range(x, na.rm = TRUE),
         main = "", xlab = "x", ylab = "Distribution function at x", ...)
}
\arguments{
  \item{x}{\code{\link{numeric}} vector of data (of which the empirical
    distribution function is to be plotted).}
  \item{do.points}{\code{\link{logical}} indicating whether points are
    to be plotted; see \code{?plot.stepfun}.}
  \item{log}{\code{\link{character}} indicating whether
     a logarithmic x-axis is used (available are only \code{""} and \code{"x"}).}
  \item{xlim}{x-axis limits; default \code{range(x)} avoids possible failure if
    \code{log = "x"} and data points are all positive
    (\code{\link{plot.stepfun}()} extends the range, possibly below 0).}
  \item{main}{title.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{plot.stepfun}()}.}
}
\value{
  Returns the return value of the underlying \code{\link{plot.stepfun}()}.
}
\author{Marius Hofert}
\examples{
## Example 1
x <- c(5, 2, 4, 2, 3, 2, 2, 2, 1, 2) # test data
edf_plot(x, verticals = FALSE) # the 'mathematical' version
edf_plot(x, do.points = FALSE) # good for many sample points
edf_plot(x) # the default
edf_plot(x, log = "x") # logarithmic
(edf_plot(x, log = "x")) # ... with return value
## Note: flat part below first jump cannot be shown in log-scale

## Example 2
data(fire)
u <- 10 # threshold
exceed <- fire[fire > u] # exceedances
excess <- exceed - u # excesses
fit <- fit_GPD_MLE(excess) # fit GPD to excesses

## Plot empirical distribution function of excesses
## (partially with fitted GPD overlaid)
edf_plot(excess)
res <- edf_plot(excess, log = "x")
z <- tail(res$t, n = -1)
lines(z, pGPD(z, shape = fit$par[["shape"]], scale = fit$par[["scale"]]))

## Plot empirical distribution function of exceedances
## (partially with fitted GPD overlaid)
edf_plot(exceed)
res <- edf_plot(exceed, log = "x")
z <- tail(res$t, n = -1)
lines(z, pGPD(z-u, shape = fit$par[["shape"]], scale = fit$par[["scale"]]))

## Note: Q-Q plots are typically more meaningful
qf <- function(p) # quantile function of df
    qGPD(p, shape = fit$par[["shape"]], scale = fit$par[["scale"]])
qq_plot(excess, FUN = qf)
qq_plot(exceed, FUN = function(p) u + qf(p))
}
\keyword{hplot}