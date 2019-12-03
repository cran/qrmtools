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
x <- c(5, 2, 4, 2, 3, 2, 2, 2, 1, 2) # example data
edf_plot(x, verticals = FALSE) # the 'mathematical' version
edf_plot(x, do.points = FALSE) # good for many sample points
edf_plot(x) # the default
edf_plot(x, log = "x") # logarithmic
(edf_plot(x, log = "x")) # ... with return value
## Note: flat part below first jump cannot be shown in log-scale
}
\keyword{hplot}