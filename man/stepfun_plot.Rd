\name{stepfun_plot}
\alias{stepfun_plot}
\alias{edf_plot}
\title{Plot of Step Functions and Empirical Distribution Functions}
\description{
  Plotting step functions and empirical distribution functions.
}
\usage{
stepfun_plot(x, yleft, do.points = NA, log = "",
             xlim = NULL, ylim = NULL, col = NULL, main = "",
             xlab = "x", ylab = "Function value at x", ...)
edf_plot(x, yleft = 0, do.points = NA, log = "", xlim = NULL, ylim = NULL,
         col = NULL, main = "", xlab = "x", ylab = "Distribution function at x", ...)
}
\arguments{
  \item{x}{
    \describe{
      \item{\code{stepfun_plot()}}{2-column \code{\link{matrix}} (with
	x- and y-values) or \code{\link{list}} of such; if a list, each element
        corresponds to one function to plot.}
      \item{\code{edf_plot()}}{\code{\link{numeric}} vector or a
	\code{\link{list}} of such; if a list, each
	element corresponds to one empirical distribution function.}
    }
  }
  \item{yleft}{y-value(s) of the graph(s) extending to the left of the
    first x-value.}
  \item{do.points}{\code{\link{logical}} (vector) indicating whether points are
    to be plotted (defaults to \code{\link{TRUE}} if and only if there
    are less than or equal to 100 data points); see \code{?plot.stepfun}.}
  \item{log}{\code{\link{character}} indicating whether
     a logarithmic x-axis is used.}
  \item{xlim}{x-axis limits.}
  \item{ylim}{y-axis limits.}
  \item{col}{(vector of) color(s).}
  \item{main}{title.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{plot.stepfun}()}.}
}
\value{
  Nothing (plot by side-effect).
}
\author{Marius Hofert}
\examples{
x <- c(5, 2, 4, 2, 3, 2, 2, 2, 1, 2) # example data
edf_plot(x) # the default
edf_plot(x, verticals = FALSE) # the 'mathematical' version
edf_plot(x, do.points = FALSE) # good for many sample points
edf_plot(x, log = "x") # logarithmic; cannot show flat part before first jump
edf_plot(list(x, x+2), col = c("black", "royalblue3")) # plots both empirical distributions
edf_plot(list(x, x[1:4]+2), log = "x", col = c("black", "royalblue3"))
}
\keyword{hplot}