\name{tail_plot}
\alias{tail_plot}
\title{Plot of a Non-Parametric Tail Estimator}
\description{
  Plot a non-parametric tail distribution, possibly overlaid
  with the Smith estimator.
}
\usage{
tail_plot(x, threshold, shape = NULL, scale = NULL, q = NULL,
          length.out = 129, lines.args = list(),
          log = "xy", xlab = "Value", ylab = "Tail probability", ...)
}
\arguments{
  \item{x}{\code{\link{numeric}} vector of data.}
  \item{threshold}{\code{numeric(1)} giving the threshold \eqn{u} above
   which the tail (starts and) is to be plotted.}
  \item{shape}{\code{NULL} or the GPD shape parameter \eqn{\xi}{xi}
    (typically obtained via \code{\link{fit_GPD_MLE}()}).}
  \item{scale}{\code{NULL} or the GPD shape parameter \eqn{\beta}{beta}
    (typically obtained via \code{\link{fit_GPD_MLE}()}).}
  \item{q}{\code{\link{numeric}} vector of evaluationn points of the
    Smith estimator (semi-parametric GPD-based tail estimator in the POT
    method). If \code{is.null(q)}, it is determined internally as an
    equidistant sequence of length \code{length.out} between the
    smallest and largest exceedance (taken equidistant in log-scale if
    \code{log} contains \code{"x"}).}
  \item{length.out}{length of \code{q}.}
  \item{lines.args}{\code{\link{list}} of arguments passed to the
    underlying \code{\link{lines}()}.}
  \item{log}{\code{\link{character}} indicating whether
    logarithmic axes are to be used.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{plot}()}.}
}
\value{
  If both \code{shape} and \code{scale} are provided, \code{tail_plot()}
  overlays the non-parametric tail estimator (evaluated at the exceedances)
  with the corresponding GPD. In this case, \code{tail_plot()} invisibly
  returns a list with two two-column matrices, once containing the x-values
  and y-values of the non-parametric tail estimator and once containing the
  x-values and y-values of the Smith estimator. If \code{shape} or
  \code{scale} are \code{NULL}, \code{tail_plot()} invisibly returns
  a two-column matrix with the x-values and y-values of the non-parametric
  tail estimator.
}
\author{Marius Hofert}
\examples{
data(fire)
u <- 10 # threshold choice
tail_plot(fire, threshold = u, log = "", type = "b") # => need log-scale
tail_plot(fire, threshold = u, type = "s") # as a step function
fit <- fit_GPD_MLE(fire[fire > u] - u) # fit GPD to excesses (POT method)
tail_plot(fire, threshold = u, # without log-scale
          shape = fit$par[["shape"]], scale = fit$par[["scale"]], log = "")
tail_plot(fire, threshold = u, # highlights linearity
          shape = fit$par[["shape"]], scale = fit$par[["scale"]])
}
\keyword{hplot}