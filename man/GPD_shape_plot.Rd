\name{GPD_shape_plot}
\alias{GPD_shape_plot}
\title{Fitted GPD Shape as a Function of the Threshold}
\description{
  Fit GPDs to various thresholds and plot the fitted GPD shape as
  a function of the threshold.
}
\usage{
GPD_shape_plot(x, thresholds = seq(quantile(x, 0.5), quantile(x, 0.99), length.out = 33),
               conf.level = 0.95, lines.args = list(lty = 2),
               xlab2 = "Excesses", xlab = "Threshold",
               ylab = "Estimated GPD shape parameter with confidence intervals", ...)
}
\arguments{
  \item{x}{\code{\link{numeric}} vector of data.}
  \item{thresholds}{\code{\link{numeric}} vector of thresholds for which
    to fit a GPD to the excesses.}
  \item{conf.level}{confidence level of the confidence intervals.}
  \item{lines.args}{\code{\link{list}} of arguments passed to
    the underlying \code{\link{lines}()} for drawing the confidence intervals.}
  \item{xlab2}{label of the secondary x-axis.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{plot}()}.}
}
\value{
  \code{\link{invisible}()}
}
\details{
  Such plots can be used in the peaks-over-threshold method for
  assessing possible threshold choices. Ideally, the plot should be
  (roughly) flat after the chosen threshold (indicating stability of the
  threshold choice).
}
\author{Marius Hofert}
\examples{
data(fire)
GPD_shape_plot(fire)
abline(v = c(10, 20)) # possible threshold choices from mean_excess_plot()
}
\keyword{parametric}
\keyword{hplot}