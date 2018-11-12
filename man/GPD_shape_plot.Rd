\name{GPD_shape_plot}
\alias{GPD_shape_plot}
\title{Fitted GPD Shape as a Function of the Threshold}
\description{
  Fit GPDs to various thresholds and plot the fitted GPD shape as
  a function of the threshold.
}
\usage{
GPD_shape_plot(x, thresholds = seq(quantile(x, 0.5), quantile(x, 0.99),
                                   length.out = 65),
               estimate.cov = TRUE, conf.level = 0.95,
               lines.args = list(lty = 2), xlab = "Threshold", ylab = NULL,
               xlab2 = "Excesses", plot = TRUE, ...)
}
\arguments{
  \item{x}{\code{\link{numeric}} vector of data.}
  \item{thresholds}{\code{\link{numeric}} vector of thresholds for which
    to fit a GPD to the excesses.}
  \item{estimate.cov}{\code{\link{logical}} indicating whether
    confidence intervals are to be computed.}
  \item{conf.level}{confidence level of the confidence intervals if
    \code{estimate.cov}.}
  \item{lines.args}{\code{\link{list}} of arguments passed to
    the underlying \code{\link{lines}()} for drawing the confidence intervals.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label (if \code{NULL}, a default is used).}
  \item{xlab2}{label of the secondary x-axis.}
  \item{plot}{\code{\link{logical}} indicating whether a plot is produced.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{plot}()}.}
}
\value{
  Invisibly returns a \code{\link{list}} containing the thresholds
  considered, the corresponding excesses and the fitted GPD
  objects as returned by the underlying \code{\link{fit_GPD_MLE}()}.
}
\details{
  Such plots can be used in the peaks-over-threshold method for
  determining the optimal threshold (as the smallest after which the plot
  is (roughly) stable).
}
\author{Marius Hofert}
\examples{
data(fire)
GPD_shape_plot(fire)
abline(v = c(10, 20), lty = 3) # possible threshold choices from mean_excess_plot()
}
\keyword{parametric}
\keyword{hplot}