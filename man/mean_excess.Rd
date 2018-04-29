\name{mean_excess}
\alias{mean_excess_np}
\alias{mean_excess_plot}
\alias{mean_excess_GPD}
\title{Mean Excess}
\description{
  Sample mean excess function, mean excess function of a GPD
  and sample mean excess plot.
}
\usage{
mean_excess_np(x, omit = 3)
mean_excess_plot(x, omit = 3,
                 xlab = "Threshold", ylab = "Mean excess over threshold", ...)
mean_excess_GPD(x, shape, scale)
}
\arguments{
  \item{x}{
    \describe{
      \item{\code{mean_excess_GPD()}}{\code{\link{numeric}} vector of
	evaluation points of the mean excess function of the GPD.}
      \item{otherwise}{\code{\link{numeric}} vector of data.}
    }
  }
  \item{omit}{number \eqn{\ge 1}{>= 1} of unique % differs from QRM::MEplot()!
    last observations to be omitted from the sorted data (as mean excess plot
    becomes unreliable for these observations as thresholds).}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{plot}()}.}
  \item{shape}{GPD shape parameter \eqn{\xi}{xi}.}
  \item{scale}{GPD scale parameter \eqn{\beta}{beta}.}
}
\value{
  \code{mean_excess_np()} returns a two-column matrix giving
  the sorted data without the \code{omit}-largest unique values
  (first column) and the corresponding values of the sample mean excess
  function (second column). It is mainly used in \code{mean_excess_plot()}.

  \code{mean_excess_plot()} returns \code{invisible()}.

  \code{mean_excess_GPD()} returns the mean excess function of a
  generalized Pareto distribution evaluated at \code{x}.
}
\details{
  Mean excess plots can be used in the peaks-over-threshold method for
  choosing a threshold. To this end, one chooses the smallest threshold
  above which the mean excess plot is roughly linear.
}
\author{Marius Hofert}
\examples{
## (Sample) mean excess function
data(fire)
ME <- mean_excess_np(fire)
stopifnot(dim(ME) == c(2164, 2),
          all.equal(ME[nrow(ME),], c(65.707491, 121.066231),
                    check.attributes = FALSE))

## A 'manual' (sample) mean excess plot
plot(ME, xlab = "Threshold", ylab = "Mean excess over threshold")

## (Sample) mean excess plot
mean_excess_plot(fire)
## => Any value in [10, 20] seems reasonable here as threshold choice
##    (one would probably go with 10 to benefit from a larger sample size).

## With mean excess functions of two fitted GPDs overlaid
u <- c(10, 20) # thresholds
fit <- lapply(u, function(u.) fit_GPD_MLE(fire[fire > u.] - u.))
q <- lapply(u, function(u.) seq(u., ME[nrow(ME),"x"], length.out = 129))
MEF.GPD <- lapply(1:2, function(k)
    mean_excess_GPD(q[[k]]-u[k], shape = fit[[k]]$par[["shape"]],
                    scale = fit[[k]]$par[["scale"]]))
mean_excess_plot(fire, ylim = range(ME, unlist(MEF.GPD)))
col <- c("royalblue3", "maroon3")
for(k in 1:2) lines(q[[k]], MEF.GPD[[k]], col = col[k])
legend("bottomright", col = rev(col), lty = rep(1, length(u)), bty = "n",
       legend = as.expression(sapply(rev(seq_along(u)),
       function(k) substitute("Threshold choice"~~u==u., list(u. = u[k])))))
}
\keyword{nonparametric}
\keyword{hplot}