\name{matrix_plot}
\alias{matrix_plot}
\title{Graphical Tool for Visualizing Matrices}
\description{
  Plot of a matrix.
}
\usage{
matrix_plot(x, ylim = rev(c(0.5, nrow(x) + 0.5)),
            xlab = "Column", ylab = "Row",
            scales = list(alternating = c(1,1), tck = c(1,0),
                          x = list(at = pretty(1:ncol(x)), rot = 90),
                          y = list(at = pretty(1:nrow(x)))),
            at = NULL, colorkey = NULL, col = c("royalblue3", "white", "maroon3"),
            col.regions = NULL, ...)
}
\arguments{
  \item{x}{\code{\link{matrix}}-like object.}
  \item{ylim}{y-axis limits in reverse order (for the rows to appear 'top down').}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{scales}{see \code{\link{levelplot}()}; if \code{\link{NULL}},
    labels and ticks are omitted.}
  \item{at}{see \code{\link{levelplot}()}. If \code{\link{NULL}}, a
    useful default is computed based on the given values in \code{x}.}
  \item{colorkey}{see \code{\link{levelplot}()}. If \code{\link{NULL}}, a
    useful default is computed based on \code{at}.}
  \item{col}{\code{\link{vector}} of length two (if all values of
    \code{x} are non-positive or all are non-negative; note that also a
    vector of length three is allowed in this case) or three (if
    \code{x} contains negative and positive values) providing the color
    key's default colors.}
  \item{col.regions}{see \code{\link{levelplot}()}. If \code{\link{NULL}}, a
    useful default is computed based on \code{at}.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{levelplot}()}.}
}
\value{
  The plot, a Trellis object.
}
\details{
  Plot of a matrix.
}
\author{Marius Hofert}
\examples{
## Generate a random correlation matrix
d <- 50
L <- diag(1:d)
set.seed(271)
L[lower.tri(L)] <- runif(choose(d,2))
Sigma <- L %*% t(L) # ... real, symmetric, positive definite (covariance) matrix Sigma
P <- cor(Sigma)

## Default
matrix_plot(P)

## Default if nonnegative
matrix_plot(abs(P))

## Without diagonal
P. <- abs(P)
diag(P.) <- NA
matrix_plot(P.)

## Default if nonpositive
matrix_plot(-abs(P))

## Extending the color key to [-1,1] with darker color for |rho| >> 0
## Note: When specifying 'at', one most likely also wants to provide 'col.regions'
matrix_plot(P, at = seq(-1, 1, length.out = 200),
            col.regions = grey(c(seq(0, 1, length.out = 100), seq(1, 0, length.out = 100))))

## An example with overlaid lines
library(lattice)
my_panel <- function(...) {
    panel.levelplot(...)
    panel.abline(h = c(10, 20), v = c(10, 20), lty = 2)
}
matrix_plot(P, panel = my_panel)
}
\keyword{hplot}