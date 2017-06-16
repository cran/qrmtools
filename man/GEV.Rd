\name{GEV}
\alias{dGEV}
\alias{pGEV}
\alias{qGEV}
\alias{rGEV}
\title{Generalized Extreme Value Distribution}
\description{
  Density, distribution function, quantile function and random variate
  generation for the generalized extreme value distribution (GEV).
}
\usage{
dGEV(x, xi, mu = 0, sigma = 1, log = FALSE)
pGEV(q, xi, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE)
qGEV(p, xi, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE)
rGEV(n, xi, mu = 0, sigma = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{xi}{GEV shape parameter, a real number.}
  \item{mu}{GEV location parameter, a real number.}
  \item{sigma}{GEV scale parameter, a positive number.}
  \item{lower.tail}{\code{\link{logical}}; if TRUE (default)
    probabilities are \eqn{P(X \le x)}{P(X <= x)} otherwise, \eqn{P(X > x)}.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
}
\value{
  \code{dGEV()} computes the density, \code{pGEV()} the distribution
  function, \code{qGEV()} the quantile function and \code{rGEV()} random
  variates of the generalized extreme value distribution.
}
\details{
  The distribution function of the generalized extreme value
  distribution is given by
  \deqn{F(x) = \cases{
    \exp(-(1-\xi(x-\mu)/\sigma)^{-1/\xi}),&if $\xi\neq 0,\ 1+\xi(x-\mu)/\sigma>0$,\cr
    \exp(-e^{-(x-\mu)/\sigma}),&if $\xi = 0$,\cr}}{%
    F(x) = exp(-(1-xi(x-mu)/sigma)^{-1/xi}) if xi != 0, 1+xi*(x-mu)/sigma>0 and
    exp(-e^{-(x-mu)/sigma}) if xi = 0,}
  where \eqn{\sigma>0}{sigma>0}.
}
\author{Marius Hofert}
\references{
   McNeil, A. J., Frey, R., and Embrechts, P. (2015).
   \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
   Princeton University Press.
}
\examples{
## Basic sanity checks
plot(pGEV(rGEV(1000, xi = 0.5), xi = 0.5)) # should be U[0,1]
curve(dGEV(x, xi = 0.5), from = -3, to = 5)
}
\keyword{distribution}