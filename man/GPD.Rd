\name{GPD}
\alias{dGPD}
\alias{pGPD}
\alias{qGPD}
\alias{rGPD}
\alias{dPar}
\alias{pPar}
\alias{qPar}
\alias{rPar}
\title{(Generalized) Pareto Distribution}
\description{
  Density, distribution function, quantile function and random variate
  generation for the (generalized) Pareto distribution (GPD).
}
\usage{
dGPD(x, xi, beta, log=FALSE)
pGPD(q, xi, beta, lower.tail=TRUE, log.p=FALSE)
qGPD(p, xi, beta, lower.tail=TRUE, log.p=FALSE)
rGPD(n, xi, beta)

dPar(x, theta, log=FALSE)
pPar(q, theta, lower.tail=TRUE, log.p=FALSE)
qPar(p, theta, lower.tail=TRUE, log.p=FALSE)
rPar(n, theta)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{xi}{GPD shape parameter, a real number.}
  \item{beta}{GPD scale parameter, a positive number.}
  \item{theta}{Pareto parameter, a positive number.}
  \item{lower.tail}{\code{\link{logical}}; if TRUE (default)
    probabilities are \eqn{P(X \le x)}{P(X <= x)} otherwise, \eqn{P(X > x)}.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
}
\value{
  \code{dGPD()} computes the density, \code{pGPD()} the distribution
  function, \code{qGPD()} the quantile function and \code{rGPD()} random
  variates of the generalized Pareto distribution.

  Similary for \code{dPar()}, \code{pPar()}, \code{qPar()} and
  \code{rPar()} for the (standard) Pareto distribution.
}
\details{
  The distribution function of the generalized Pareto distribution is given by
  \deqn{F(x)=\cases{
    1-(1+\xi x/\beta)^{-1/\xi},&if $\xi\neq 0$,\cr
    1-\exp(-x/\beta),&if $\xi=0$,\cr}}{%
    F(x)=1-(1+xi x/beta)^{-1/xi} if xi !=0 and
    1-exp(-x/beta) if xi=0,}
  where \eqn{\beta>0}{beta>0} and \eqn{x\ge0}{x>=0} if \eqn{\xi\ge 0}{xi>=0}
  and \eqn{x\in[0,-\beta/\xi]}{x in [0,-beta/xi]} if \eqn{\xi<0}{xi<0}.

  The distribution function of the (standard) Pareto distribution is given by
  \deqn{F(x)=1-(1+x)^{-\theta},\ x\ge 0,}{F(x)=1-(1+x)^{-theta}, x>=0,}
  where \eqn{\theta>0}{theta>0}.

  In contrast to \code{dGPD()}, \code{pGPD()}, \code{qGPD()} and
  \code{rGPD()}, the functions \code{dPar()}, \code{pPar()},
  \code{qPar()} and \code{rPar()} are vectorized in both their main
  argument and \eqn{\theta}{theta}.
}
\author{Marius Hofert}
\references{
   McNeil, A. J., Frey, R., and Embrechts, P. (2015).
   \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
   Princeton University Press.
}
\examples{
## Basic sanity checks
plot(pGPD(rGPD(1000, xi=0.5, beta=3), xi=0.5, beta=3)) # should be U[0,1]
curve(dGPD(x, xi=0.5, beta=3), from=-1, to=5)
}
\keyword{distribution}