\name{risk_measures}
\alias{VaR_t}
\alias{ES_t}
\alias{VaR_Par}
\alias{ES_Par}
\title{Risk Measures}
\description{
  Computing risk measures for various distributions.
}
\usage{
VaR_t(alpha, mu=0, sigma=1, df=Inf)
ES_t(alpha, mu=0, sigma=1, df=Inf)
VaR_Par(alpha, theta)
ES_Par(alpha, theta)
}
\arguments{
  \item{alpha}{confidence level (in [0,1]).}
  \item{mu}{location parameter.}
  \item{sigma}{scale parameter, a positive number.}
  \item{df}{degrees of freedom, a positive number; choose \code{df=Inf}
    for the normal distribution.}
  \item{theta}{Pareto parameter, a positive number.}
}
\value{
  \code{VaR_t()}, \code{ES_t()} compute Value-at-Risk and expected
  shortfall for the \eqn{t} (or normal) distribution.

  \code{VaR_Par()}, \code{ES_Par()} compute Value-at-Risk and expected
  shortfall for the Pareto distribution.
}
\details{
  The distribution function of the (standard) Pareto distribution is given by
  \deqn{F(x)=1-(1+x)^{-\theta},\ x\ge 0,}{F(x)=1-(1+x)^{-theta}, x>=0,}
  where \eqn{\theta>0}{theta>0}.
}
\author{Marius Hofert}
\references{
   McNeil, A. J., Frey, R., and Embrechts, P. (2015).
   \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
   Princeton University Press.
}
\examples{
## Example 2.16 in McNeil, Frey, Embrechts (2015)
V <- 10000 # value of the portfolio today
sig <- 0.2/sqrt(250) # daily volatility (annualized volatility of 20%)
nu <- 4 # degrees of freedom for the t distribution
alpha <- seq(0.001, 0.999, length.out=256) # confidence levels
VaRnorm <- VaR_t(alpha, sigma=V*sig, df=Inf)
VaRt4 <- VaR_t(alpha, sigma=V*sig*sqrt((nu-2)/nu), df=nu)
ESnorm <- ES_t(alpha, sigma=V*sig, df=Inf)
ESt4 <- ES_t(alpha, sigma=V*sig*sqrt((nu-2)/nu), df=nu)
ran <- range(VaRnorm, VaRt4, ESnorm, ESt4)
plot(alpha, VaRnorm, type="l", ylim=ran, xlab=expression(alpha), ylab="")
lines(alpha, VaRt4, col="royalblue3")
lines(alpha, ESnorm, col="darkorange2")
lines(alpha, ESt4, col="maroon3")
legend("bottomright", bty="n", lty=rep(1,4), col=c("black",
       "royalblue3", "darkorange3", "maroon3"),
       legend=c(expression(VaR[alpha]~~"for normal model"),
                expression(VaR[alpha]~~"for "*t[4]*" model"),
                expression(ES[alpha]~~"for normal model"),
                expression(ES[alpha]~~"for "*t[4]*" model")))
}
\keyword{models}