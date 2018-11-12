\name{risk_measures}
\alias{VaR_np}
\alias{VaR_t}
\alias{VaR_GPD}
\alias{VaR_Par}
\alias{VaR_GPDtail}
\alias{ES_np}
\alias{ES_t}
\alias{ES_GPD}
\alias{ES_Par}
\alias{ES_GPDtail}
\alias{gVaR}
\alias{gEX}
\title{Risk Measures}
\description{
  Computing risk measures.
}
\usage{
## Value-at-risk
VaR_np(x, level, names = FALSE, type = 1, ...)
VaR_t(level, loc = 0, scale = 1, df = Inf)
VaR_GPD(level, shape, scale)
VaR_Par(level, shape, scale = 1)
VaR_GPDtail(level, threshold, p.exceed, shape, scale)

## Expected shortfall
ES_np(x, level, method = c(">", ">="), verbose = FALSE, ...)
ES_t(level, loc = 0, scale = 1, df = Inf)
ES_GPD(level, shape, scale)
ES_Par(level, shape, scale = 1)
ES_GPDtail(level, threshold, p.exceed, shape, scale)

## Multivariate geometric value-at-risk and expectiles
gVaR(x, level, start = colMeans(x),
     method = if(length(level) == 1) "Brent" else "Nelder-Mead", ...)
gEX(x, level, start = colMeans(x),
    method = if(length(level) == 1) "Brent" else "Nelder-Mead", ...)
}
\arguments{
  \item{x}{
    \describe{
      \item{\code{gVaR()}, \code{gEX()}}{\code{\link{matrix}} of
	(rowwise) multivariate losses.}
      \item{otherwise}{\code{\link{vector}} of losses.}
    }
  }
  \item{level}{
    \describe{
      \item{\code{gVaR()}, \code{gEX()}}{\code{\link{vector}} or
	\code{\link{matrix}} of	(rowwise) confidence levels \eqn{\alpha}{alpha}
	(all in \eqn{[0,1]}).}
      \item{otherwise}{confidence level \eqn{\alpha\in[0,1]}{alpha in
	  [0,1]}.}
    }
  }
  \item{names}{see \code{?\link{quantile}}.}
  \item{type}{see \code{?\link{quantile}}.}
  \item{loc}{location parameter \eqn{\mu}{mu}.}
  \item{shape}{
    \describe{
      \item{\code{VaR_GPD()}, \code{ES_GPD()}}{GPD shape parameter
	\eqn{\xi}{xi}, a real number.}
      \item{\code{VaR_Par()}, \code{ES_Par()}}{Pareto shape parameter
	\eqn{\theta}{theta}, a positive number.}
    }
  }
  \item{scale}{
    \describe{
      \item{\code{VaR_t()}, \code{ES_t()}}{\eqn{$t$}{t} scale parameter
	\eqn{\sigma}{sigma}, a positive number.}
      \item{\code{VaR_GPD()}, \code{ES_GPD()}}{GPD scale parameter
	\eqn{\beta}{beta}, a positive number.}
      \item{\code{VaR_Par()}, \code{ES_Par()}}{Pareto scale parameter
	\eqn{\kappa}{kappa}, a positive number.}
    }
  }
  \item{df}{degrees of freedom, a positive number; choose \code{df = Inf}
    for the normal distribution.}
  \item{threshold}{threhold \eqn{u}{u} (used to estimate the exceedance
    probability based on the data \code{x}).}
  \item{p.exceed}{exceedance probability; typically \code{mean(x > threshold)}
    for \code{x} being the data modeled with the peaks-over-threshold
    (POT) method.}
  \item{start}{\code{\link{vector}} of initial values for the underlying
    \code{\link{optim}()}.}
  \item{method}{
    \describe{
      \item{\code{ES_np()}}{\code{\link{character}} string indicating the method
	for computing expected shortfall.}
      \item{\code{gVaR()}, \code{gEX()}}{the optimization method passed to the
	underlying \code{\link{optim}()}.}
    }
  }
  \item{verbose}{\code{\link{logical}} indicating whether verbose
    output is given (in case the mean is computed over (too) few
    observations).}
  \item{\dots}{
    \describe{
      \item{\code{VaR_np()}}{additional arguments passed to the
	underlying \code{\link{quantile}()}.}
      \item{\code{ES_np()}}{additional arguments passed to
	\code{VaR_np()}.}
      \item{\code{gVaR()}, \code{gEX()}}{additional arguments passed to
	the underlying \code{\link{optim}()}.}
    }
  }
}
\value{
  \code{VaR_np()}, \code{ES_np()} estimate value-at-risk and expected
  shortfall non-parametrically. For the latter, the mean over all losses
  (strictly) beyond value-at-risk is computed. If \code{method = ">="},
  there is always at least one such loss, whereas if \code{method =
    ">"}, there might be no such loss, in which case \code{NaN} is returned.

  \code{VaR_t()}, \code{ES_t()} compute value-at-risk and expected
  shortfall for the \eqn{t} (or normal) distribution.

  \code{VaR_GPD()}, \code{ES_GPD()} compute value-at-risk and expected
  shortfall for the generalized Pareto distribution (GPD).

  \code{VaR_Par()}, \code{ES_Par()} compute value-at-risk and expected
  shortfall for the Pareto distribution.

  \code{gVaR()}, \code{gEX()} compute the multivariate geometric
  value-at-risk and expectiles suggested by Chaudhuri (1996) and
  Herrmann et al. (2018), respectively.
}
\details{
  The distribution function of the Pareto distribution is given by
  \deqn{F(x) = 1-(\kappa/(\kappa+x))^{\theta},\ x\ge 0,}{F(x) = 1-(kappa/(kappa+x))^{theta}, x
  >= 0,} where \eqn{\theta > 0}{theta > 0}, \eqn{\kappa > 0}{kappa > 0}.
}
\author{Marius Hofert}
\references{
   McNeil, A. J., Frey, R. and Embrechts, P. (2015).
   \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
   Princeton University Press.

   Chaudhuri, P. (1996).
   On a geometric notion of quantiles for multivariate data.
   \emph{Journal of the American Statistical Assosiation} 91(434),
   862--872.

   Herrmann, K., Hofert, M. and Mailhot, M. (2018).
   Multivariate geometric expectiles.
   \emph{Scandinavian Actuarial Journal}, 2018(7), 629--659.
}
\examples{
### 1 Univariate measures ######################################################

## Generate some losses and (non-parametrically) estimate VaR_alpha and ES_alpha
set.seed(271)
L <- rlnorm(1000, meanlog = -1, sdlog = 2) # L ~ LN(mu, sig^2)
## Note: - meanlog = mean(log(L)) = mu, sdlog = sd(log(L)) = sig
##       - E(L) = exp(mu + (sig^2)/2), var(L) = (exp(sig^2)-1)*exp(2*mu + sig^2)
##         To obtain a sample with E(L) = a and var(L) = b, use:
##         mu = log(a)-log(1+b/a^2)/2 and sig = sqrt(log(1+b/a^2))
VaR_np(L, level = 0.99)
ES_np(L,  level = 0.99)

## Example 2.16 in McNeil, Frey, Embrechts (2015)
V <- 10000 # value of the portfolio today
sig <- 0.2/sqrt(250) # daily volatility (annualized volatility of 20%)
nu <- 4 # degrees of freedom for the t distribution
alpha <- seq(0.001, 0.999, length.out = 256) # confidence levels
VaRnorm <- VaR_t(alpha, scale = V*sig, df = Inf)
VaRt4 <- VaR_t(alpha, scale = V*sig*sqrt((nu-2)/nu), df = nu)
ESnorm <- ES_t(alpha, scale = V*sig, df = Inf)
ESt4 <- ES_t(alpha, scale = V*sig*sqrt((nu-2)/nu), df = nu)
ran <- range(VaRnorm, VaRt4, ESnorm, ESt4)
plot(alpha, VaRnorm, type = "l", ylim = ran, xlab = expression(alpha), ylab = "")
lines(alpha, VaRt4, col = "royalblue3")
lines(alpha, ESnorm, col = "darkorange2")
lines(alpha, ESt4, col = "maroon3")
legend("bottomright", bty = "n", lty = rep(1,4), col = c("black",
       "royalblue3", "darkorange3", "maroon3"),
       legend = c(expression(VaR[alpha]~~"for normal model"),
                  expression(VaR[alpha]~~"for "*t[4]*" model"),
                  expression(ES[alpha]~~"for normal model"),
                  expression(ES[alpha]~~"for "*t[4]*" model")))


### 2 Multivariate measures ####################################################

## Setup
library(copula)
n <- 1e4 # MC sample size
nu <- 3 # degrees of freedom
th <- iTau(tCopula(df = nu), tau = 0.5) # correlation parameter
cop <- tCopula(param = th, df = nu) # t copula
set.seed(271) # for reproducibility
U <- rCopula(n, cop = cop) # copula sample
theta <- c(2.5, 4) # marginal Pareto parameters
stopifnot(theta > 2) # need finite 2nd moments
X <- sapply(1:2, function(j) qPar(U[,j], shape = theta[j])) # generate X
N <- 17 # number of angles (rather small here because of run time)
phi <- seq(0, 2*pi, length.out = N) # angles
r <- 0.98 # radius
alpha <- r * cbind(alpha1 = cos(phi), alpha2 = sin(phi)) # vector of confidence levels

## Compute geometric value-at-risk
system.time(res <- gVaR(X, level = alpha))
gvar <- t(sapply(seq_len(nrow(alpha)), function(i) {
    x <- res[[i]]
    if(x[["convergence"]] != 0) # 0 = 'converged'
        warning("No convergence for alpha = (", alpha[i,1], ", ", alpha[i,2],
                ") (row ", i, ")")
    x[["par"]]
})) # (N, 2)-matrix

## Compute geometric expectiles
system.time(res <- gEX(X, level = alpha))
gex <- t(sapply(seq_len(nrow(alpha)), function(i) {
    x <- res[[i]]
    if(x[["convergence"]] != 0) # 0 = 'converged'
        warning("No convergence for alpha = (", alpha[i,1], ", ", alpha[i,2],
                ") (row ", i, ")")
    x[["par"]]
})) # (N, 2)-matrix

## Plot geometric VaR and geometric expectiles
plot(gvar, type = "b", xlab = "Component 1 of geometric VaRs and expectiles",
     ylab = "Component 2 of geometric VaRs and expectiles",
     main = "Multivariate geometric VaRs and expectiles")
lines(gex, type = "b", col = "royalblue3")
legend("bottomleft", lty = 1, bty = "n", col = c("black", "royalblue3"),
       legend = c("geom. VaR", "geom. expectile"))
lab <- substitute("MC sample size n ="~n.*","~t[nu.]~"copula with Par("*th1*
                  ") and Par("*th2*") margins",
                  list(n. = n, nu. = nu, th1 = theta[1], th2 = theta[2]))
mtext(lab, side = 4, line = 1, adj = 0)
}
\keyword{models}