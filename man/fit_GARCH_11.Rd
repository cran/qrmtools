\name{fit_GARCH_11}
\alias{fit_GARCH_11}
\title{Fast(er) and Numerically More Robust Fitting of GARCH(1,1) Processes}
\description{
  Fast(er) and numerically more robust fitting of GARCH(1,1) processes
  according to Zumbach (2000).
}
\usage{
fit_GARCH_11(x, init = NULL, sig2 = mean(x^2), delta = 1,
             distr = c("norm", "st"), control = list(), ...)
}
\arguments{
  \item{x}{vector of length \eqn{n} containing the data (typically
    log-returns) to be fitted a GARCH(1,1) to.}
  \item{init}{vector of length 2 giving the initial values for the
    likelihood fitting. Note that these are initial values for
    \eqn{z_{corr}}{z[corr]} and \eqn{z_{ema}}{z[ema]}
    as in Zumbach (2000).}
  \item{sig2}{annualized variance (third parameter of the
    reparameterization according to Zumbach (2000)).}
  \item{delta}{unit of time (defaults to 1 meaning daily data;
    for yearly data, use 250).}
  \item{distr}{character string specifying the innovation distribution
    (\code{"norm"} for N(0,1) or \code{"st"} for a standardized \eqn{t}
    distribution).}
  \item{control}{see \code{?\link{optim}()}.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{optim}()}.}
}
\value{
  A list with components
  \describe{
    \item{coef}{estimated coefficients \eqn{\alpha_0}{alpha[0]},
      \eqn{\alpha_1}{alpha[1]}, \eqn{\beta_1}{beta[1]} and, if
      \code{distr == "st"} the estimated degrees of freedom.}
     \item{logLik}{maximized log-likelihood.}
     \item{counts}{number of calls to the objective function; see
       \code{?\link{optim}}.}
     \item{convergence}{convergence code ('0' indicates successful
       completion; see \code{?\link{optim}}.}
     \item{message}{see \code{?\link{optim}}.}
     \item{sig.t}{vector of length \eqn{n} giving the conditional
       volatility.}
     \item{Z.t}{vector of length \eqn{n} giving the standardized
       residuals.}
  }
}
\author{Marius Hofert (based on an early version by Marcel Braeutigam)}
\references{
  Zumbach, G. (2000). The pitfalls in fitting GARCH (1,1) processes.
  \emph{Advances in Quantitative Asset Management} \bold{1}, 179--200.
}
\seealso{
  \code{\link{fit_ARMA_GARCH}()} based on \pkg{rugarch}.
}
\examples{
### Example 1: N(0,1) innovations ##############################################

## Generate data from a GARCH(1,1) with N(0,1) innovations
library(rugarch)
uspec <- ugarchspec(variance.model = list(model = "sGARCH",
                                          garchOrder = c(1, 1)),
                    distribution.model = "norm",
                    mean.model = list(armaOrder = c(0, 0)),
                    fixed.pars = list(mu = 0,
                                      omega = 0.1, # alpha_0
                                      alpha1 = 0.2, # alpha_1
                                      beta1 = 0.3)) # beta_1
X <- ugarchpath(uspec, n.sim = 1e4, rseed = 271) # sample (set.seed() fails!)
X.t <- as.numeric(X@path$seriesSim) # actual path (X_t)

## Fitting via ugarchfit()
uspec. <- ugarchspec(variance.model = list(model = "sGARCH",
                                           garchOrder = c(1, 1)),
                     distribution.model = "norm",
                     mean.model = list(armaOrder = c(0, 0)))
fit <- ugarchfit(uspec., data = X.t)
coef(fit) # fitted mu, alpha_0, alpha_1, beta_1
Z <- fit@fit$z # standardized residuals
stopifnot(all.equal(mean(Z), 0, tol = 1e-2),
          all.equal(var(Z),  1, tol = 1e-3))

## Fitting via fit_GARCH_11()
fit. <- fit_GARCH_11(X.t)
fit.$coef # fitted alpha_0, alpha_1, beta_1
Z. <- fit.$Z.t # standardized residuals
stopifnot(all.equal(mean(Z.), 0, tol = 5e-3),
          all.equal(var(Z.),  1, tol = 1e-3))

## Compare
stopifnot(all.equal(fit.$coef, coef(fit)[c("omega", "alpha1", "beta1")],
                    tol = 5e-3, check.attributes = FALSE)) # fitted coefficients
summary(Z. - Z) # standardized residuals


### Example 2: t_nu(0, (nu-2)/nu) innovations ##################################

## Generate data from a GARCH(1,1) with t_nu(0, (nu-2)/nu) innovations
uspec <- ugarchspec(variance.model = list(model = "sGARCH",
                                          garchOrder = c(1, 1)),
                    distribution.model = "std",
                    mean.model = list(armaOrder = c(0, 0)),
                    fixed.pars = list(mu = 0,
                                      omega = 0.1, # alpha_0
                                      alpha1 = 0.2, # alpha_1
                                      beta1 = 0.3, # beta_1
                                      shape = 4)) # nu
X <- ugarchpath(uspec, n.sim = 1e4, rseed = 271) # sample (set.seed() fails!)
X.t <- as.numeric(X@path$seriesSim) # actual path (X_t)

## Fitting via ugarchfit()
uspec. <- ugarchspec(variance.model = list(model = "sGARCH",
                                           garchOrder = c(1, 1)),
                     distribution.model = "std",
                     mean.model = list(armaOrder = c(0, 0)))
fit <- ugarchfit(uspec., data = X.t)
coef(fit) # fitted mu, alpha_0, alpha_1, beta_1, nu
Z <- fit@fit$z # standardized residuals
stopifnot(all.equal(mean(Z), 0, tol = 1e-2),
          all.equal(var(Z),  1, tol = 5e-2))

## Fitting via fit_GARCH_11()
fit. <- fit_GARCH_11(X.t, distr = "st")
c(fit.$coef, fit.$df) # fitted alpha_0, alpha_1, beta_1, nu
Z. <- fit.$Z.t # standardized residuals
stopifnot(all.equal(mean(Z.), 0, tol = 2e-2),
          all.equal(var(Z.),  1, tol = 2e-2))

## Compare
fit.coef <- coef(fit)[c("omega", "alpha1", "beta1", "shape")]
fit..coef <- c(fit.$coef, fit.$df)
stopifnot(all.equal(fit.coef, fit..coef, tol = 7e-2, check.attributes = FALSE))
summary(Z. - Z) # standardized residuals
}
\keyword{ts}