### Fitting and predicting VaR based on an ARMA(p,q)-GARCH(p.,q.) process

require(rugarch)


### 1) Simulate (log-return) data (X_t) from an ARMA(1,1)-GARCH(1,1) process ###

## ARMA(1,1)-GARCH(1,1) specification
fixed.p <- list(mu=0, ar1=0.5, ma1=0.3, omega=4, alpha1=0.4, beta1=0.2)
armaOrder <- c(1,1)
garchOrder <- c(1,1)
varModel <- list(model="sGARCH", garchOrder=garchOrder)
distr <- "norm" # or use "std"
spec <- ugarchspec(varModel, mean.model=list(armaOrder=armaOrder),
                   fixed.pars=fixed.p, distribution.model=distr)

## simulate one path (our example data X here)
n <- 1000 # sample size
set.seed(271) # reproducibility
x <- ugarchpath(spec, n.sim=n, m.sim=1) # n = length of simulated path; m = number of paths
if(FALSE)
    str(x, max.level=2) # components path, model, seed
X <- x@path$seriesSim # simulated process
sig <- x@path$sigmaSim # standard deviation
eps <- x@path$residSim # residuals

## sanity check
if(FALSE) {
    plot(X, type="l", xlab="t", ylab=expression("Simulated process"~X[t]))
    plot(sig, type="l", xlab="t", ylab=expression("Conditional standard deviation"~sigma[t]))
    plot(eps, type="l", xlab="t", ylab=expression("Residuals"~epsilon[t]))
}


### 2) Fit an ARMA-GARCH model to the (simulated) data #########################

## fit an ARMA(p,q)-GARCH(p.,q.) process to X
## (with the correct, known orders (here))
mod <- list(model="sGARCH", garchOrder=garchOrder)
spec <- ugarchspec(mod, mean.model=list(armaOrder=armaOrder),
                   distribution.model=distr)
fit <- ugarchfit(spec, data=X) # components fit, model
mu <- fit@fit$fitted.values # \mu_t = \hat{X}_t (as E[Z]=0)

## sanity check
if(FALSE) {
    plot(X, type="l", xlab="t",
         ylab=expression("Data"~X[t]~"and fitted values"~hat(mu)[t])) # data
    lines(mu, col=adjustcolor("blue", alpha.f=0.4)) # fitted values
    plot(fit@fit$residuals, type="l", xlab="t", ylab=expression(epsilon[t])) # check residuals
    Z <- fit@fit$z
    qqplot(qnorm(ppoints(length(Z))), Z, xlab="N(0,1) quantiles",
           ylab="Z quantiles") # check distribution of Z
    qqline(Z, distribution=qnorm)
}


### 3) Calculate the VaR time series ###########################################

## compute VaR estimate
alpha <- 0.99
VaR <- as.numeric(quantile(fit, probs=alpha)) # a vector (since fit is a rugarch object)
## => computes \hat{mu}_t+\hat{sigma}_t*q_Z(alpha)

## sanity check
if(FALSE) {
    sig <- fit@fit$sigma # extract sigma_t
    mu <- fit@fit$fitted.values # extract mu_t; note: E[Z] = 0
    VaR. <- mu+sig*qnorm(alpha) # VaR_alpha computed by hand
    stopifnot(VaR. == VaR)
    ## conclusion: quantile("rugarch object", probs=alpha) provides VaR_alpha
}


### 4) Backtesting via randomness check ########################################

btest <- VaRTest(alpha, actual=X, VaR=VaR, conf.level=0.95) # backtest for VaR_0.95
if(FALSE) {
    btest$expected.exceed
    btest$actual.exceed
    btest$uc.Decision # unconditional test decision (note: cc.Decision is NA here)
}


### 5) Predict VaR based on fitted model #######################################

m <- ceiling(n/10) # number of steps to forecast; => roll m-1 times with frequency 1
fspec <- getspec(fit) # specification of the fitted process
setfixed(fspec) <- as.list(coef(fit))
pred <- ugarchforecast(fspec, data=X, n.ahead=1, n.roll=m-1, out.sample=m) # predict from the fitted process
mu.predict <- pred@forecast$seriesFor # extract predicted mu_t; note: E[Z] = 0
VaR.predict <- as.numeric(quantile(pred, probs=alpha))

## sanity check
if(FALSE) {
    sig <- pred@forecast$sigmaFor # extract predicted sigma_t
    VaR. <- mu+sig*qnorm(alpha) # VaR_alpha computed by hand
    stopifnot(VaR. == VaR.predict)
    ## conclusion: quantile("rugarch object", probs=alpha) provides VaR_alpha
}


### 6) Simulate future trajectories of (X_t) and compute corresponding VaRs ####

## simulated paths
B <- 1000 # bootstrap B paths
X.boot <- ugarchpath(fspec, n.sim=m, m.sim=B) # simulate future paths
if(FALSE)
    str(X.boot, max.level=2) # components path, model, seed

## estimate VaR for each simulated path
## note: quantile() can't be used here => we have to construct VaR 'by hand'
if(FALSE)
    str(X.boot@path) # each series is now an (m, B) matrix (each row is one path)
sig.t.boot <- X.boot@path$sigmaSim # extract sigma_t
mu.t.boot <- X.boot@path$seriesSim # extract mu_t; note: E[Z] = 0
VaR.boot <- mu.t.boot+sig.t.boot*qnorm(alpha) # bootstrapped VaR_alpha computed by hand; (m, B) matrix

## compute bootstrapped two-sided 95%-confidence intervals for VaR_alpha
VaR.CI <- apply(VaR.boot, 1, function(x) quantile(x, probs=c(0.025, 0.975)))


### 7) Plot ####################################################################

## plot
yran <- range(X, mu, VaR, mu.predict, VaR.predict, VaR.CI)
xran <- c(1, length(X) + m)
## past data (simulated time series X_t)
plot(X, type="l", xlim=xran, ylim=yran, ylab="",
     main="Data, fitted ARMA-GARCH process, VaR, VaR predictions and VaR CIs",
     sub=paste0("Expected exceedances: ", btest$expected.exceed, "   Actual exceedances: ",
                btest$actual.exceed, "   Test decision: ", btest$uc.Decision))
lines(mu, col=adjustcolor("darkblue", alpha.f=0.4)) # \hat{\mu}_t
lines(VaR, col="darkred") # estimated VaR_alpha
## predictions
t. <- length(X) + seq_len(m) # future time points
lines(t., mu.predict, col="blue") # predicted process X_t (or \mu_t)
lines(t., VaR.predict, col="red") # predicted VaR_alpha
lines(t., VaR.CI[1,], col="orange") # lower 95%-CI for VaR_alpha
lines(t., VaR.CI[2,], col="orange") # upper 95%-CI for VaR_alpha
## legend
legend("topright", bty="n", lty=rep(1, 6), lwd=1.6,
       col=c("black", "darkblue", "blue", "darkred", "red", "orange"),
       legend=c(expression(X[t]), expression(hat(mu)[t]),
                expression("predicted"~mu[t]~"(or"~X[t]*")"),
                substitute(widehat(VaR)[a], list(a=alpha)),
                substitute("predicted"~VaR[a], list(a=alpha)),
                substitute("95%-CI for"~VaR[a], list(a=alpha))))
