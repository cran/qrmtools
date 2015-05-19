### Worst Value-at-Risk under known margins ####################################

require(qrmtools)
doPDF <- !dev.interactive(orNone=TRUE)


### 1) Homogeneous case (all margins equal) ####################################

## Setup
alpha <- 0.99 # confidence level
d <- 8 # dimension (affects file names below)
th <- 2 # Pareto parameter (affects file names below)
qF <- function(p) qPar(p, theta=th) # Pareto quantile function
pF <- function(q) pPar(q, theta=th) # Pareto distribution function


### 1.1) Checks for method="dual" ##############################################

## Investigating h(s,t) = D(s,t) - ..., the function for the inner
## root-finding to compute D(s); see dual_bound()
s <- c(1, 5, 10, 100, 500, 1000)
t <- sapply(seq_along(s), function(i) {
    res <- exp(seq(log(1e-3), log(s[i]/d), length.out=257))
    res[length(res)] <- s[i]/d # to avoid numerical issues (t > s/d)
    res
})
f <- sapply(seq_along(s), function(i)
            sapply(t[,i], function(t.)
                   qrmtools:::dual_bound_2(s[i], t=t., d=d, pF=pF) -
                   qrmtools:::dual_bound_2_deriv_term(s[i], t=t., d=d, pF=pF)))
palette <- colorRampPalette(c("red", "orange", "blue"), space="Lab")
cols <- palette(6)
if(doPDF)
    pdf(file=(file <- paste0("fig_worst_VaR_",alpha,"_hom_dual_h_Par=",th,"_d=",d,".pdf")),
        width=6, height=6)
par(pty="s")
plot(t[,1], f[,1], type="l", log="x", xlim=range(t), ylim=range(f), col=cols[1],
     xlab="t", ylab=expression("h(s,t) for d = 8 and F being Par(2)"))
lines(t[,2], f[,2], col=cols[2])
lines(t[,3], f[,3], col=cols[3])
lines(t[,4], f[,4], col=cols[4])
lines(t[,5], f[,5], col=cols[5])
lines(t[,6], f[,6], col=cols[6])
abline(h=0, lty=2)
legend("topright", inset=0.02, lty=rep(1,6), col=cols,
       bty="n", legend=as.expression(lapply(1:6,
           function(i) substitute(s==s., list(s.=s[i])))), y.intersp=1.2)
if(doPDF) dev.off.pdf(file=file)
## Conclusion: As we know, h(s, s/d) = 0. We also see that s has to be
##             sufficiently large in order to find a root h(s, t) = 0 for t < s/d

## Plot dual bound D(s) for various theta (concerns the outer root-finding)
theta <- c(0.5, 1, 2, 4) # theta values
s <- seq(48, 2000, length.out=257) # s values
D <- sapply(theta, function(th)
            sapply(s, function(s.)
                   dual_bound(s., d=8, pF=function(q) pPar(q, theta=th)))) # (s, theta) matrix
if(doPDF)
    pdf(file=(file <- paste0("fig_worst_VaR_",alpha,"_hom_dual_D_s_Par=",
                             paste0(theta, collapse="_"),"_d=",d,".pdf")),
        width=6, height=6)
par(pty="s")
plot(s, D[,1], type="l", ylim=range(D),
     ylab=expression("Dual bound D(s) for d = 8 and F being Par("*theta*")"))
lines(s, D[,2], col="blue")
lines(s, D[,3], col="orange")
lines(s, D[,4], col="red")
legend("topright", inset=0.02, lty=rep(1,4), y.intersp=1.2,
       col=c("black", "blue", "orange", "red"),
       bty="n", legend=as.expression(lapply(1:4,
           function(i) substitute(theta==i, list(i=theta[i])))))
if(doPDF) dev.off.pdf(file=file)


### 1.2) Checks for method="Wang"/"Wang.Par" ###################################

### Check *with* numerical integration

## Check Wang_Ibar()
c <- seq(0, (1-alpha)/d, length.out=129) # initial interval for root finding
yc <- sapply(c, function(c.) qrmtools:::Wang_Ibar(c., alpha=alpha, d=d, qF=qF))
par(mar=c(5.1, 6.1, 4.1, 2.1)) # more space for the y-axis label
plot(c, yc, type="l", xlab="c (in initial interval)",
     ylab=expression(bar(I)(c)))

## Check Wang_h_aux()
yc. <- qrmtools:::Wang_h_aux(c=c, alpha=alpha, d=d, qF=qF)
plot(c, yc., type="l", xlab="c (in initial interval)",
     ylab=expression(frac(d-1,d)~{F^{-1}}(a[c])+frac(1,d)~{F^{-1}}(b[c])))

## Check objective function h(c) (Wang_h()) for the chosen theta
yc.. <- sapply(c, function(c.) qrmtools:::Wang_h(c., alpha=alpha, d=d, qF=qF))
if(doPDF)
    pdf(file=(file <- paste0("fig_worst_VaR_",alpha,"_hom_Wang_h_Par=",th,"_d=",d,".pdf")),
        width=6.5, height=6.5)
par(pty="s")
plot(c, yc.., type="l", xlab="c (in initial interval)",
     ylab=expression("h(c) for"~~alpha~"= 0.99, d = 8 and F being Par(2)"))
abline(h=0, lty=2)
if(doPDF) dev.off.pdf(file=file)

## Check endpoints of objective function for root-finding
sapply(c(0, (1-alpha)/d), function(c.)
       qrmtools:::Wang_h(c., alpha=alpha, d=d, qF=qF)) # -Inf, 0
## -Inf is not a problem for root finding (for theta > 1; for theta <= 1 it is NaN),
## but the 0 at the right endpoint is. We take care of this by adjusting
## f.upper in the root finding; see worst_VaR_hom()


### Check *without* numerical integration

## Check objective function h(c) (Wang_h()) for theta = 1/2 (remains so flat for d=100)
yc... <- sapply(c, function(c.) qrmtools:::Wang_h(c., alpha=alpha, d=d, method="Wang.Par", theta=1/2))
if(doPDF)
    pdf(file=(file <- paste0("fig_worst_VaR_",alpha,"_hom_Wang_h_Par=",0.5,"_d=",d,".pdf")),
        width=6.5, height=6.5)
par(pty="s")
plot(c, yc..., type="l", xlab="c (in initial interval)",
     ylab=expression("h(c) for"~~alpha~"= 0.99, d = 8 and F being Par(1/2)"))
abline(h=0, lty=2)
if(doPDF) dev.off.pdf(file=file)

## Compute worst VaR_alpha for d=8 and various Par(theta)
## Note: for theta in (0, 1], c_l has to be > 0
theta. <- c(0.5, 1, 2, 4) # theta values
alpha. <- rev(1-2^seq(-10, -0.001, length.out=128)) # alpha values (in (0,1); concentrated near 1)
d <- 8 # or d=100
worst.VaR.Wang  <- sapply(theta., function(th)
                          sapply(alpha., function(a) {
                              worst_VaR_hom(a, d=d, interval=c(if(th <= 1) .Machine$double.eps else 0,
                                                               (1-a)/d), # as worst_VaR_hom()
                                            method="Wang.Par", theta=th) #, tol=1e-100) (no improvement)
                                            ## numerical integration fails here:
                                            ## method="Wang", qF=function(p) qPar(p, theta=th))
                          })) # (alpha, theta) matrix

## Plot
if(doPDF)
    pdf(file=(file <- paste0("fig_worst_VaR_hom_Wang_Par_d=",d,".pdf")),
        width=6.5, height=6.5)
par(pty="s")
plot(1-alpha., worst.VaR.Wang[,1], type="l", log="xy", ylim=range(worst.VaR.Wang),
     xlab=expression(1-alpha),
     ylab=substitute(bar(VaR)[alpha]*group("(",L^{"+"},")")~
     "for d ="~d.~"and F being Par("*theta*")", list(d.=d)))
lines(1-alpha., worst.VaR.Wang[,2], col="blue")
lines(1-alpha., worst.VaR.Wang[,3], col="orange")
lines(1-alpha., worst.VaR.Wang[,4], col="red")
legend("topright", inset=0.02, y.intersp=1.2, lty=rep(1,4),
       col=c("black", "blue", "orange", "red"),
       bty="n", legend=as.expression(lapply(1:4,
           function(i) substitute(theta==i, list(i=theta.[i])))))
if(doPDF) dev.off.pdf(file=file)


### 2) Compare 'crude', "Wang" (numerical), "Wang.Par", "dual", 'RA' ###########

### 2.1) Addressing numerical issues for theta <= 1 ############################

## Setup (as before)
alpha <- 0.99 # confidence level
d <- 8 # dimension

## Why we can't use interval[1]=0 if theta in (0,1]
interval <- c(0, (1-alpha)/d)
c <- interval[1]
method <- "Wang.Par"
th <- 0.99
qrmtools:::Wang_h(c, alpha=alpha, d=d, method=method, theta=th) # NaN => uniroot() fails
qrmtools:::Wang_Ibar(c, alpha=alpha, d=d, method=method, theta=th) # Inf
qrmtools:::Wang_h_aux(c, alpha=alpha, d=d, method=method, theta=th) # Inf
## note: Wang_h() is NaN for c <= 1e-17
qrmtools:::Wang_h(.Machine$double.eps, alpha=alpha, d=d, method=method, theta=th) # okay

## Check the default interval[1]=.Machine$double.eps if theta in (0,1]
n.th <- 128 # number of thetas
th <- 2^seq(-4, 2, length.out=n.th) # thetas
VaR  <- sapply(th, function(th.) worst_VaR_hom(alpha, d=d, method="Wang.Par", theta=th.))
plot(th, VaR,  type="l", log="y")
## VaR. <- sapply(th, function(th.) worst_VaR_hom(alpha, d=d, method="Wang.Par", theta=th.,
##                                                tol=.Machine$double.eps))
## plot(th, VaR., type="l", log="y")
## => no improvement around theta = 0.5


### 2.2) Graphical comparison ##################################################

## Setup
alpha <- 0.99 # confidence level
d <- 8 # dimension
n.th <- 64 # number of thetas
th.low <- 1.05 # or 0.05 (to see how much the methods differ for theta in (0,1])
th.up <- 5
th <- seq(th.low, th.up, length.out=n.th) # thetas
qFs <- lapply(th, function(th.) {th.; function(p) qPar(p, theta=th.)}) # n.th-vector of Pareto quantile functions
pFs <- lapply(th, function(th.) {th.; function(q) pPar(q, theta=th.)}) # n.th-vector of Pareto dfs
N <- 1e4 # number of discretization points for RA(); N=1e5 does not improve the situation

## Compute values
res <- matrix(, nrow=n.th, ncol=8)
colnames(res) <- c("crude.low", "crude.up", "Wang", "Wang.Par",
                   "Wang.Par.uniroot.tol", "dual", "RA.low", "RA.up")
## ~= 1--2min
for(i in seq_len(n.th)) {
    ## crude bounds (also used as initial interval for method "dual" below)
    I <- crude_VaR_bounds(alpha, d=d, qF=qFs[[i]])
    res[i,"crude.low"] <- I[1]
    res[i,"crude.up"] <- I[2]
    ## method "Wang" (numerical integration critical for theta > 1 small)
    Wang.num.res <- tryCatch(worst_VaR_hom(alpha, d=d, qF=qFs[[i]]), error=function(e) e)
    if(is(Wang.num.res, "simpleError")) {
        warning("there was an error: ", conditionMessage(Wang.num.res), " (will use NA as result)")
        Wang.num.res <- NA
    }
    res[i,"Wang"] <- Wang.num.res
    ## method "Wang.Par"
    res[i,"Wang.Par"] <- worst_VaR_hom(alpha, d=d, method="Wang.Par", theta=th[i])
    ## method "Wang.Par" (with uniroot()'s default tolerance)
    res[i,"Wang.Par.uniroot.tol"] <- worst_VaR_hom(alpha, d=d, method="Wang.Par", theta=th[i],
                                                   tol=.Machine$double.eps^0.25)
    ## method "dual"
    res[i,"dual"] <- worst_VaR_hom(alpha, d=d, method="dual", interval=I,
                                   pF=pFs[[i]])
    ## Rearrangement Algorithm
    set.seed(271) # use the same sampling for each theta
    RA. <- RA(alpha, d=d, qF=qFs[[i]], N=N, abs.err=0.001)
    res[i,"RA.low"] <- RA.$bounds[1]
    res[i,"RA.up"]  <- RA.$bounds[2]
}

## Plot
res. <- res[,c(-1,-2)] # omit crude bounds, too crude
res. <- res.[,1:5]/res.[,6] # standardize w.r.t. RA.up
if(doPDF) {
    file <- paste0("fig_worst_VaR_",alpha,"_hom_comparison_d=",d,"_N=",N,"_theta_low=",th.low,".pdf")
    pdf(file=(file), width=6.5, height=6.5)
}
par(pty="s")
plot(th, res.[,"Wang"], type="l", ylim=range(res., na.rm=TRUE),
     xlab=expression(theta), ylab=expression(bar(VaR)[0.99]*group("(",L^{"+"},")")~
     "standardized by the upper RA bound for d = 8 and F being Par("*theta*")"),
     col="gray80", lty=2, lwd=5)
lines(th, res.[,"Wang.Par"], col="gray50", lty=2, lwd=3) # ~= res.[,"Wang"]
lines(th, res.[,"dual"], col="gray20", lty=2, lwd=1) # ~= res.[,"Wang"]
lines(th, res.[,"RA.low"], col="black")
lines(th, res.[,"Wang.Par.uniroot.tol"], col="red")
legend("topright", inset=0.02, y.intersp=1.2, bty="n",
       col=c("gray80", "gray50", "gray20", "black", "red"),
       lty=c(2,2,2,1,1), lwd=c(5,3,1,1,1),
       legend=c("Wang (num. int.)", "Wang", "Dual bound",
                "lower RA bound", "Wang (uniroot() tol.)"))
if(doPDF) dev.off.pdf(file=file)
