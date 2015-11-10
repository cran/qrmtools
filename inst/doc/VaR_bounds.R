## ---- message=FALSE------------------------------------------------------
require(qrmtools)
require(QRM)
require(copula)
require(combinat)
require(sfsmisc)
doPDF <- FALSE

## ------------------------------------------------------------------------
qF <- function(p, th=2) qPar(p, theta=th) # Pareto quantile function
pF <- function(q, th=2) pPar(q, theta=th) # Pareto distribution function
dim <- 8 # variable dimension (we use 8 or 100 here)

## ------------------------------------------------------------------------
d <- 8 # dimension
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
palette <- colorRampPalette(c("maroon3", "darkorange2", "royalblue3"), space="Lab")
cols <- palette(6)
if(doPDF)
    pdf(file=(file <- paste0("fig_worst_VaR_hom_dual_h_Par=2_d=",d,".pdf")),
        width=6, height=6)
plot(t[,1], f[,1], type="l", log="x", xlim=range(t), ylim=range(f), col=cols[1],
     xlab="t", ylab=expression("h(s,t) for d = 8 and F being Par(2)"))
lines(t[,2], f[,2], col=cols[2])
lines(t[,3], f[,3], col=cols[3])
lines(t[,4], f[,4], col=cols[4])
lines(t[,5], f[,5], col=cols[5])
lines(t[,6], f[,6], col=cols[6])
abline(h=0, lty=2)
legend("topright", lty=rep(1,6), col=cols,
       bty="n", legend=as.expression(lapply(1:6,
           function(i) substitute(s==s., list(s.=s[i])))))
if(doPDF) dev.off()

## ------------------------------------------------------------------------
theta <- c(0.5, 1, 2, 4) # theta values
s <- seq(48, 2000, length.out=257) # s values
D <- sapply(theta, function(th)
            sapply(s, function(s.)
                   dual_bound(s., d=d, pF=function(q) pPar(q, theta=th)))) # (s, theta) matrix
if(doPDF)
    pdf(file=(file <- paste0("fig_worst_VaR_hom_dual_D_s_Par=",
                             paste0(theta, collapse="_"),"_d=",d,".pdf")),
        width=6, height=6)
plot(s, D[,1], type="l", ylim=range(D), col="maroon3",
     ylab=substitute("Dual bound D(s) for d ="~d.~"Par("*theta*") margins", list(d.=d)))
lines(s, D[,2], col="darkorange2")
lines(s, D[,3], col="royalblue3")
lines(s, D[,4], col="black")
legend("topright", lty=rep(1,4),
       col=c("maroon3", "darkorange2", "royalblue3", "black"),
       bty="n", legend=as.expression(lapply(1:4,
           function(j) substitute(theta==j, list(j=theta[j])))))
if(doPDF) dev.off()

## ------------------------------------------------------------------------
d <- 8 # dimension
alpha <- 0.99 # confidence level
c <- seq(0, (1-alpha)/d, length.out=129) # domain of h
h.aux <- qrmtools:::Wang_h_aux(c, alpha=alpha, d=d, qF=qF)
par(mar=c(5, 4+1, 4, 2) + 0.1) # increase space (for y axis label)
plot(c, h.aux, type="l", xlab="c (in initial interval)",
     ylab=expression(frac(d-1,d)~{F^{-1}}(a[c])+frac(1,d)~{F^{-1}}(b[c])))

## ------------------------------------------------------------------------
h <- sapply(c, function(c.) qrmtools:::Wang_h(c., alpha=alpha, d=d, qF=qF))
if(doPDF)
    pdf(file=(file <- paste0("fig_worst_VaR_",alpha,"_hom_Wang_h_Par=2_d=",d,"_num.pdf")),
        width=6, height=6)
plot(c, h, type="l", xlab="c (in initial interval)",
     ylab=substitute("h(c) for"~~alpha~"= 0.99 and d ="~d.~"Par(2) margins",
                     list(d.=d)))
abline(h=0, lty=2)
if(doPDF) dev.off()

## ------------------------------------------------------------------------
sapply(c(0, (1-alpha)/d), function(c.)
       qrmtools:::Wang_h(c., alpha=alpha, d=d, qF=qF)) # -Inf, 0

## ------------------------------------------------------------------------
method <- "Wang.Par" # this also holds for (the numerical) method="Wang"
th <- 0.99
qrmtools:::Wang_h(0, alpha=alpha, d=d, method=method, theta=th) # NaN => uniroot() fails
## Note: Wang_h() is actually already NaN for c <= 1e-17
qrmtools:::Wang_h_aux(0, alpha=alpha, d=d, method=method, theta=th) # Inf

## ------------------------------------------------------------------------
d <- dim # dimension
alpha <- 0.99 # confidence level
theta <- c(0.1, 0.5, 1, 5, 10, 50) # theta values
palette <- colorRampPalette(c("darkorange2", "maroon3", "royalblue3", "black"), space="Lab")
cols <- palette(length(theta))
c <- seq(0, (1-alpha)/d, length.out=2^13+1)
## => They all go further down to 0 if length.out is increased.
##    Smaller theta thus corresponds to a larger derivative in the root
##    Root-finding thus requires higher precision for smaller theta
h <- matrix(, nrow=length(c), ncol=length(theta))
for(j in 1:length(theta))
    h[,j] <- sapply(c, function(c.)
        qrmtools:::Wang_h(c., alpha=alpha, d=d, method="Wang.Par", theta=theta[j]))
z <- h
z[z <= 0] <- NA # > 0 => makes log-scale possible

## ------------------------------------------------------------------------
if(doPDF)
    pdf(file=(file <- paste0("fig_worst_VaR_",alpha,"_hom_Wang_h_Par_d=",d,".pdf")),
        width=6, height=6)
plot(NA, xlim=range(c), ylim=range(z, na.rm=TRUE), log="y", yaxt="n", xlab="c",
     ylab=substitute("h(c) for"~~alpha~"= 0.99, d ="~d.~"Par("*theta*") margins",
                     list(d.=d)))
eaxis(2)
for(k in 1:length(theta))
    lines(c, z[,k], col=cols[k])
legend("topleft", bty="n", lty=rep(1,length(theta)), col=cols,
       legend=as.expression(lapply(1:length(theta),
       function(k) substitute(theta==k, list(k=theta[k])))))
if(doPDF) dev.off()

## ------------------------------------------------------------------------
d <- dim # dimension
alpha <- 1-2^seq(-0.001, -10, length.out=128) # confidence levels; concentrated near 1
theta <- c(0.1, 0.5, 1, 5, 10, 50) # theta values
VaR <- simplify2array(sapply(alpha, function(a)
    sapply(theta, function(th) VaR_bounds_hom(a, d=d, method="Wang.Par",
                                              theta=th)), simplify=FALSE))
## => (best/worst VaR, theta, alpha)-matrix

## ------------------------------------------------------------------------
if(doPDF)
    pdf(file=(file <- paste0("fig_worst_VaR_hom_Wang_Par_d=",d,".pdf")),
        width=6, height=6)
plot(NA, xlim=range(alpha), ylim=range(VaR), log="xy", yaxt="n",
     xlab=expression(1-alpha),
     ylab=as.expression(substitute(VaR[alpha]~"for d ="~d.~"and Par("*theta*") margins",
                                   list(d.=d))))
eaxis(2)
for(k in 1:length(theta)) {
    lines(1-alpha, VaR[2,k,], col=cols[k]) # worst VaR
    lines(1-alpha, VaR[1,k,], col=cols[k], lty=2) # best VaR
}
legend("topright", bty="n", lty=rep(1,length(theta)), col=cols,
       legend=as.expression(lapply(1:length(theta),
       function(k) substitute(theta==k, list(k=theta[k])))))
mtext("Solid line: worst; dashed line: best VaR",
      side=4, line=1, adj=0)
if(doPDF) dev.off()

## ------------------------------------------------------------------------
d <- seq(2, 1002, by=20) # dimensions
alpha <- 0.99 # confidence level
theta <- c(0.1, 0.5, 1, 5, 10, 50) # theta values
VaR <- simplify2array(sapply(d, function(d.)
    sapply(theta, function(th) VaR_bounds_hom(alpha, d=d., method="Wang.Par",
                                              theta=th)), simplify=FALSE))
## => (best/worst VaR, theta, alpha)-matrix

## ------------------------------------------------------------------------
if(doPDF)
    pdf(file=(file <- paste0("fig_worst_VaR_hom_Wang_Par_alpha=",alpha,".pdf")),
        width=6, height=6)
plot(NA, xlim=range(d), ylim=range(VaR), log="xy", yaxt="n",
     xlab=expression(d),
     ylab=as.expression(substitute(VaR[a]~"for Par("*theta*") margins",
                                   list(a=alpha))))
eaxis(2)
for(k in 1:length(theta)) {
    lines(d, VaR[2,k,], col=cols[k]) # worst VaR
    lines(d, VaR[1,k,], col=cols[k], lty=2) # best VaR
}
legend("topleft", bty="n", lty=rep(1,length(theta)), col=cols,
       legend=as.expression(lapply(1:length(theta),
       function(k) substitute(theta==k, list(k=theta[k])))))
mtext("Solid line: worst; dashed line: best VaR",
      side=4, line=1, adj=0)
if(doPDF) dev.off()

## ------------------------------------------------------------------------
## Initial interval for the root finding in case of worst VaR
init_interval <- function(alpha, d, theta, trafo=FALSE, adjusted=FALSE)
{
    if(trafo) {
        low <- if(theta == 1) {
            d/2
        } else {
            (d-1)*(1+theta)/(d-1+theta)
        }
        up <- if(theta > 1) {
            r <- (1+d/(theta-1))^theta
            if(adjusted) 2*r else r
        } else if(theta == 1) {
            e <- exp(1)
            (d+1)^(e/(e-1))
        } else {
            d*theta/(1-theta)+1
        }
        c(low, up)
    } else {
        low <- if(theta > 1) {
            r <- (1-alpha)/((d/(theta-1)+1)^theta + d-1)
            if(adjusted) r/2 else r
        } else if(theta == 1) {
            e <- exp(1)
            (1-alpha)/((d+1)^(e/(e-1))+d-1)
        } else {
            r <- (1-theta)*(1-alpha)/d
            if(adjusted) r/2 else r
        }
        up <- if(theta == 1) (1-alpha)/(3*d/2-1)
              else (1-alpha)*(d-1+theta)/((d-1)*(2*theta+d))
        c(low, up)
    }
}

## Function to compute the best/worst Value-at-Risk in the homogeneous case with
## Par(theta) margins
VaR_hom_Par <- function(alpha, d, theta, method=c("worst", "best"),
                        trafo=FALSE, interval=NULL, adjusted=FALSE,
                        avoid.cancellation=FALSE, ...)
{
    ## Pareto quantile function
    qF <- function(p) (1 - p)^(-1/theta) - 1

    ## Compute \bar{I}
    Ibar <- function(a, b, alpha, d, theta)
    {
        if(theta == 1) log((1-a)/(1-b))/(b-a) - 1
        else (theta/(1-theta))*((1-b)^(1-1/theta)-(1-a)^(1-1/theta))/(b-a) - 1
    }

    ## Main
    method <- match.arg(method)
    switch(method,
    "worst" = {

        ## Distinguish according to whether we optimize the auxiliary function
        ## on a transformed scale
        h <- if(trafo) {
            ## Auxiliary function to find the root of on (1, Inf)
            if(theta == 1) {
                function(x) x^2 + x*(-d*log(x)+d-2)-(d-1)
            } else {
                function(x)
                (d/(1-theta)-1)*x^(-1/theta + 1) - (d-1)*x^(-1/theta) + x - (d*theta/(1-theta) + 1)
            }
        } else {
            ## Auxiliary function to find the root of on (0, (1-alpha)/d)
            function(c) {
                a <- alpha+(d-1)*c
                b <- 1-c
                Ib <- if(c == (1-alpha)/d) { # Properly deal with limit c=(1-alpha)/d
                    ((1-alpha)/d)^(-1/theta) - 1
                } else {
                    Ibar(a=a, b=b, alpha=alpha, d=d, theta=theta)
                }
                Ib - (qF(a)*(d-1)/d + qF(b)/d)
            }
        }

        ## Do the optimization
        if(is.null(interval)) interval <- init_interval(alpha, d, theta,
                                                        trafo=trafo, adjusted=adjusted)
        c <- uniroot(h, interval=interval, ...)$root
        if(trafo) # convert value back to the right scale (c-scale)
            c <- (1-alpha)/(c+d-1)
        if(avoid.cancellation) {
            t1 <- (1-alpha)/c-(d-1)
            d * ((c^(-1/theta)/d) * ((d-1)*t1^(-1/theta) + 1) - 1) # = qF(a)*(d-1) + qF(b)
        } else {
            a <- alpha+(d-1)*c
            b <- 1-c
            qF(a)*(d-1) + qF(b)
        }
    },
    "best" = {
        max((d-1)*0 + (1-alpha)^(-1/theta)-1, # Note: Typo in Wang, Peng, Yang (2013)
            d*Ibar(a=0, b=alpha, alpha=alpha, d=d, theta))
    },
    stop("Wrong 'method'"))
}

## ------------------------------------------------------------------------
alpha <- 0.99 # confidence level
d <- dim # dimension
n.th <- 32 # number of thetas
th <- seq(0.2, 5, length.out=n.th) # thetas
qFs <- lapply(th, function(th.) {th.; function(p) qPar(p, theta=th.)}) # n.th-vector of Pareto quantile functions
pFs <- lapply(th, function(th.) {th.; function(q) pPar(q, theta=th.)}) # n.th-vector of Pareto dfs
N <- 1e4 # number of discretization points for RA(); N=1e5 does not improve the situation

## ---- results="hide", warning=FALSE--------------------------------------
res <- matrix(, nrow=n.th, ncol=7)
colnames(res) <- c("Wang", "straightforward", "transformed", "Wang.Par",
                   "dual", "RA.low", "RA.up")
pb <- txtProgressBar(max=n.th, style=if(isatty(stdout())) 3 else 1) # setup progress bar
on.exit(close(pb)) # on exit, close progress bar
for(i in seq_len(n.th)) {
    ## "Wang" (numerical integration with smaller uniroot() tolerance; still
    ## numerically critical -- we catch "the integral is probably divergent"-errors here)
    Wang.num.res <- tryCatch(VaR_bounds_hom(alpha, d=d, qF=qFs[[i]])[2], error=function(e) e)
    res[i,"Wang"] <- if(is(Wang.num.res, "simpleError")) NA else Wang.num.res
    ## Our straightforward implementation
    res[i,"straightforward"] <- VaR_hom_Par(alpha, d=d, theta=th[i])
    ## Our straightforward implementation based on the transformed auxiliary function
    res[i,"transformed"] <- VaR_hom_Par(alpha, d=d, theta=th[i], trafo=TRUE)
    ## "Wang.Par" (internally using a transformed objective function and
    ## avoids cancellation)
    res[i,"Wang.Par"] <- VaR_bounds_hom(alpha, d=d, method="Wang.Par", theta=th[i])[2]
    ## "dual" (with uniroot()'s default tolerance)
    res[i,"dual"] <- VaR_bounds_hom(alpha, d=d, method="dual",
                                    interval=crude_VaR_bounds(alpha, qF=rep(qFs[i], d)),
                                    pF=pFs[[i]])[2]
    ## Rearrangement Algorithm
    set.seed(271) # use the same random permutation for each theta
    RA. <- RA(alpha, qF=rep(qFs[i], d), N=N)
    res[i,"RA.low"] <- RA.$bounds[1]
    res[i,"RA.up"]  <- RA.$bounds[2]
    ## Progress
    setTxtProgressBar(pb, i) # update progress bar
}

## ---- fig.width=10, fig.height=7-----------------------------------------
res. <- res/res[,"dual"] # standardize
ylim <- range(res., na.rm=TRUE)
if(doPDF)
    pdf(file=(file <- paste0("fig_worst_VaR_",alpha,"_hom_comparison_d=",
                             d,"_N=",N,".pdf")), width=6, height=6)
plot(th, res.[,"Wang"], type="l", ylim=ylim,
     xlab=expression(theta), ylab=substitute("Standardized (by dual method) worst"~
     VaR[0.99]~"for d ="~d.~"Par("*theta*") margins", list(d.=d)),
     col="gray60", lty=2, lwd=5.5) # Wang (with numerical integration)
lines(th, res.[,"straightforward"], col="maroon3", lty=1, lwd=1) # still bad (although we have an initial interval)
lines(th, res.[,"transformed"], col="black", lty=1, lwd=1) # okay
lines(th, res.[,"Wang.Par"], col="royalblue3", lty=2, lwd=2.5) # Wang Pareto (wo num. integration)
lines(th, res.[,"RA.low"], col="black", lty=3, lwd=1)
lines(th, res.[,"RA.up"],  col="black", lty=2, lwd=1)
legend("topright", bty="n",
       col=c("gray60", "maroon3", "black", "royalblue3", "black", "black"),
       lty=c(2,1,1,2,3,2), lwd=c(5.5,1,1,2.5,1,1),
       legend=c("Wang (num. int.)", "Wang Pareto (straightforward)",
                "Wang Pareto (transformed)", "Wang Pareto (wo num. int.)",
                "Lower RA bound", "Upper RA bound"))
if(doPDF) dev.off()

## ------------------------------------------------------------------------
tol <- 2.2204e-16
wVaR.tol <- sapply(th, function(th.)
    VaR_hom_Par(alpha=alpha, d=d, theta=th., tol=tol))
plot(th, wVaR.tol/res[,"dual"], type="l", ylim=ylim,
     xlab=expression(theta), ylab="Wang's approach (straightforward) but with smaller tol")

## ---- warning=FALSE------------------------------------------------------
alpha <- 0.99 # confidence level
d <- seq(2, 1002, by=20) # dimensions
theta <- c(0.1, 0.5, 1, 5, 10, 50) # theta values
VaR <- simplify2array(sapply(d, function(d.)
    sapply(theta, function(th) {
        res <- tryCatch(VaR_hom_Par(alpha, d=d., theta=th, tol=tol),
                        error=function(e) e)
        if(is(res, "simpleError")) { warning(conditionMessage(res)); NA }
        else res
    }), simplify=FALSE))

## ------------------------------------------------------------------------
VaR <- simplify2array(sapply(d, function(d.)
    sapply(theta, function(th) VaR_hom_Par(alpha, d=d., theta=th, tol=tol, adjusted=TRUE)),
    simplify=FALSE))

## ------------------------------------------------------------------------
d <- 500
th <- 20
VaR_hom_Par(alpha, d=d, theta=th, trafo=TRUE) # Inf

## ------------------------------------------------------------------------
h <- function(x)
     (d/(1-th)-1)*x^(-1/th + 1) - (d-1)*x^(-1/th) + x - (d*th/(1-th) + 1)
interval <- init_interval(alpha, d, th, trafo=TRUE)
x <- uniroot(h, interval=interval)$root
(c <- (1-alpha)/(x+d-1)) # convert back to c-scale
a <- alpha+(d-1)*c
b <- 1-c
qPar(a, theta=th)*(d-1) + qPar(b, theta=th) # Inf
stopifnot(b==1) # => b is 1 => qPar(b, theta=th) = Inf

## ------------------------------------------------------------------------
qPar(a, theta=th)*(d-1) + c^(-1/th)-1

## ------------------------------------------------------------------------
d <- seq(2, 1002, by=20) # dimensions
theta <- c(0.1, 0.5, 1, 5, 10, 50) # theta values
VaR. <- simplify2array(sapply(d, function(d.)
    sapply(theta, function(th) VaR_hom_Par(alpha, d=d., theta=th,
                                           trafo=TRUE, adjusted=TRUE,
                                           avoid.cancellation=TRUE)),
    simplify=FALSE))

## ------------------------------------------------------------------------
ylim <- range(VaR, VaR.)
if(doPDF)
    pdf(file=(file <- paste0("fig_worst_VaR_",alpha,"_hom_comparison_num_problems.pdf")),
        width=6, height=6)
plot(d, VaR[1,], type="l", ylim=ylim, log="y", yaxt="n",
     xlab="d", ylab=substitute("Worst"~VaR[a]~"for Par("*theta*") margins",
                               list(a=alpha)), col=cols[1])
eaxis(2)
for(k in 2:length(theta)) lines(d, VaR [k,], col=cols[k])
for(k in 1:length(theta)) lines(d, VaR.[k,], col=cols[k], lty=2, lwd=2)
legend("right", bty="n",
       col=cols, lty=rep(1, length(theta)),
       legend=as.expression( lapply(1:length(theta), function(k)
           substitute(theta==th, list(th=theta[k]))) ))
mtext("Solid line: direct; dashed line: transformed h",
      side=4, line=1, adj=0)
if(doPDF) dev.off()

## ------------------------------------------------------------------------
rearrange_ <- function(X, tol=0)
{
    N <- nrow(X)
    d <- ncol(X)
    m.rs.old <- min(rowSums(X))
    while (TRUE) {
        Y <- X
        for(j in 1:d)
            Y[,j] <- sort(Y[,j])[N+1-rank(rowSums(Y[,-j, drop=FALSE]), ties.method="first")]
        Y.rs <- rowSums(Y)
        m.rs.new <- min(Y.rs)
        tol. <- abs((m.rs.new - m.rs.old)/m.rs.old)
        tol.reached <- if(is.null(tol)) {
            identical(Y, X)
        } else { tol. <= tol }
        if(tol.reached) {
            break
        } else {
            m.rs.old <- m.rs.new
            X <- Y
        }
    }
}

## ------------------------------------------------------------------------
## Build the input matrix (note that the columns have to be sorted in increasing order!)
alpha <- 0.99 # confidence level
N <- 2^9 # 512
p <- alpha + (1-alpha)*0:(N-1)/N
if(FALSE)
    d <- 2^(4:10) # 16, ..., 1024
d <- 2^(4:8) # 16, ..., 256 (to save time here)
res <- matrix(, nrow=length(d), ncol=2)
colnames(res) <- c("slow", "fast")
qF <- function(p, theta=2) (1-p)^(-1/theta)-1

## For each d, measure the run time
for(i in seq_along(d)) {
    ## cat("Working on d =",d[i],"\n")
    X <- sapply(rep(list(qF), d[i]), function(qF.) qF.(p))
    res[i,] <- c(system.time(rearrange_(X))[["elapsed"]],
                 system.time(rearrange(X, sample=FALSE, is.sorted=TRUE))[["elapsed"]])
}

## ------------------------------------------------------------------------
if(doPDF)
    pdf(file=(file <- paste0("fig_ARA_speed-up.pdf")),
        width=6, height=6)
plot(d, (1-res[,2]/res[,1])*100, type="b", log="x", ylim=c(0,100),
     xlab="d", ylab="Relative speed-up (in %) of implemented ARA() over a basic version")
if(doPDF) dev.off()

## ------------------------------------------------------------------------
A <- matrix(c(1:4, 2*(1:4)-1, 2^(0:3)), ncol=3)
rearrange(A, tol=NULL, sample=FALSE, is.sorted=TRUE, trace=TRUE)

## ------------------------------------------------------------------------
B <- matrix(rep(1:3, 3), ncol=3)
rearrange(B, tol=NULL, sample=FALSE, is.sorted=TRUE, trace=TRUE)

## ------------------------------------------------------------------------
### Investigate the probability of certain rearranged matrices to appear #######

##' @title Create a list of all possible input matrices for d=3
##' @param N The number of discretization points
##' @param method The method how to create the list
##' @return A N!^2-list (containing all possible (N,d) input matrices with 1:N in
##'         the first column)
##' @author Marius Hofert
create_mat <- function(N)
{
    x_perm <- permn(N) # N!-long list containing all possible permutations of 1:N
    N. <- factorial(N)
    lst <- vector("list", length=N.^2)
    cnt <- 0
    for (i in 1:N.) { # double 'for' not time critical here
        for (j in 1:N.){
            cnt <- cnt+1
            lst[[cnt]] <- cbind(1:N, x_perm[[i]], x_perm[[j]])
        }
    }
    lst
}


## Main
N <- 5 # chosen N (<= 6 due to extensive run time)
system.time(mat  <- create_mat(N)) # create the list of matrices (N!^2-many)

## Rearrange all of them
N. <- factorial(N)^2
system.time(mat.ra <- lapply(1:N., function(i) {
                ## if(i %% (N./20) == 0) cat(round(100*i/N.),"% done!\n", sep="");
                rearrange(mat[[i]], tol=NULL, sample=FALSE)$X.rearranged
            }))

## Go through all of the rearranged matrices, 'uniquify' them and count
lst <- list(mat.ra[[1]])
freq <- c(1)
for(i in 2:length(mat.ra))
{
    bool <- sapply(1:length(lst),
                   function(i) { identical(lst[[i]], mat.ra[[i]]) }) # is mat.ra in lst?
    if(any(bool)) { # mat.ra has been observed before
        freq[min(which(bool))] <- freq[min(which(bool))] + 1
    } else { # mat.ra has not been observed before
        lst <- c(lst, mat.ra[[i]]) # append it to lst()
        freq <- c(freq, 1) # append its frequency
    }
}

## Result
lst # final rearranged matrices
freq/N. # their corresponding probability

## ------------------------------------------------------------------------
data(SMI.12)
L <- -apply(log(SMI.12), 2, diff) # more sophisticated methods are available
n <- nrow(L)
d <- ncol(L)

## ------------------------------------------------------------------------
res <- vector("list", length=d)
names(res) <- colnames(L)
for(k in seq_len(d)) {

    ## Determine the threshold for company k
    L. <- L[,k]
    u <- quantile(L., probs=0.8, names=FALSE) # threshold

    ## Fit a GPD to the excesses
    fit <- fit.GPD(L., threshold=u)
    stopifnot(fit$converged)
    xi <- fit$par.ests[["xi"]] # fitted xi
    beta <- fit$par.ests[["beta"]] # fitted beta
    stopifnot(is.numeric(xi), is.numeric(beta))

    ## Graphical goodness-of-fit check for the GPD fit
    excess <- L.[L.>u]-u # excesses
    if(FALSE) {
        excess. <- sort(excess) # sorted data
        qF <- function(p) qGPD(p, xi=xi, beta=beta)
        qF. <- qF(ppoints(length(excess.))) # theoretical quantiles
        stock <- names(res)[k]
        plot(qF., excess., xlab="Theoretical quantiles",
             ylab="Sample quantiles", main=paste0("Q-Q plot for the fitted GPD(",
             round(xi, 2),", ",round(beta, 2),") distribution for ", stock))
        qqline(y=as.numeric(excess.), distribution=qF)
    }

    ## Update res
    res[[k]] <- list(loss=L., excess=excess, u=u, xi=xi, beta=beta)

}

## ------------------------------------------------------------------------
xi. <- sapply(res, function(x) x$xi) # all fitted xi's
beta. <- sapply(res, function(x) x$beta) # all fitted beta's
alpha <- 0.999 # confidence level
qF <- lapply(res, function(r) { function(p) qGPD(p, xi=r$xi, beta=r$beta) }) # list quantile functions
set.seed(271) # set a seed (for reproducibility)
res.ARA <- ARA(alpha, qF=qF) # apply ARA()
stopifnot(res.ARA$converged) # check convergence
res.ARA$bounds # ARA bounds on worst VaR
res.RA <- RA(alpha, qF=qF, N=res.ARA$N.used, abstol=0) # apply RA() (with the same N)
stopifnot(res.RA$converged) # check convergence
res.RA$bounds # RA bounds on worst VaR

## ------------------------------------------------------------------------
## Build the input matrices for rearrange()
N <- res.ARA$N.used
p.low <- alpha + (1-alpha)*(0:(N-1))/N # probability points for worst VaR (lower bound)
p.up  <- alpha + (1-alpha)*(1:N)/N # probability points for worst VaR (upper bound)
X.low <- sapply(qF, function(qF) qF(p.low)) # input matrix (lower bound)
X.up  <- sapply(qF, function(qF) qF(p.up)) # input matrix (upper bound)
X.up[N,] <- sapply(1:d, function(j) if(is.infinite(X.up[N,j]))
                  qF[[j]](alpha+(1-alpha)*(1-1/(2*N))) else X.up[N,j])

## Apply rearrange() for various tolerances (incl. NULL)
tol <- seq(0, 0.02, length.out=21) # considered tolerances
res.low  <- sapply(tol, function(t) rearrange(X.low,  tol=t, sample=FALSE, is.sorted=TRUE)$bound)
res.low. <- rearrange(X.low, tol=NULL, sample=FALSE, is.sorted=TRUE)$bound # for tol=NULL
res.up  <- sapply(tol, function(t) rearrange(X.up, tol=t, sample=FALSE, is.sorted=TRUE)$bound)
res.up. <- rearrange(X.up, tol=NULL, sample=FALSE, is.sorted=TRUE)$bound # for tol=NULL

## ------------------------------------------------------------------------
## Plot the lower and upper bound on worst VaR as a function in the chosen tol
if(doPDF)
    pdf(file=(file <- paste0("fig_worst_VaR_bounds_application.pdf")),
        width=6, height=6)
plot(tol, res.low, type="b", log="y", ylim=range(res.low, res.low., res.up, res.up.),
     col="royalblue3", xlab="Relative tolerance tol of rearrange()",
     ylab=substitute("Bounds on worst"~VaR[a], list(a=alpha)))
points(0, res.low., pch=3, col="royalblue3") # draw tol=NULL result at 0, too (as '+')
lines(tol, res.up, type="b")
points(0, res.up., pch=3) # draw tol=NULL result at 0, too (as '+')
legend("topright", bty="n", lty=rep(1,2),
       col=c("black", "royalblue3"), legend=c("upper bound", "lower bound"))
if(doPDF) dev.off()
