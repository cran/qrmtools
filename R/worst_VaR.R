### Tools for computing the worst VaR_alpha for given margins ##################

### 1) Crude VaR bounds (for both best and worst VaR) ##########################

##' @title Crude bounds for any VaR_alpha
##' @param alpha confidence level
##' @param d dimension
##' @param qF (list of) marginal quantile functions
##' @param ... ellipsis argument passed to qF()
##' @return 2-vector containing crude VaR_alpha bounds
##' @author Marius Hofert
crude_VaR_bounds <- function(alpha, d, qF, ...)
{
    if(is.function(qF))
        d * c(qF(alpha/d, ...), qF((d-1+alpha)/d, ...))
    else { # ... are passed to *all* qF()
        if(!is.list(qF))
            stop("qF has to be a (quantile) function or list of such")
        stopifnot(length(qF) == d)
        qF.low <- sapply(qF, function(qF.) qF.(alpha/d, ...))
        qF.up  <- sapply(qF, function(qF.) qF.((d-1+alpha)/d, ...))
        d * c(min(qF.low), max(qF.up))
    }
}


### 2) Explicit worst VaR in the homogeneous case ##############################

### Dual bound #################################################################

##' @title D(s,t) = d \int_{t}^{s-(d-1)t} \bar{F}(x) dx / (s-dt)
##' @param s real number
##' @param t real number < s/d
##' @param d dimension (integer > 2)
##' @param pF marginal distribution function (same for all d)
##' @param ... ellipsis argument passed to integrate()
##' @return D(s,t)
##' @author Marius Hofert
##' @note If t -> s/d-, l'Hospital's Rule shows that D(s, s/d) = d\bar{F}(s/d)
dual_bound_2 <- function(s, t, d, pF, ...)
{
    stopifnot(length(t) == 1)
    if(t > s/d) stop("t must be <= s/d")
    ## use D(s,t) = d( 1-\int_{t}^{s-(d-1)t} F(x) dx/(s-d*t) ) in this case
    if(t == s/d) d*(1-pF(s/d)) else
    d * (1 - (1/(s-d*t)) * integrate(pF, lower=t, upper=s-(d-1)*t, ...)$value)
}

##' @title Auxiliary function \bar{F}(t) + (d-1) * \bar{F}(s-(d-1)*t)
##' @param s real number
##' @param t real number < s/d
##' @param d dimension (integer > 2)
##' @param pF marginal distribution function (same for all d)
##' @return \bar{F}(t) + (d-1) * \bar{F}(s-(d-1)*t)
##' @author Marius Hofert
dual_bound_2_deriv_term <- function(s, t, d, pF)
    1-pF(t) + (d-1)*(1-pF(s-(d-1)*t))

##' @title Dual bound D(s)
##' @param s real number
##' @param d dimension (integer > 2)
##' @param pF marginal distribution function (same for all d)
##' @param ... ellipsis argument passed to dual_bound_2()'s integrate()
##' @return D(s)
##' @author Marius Hofert
##' @note The "first-order condition" (second equality in (14) in 2)) comes from the
##'       fact that
##'       (d/dt) D(s,t) = [ (-d)[\bar{F}(s-(d-1)t)(d-1)+\bar{F}(t)](s-dt) +
##'                         d^2 \int_{t}^{s-(d-1)t} \bar{F}(x) dx ] / (s-dt)^2 = 0
##'       if and only if
##'       d (\int_{t}^{s-(d-1)t} \bar{F}(x) dx) / (s-dt) = \bar{F}(s-(d-1)t)(d-1)-\bar{F}(t)
##'       => solving d (\int_{t}^{s-(d-1)t} \bar{F}(x) dx) / (s-dt) -
##'                  (\bar{F}(s-(d-1)t)(d-1)-\bar{F}(t)) = 0
##'          as a function in t for sufficiently large s leads to D(s)
dual_bound <- function(s, d, pF, ...)
{
    stopifnot(length(s) == 1, s >= 0)
    if(s > 0) {
        ## h(s, t)
        h <- function(t) dual_bound_2(s, t=t, d=d, pF=pF, ...) -
            dual_bound_2_deriv_term(s, t=t, d=d, pF=pF)
        ## note: f(t) -> 0 for t -> s/d- this is bad for uniroot() which will
        ##       simply stop with the root t=s/d => we thus set f.upper > 0
        ##       alternatively, one could use (the more invasive)
        ##       interval = c(0, s/d - .Machine$double.eps^0.25)
        h.up <- sign(-h(0)) * .Machine$double.xmin # guarantee that uniroot() doesn't fail due to root s/d
        t. <- uniroot(h, interval=c(0, s/d), f.upper=h.up)$root # optimal t in Equation (12) [= arginf]
        dual_bound_2_deriv_term(s, t=t., d=d, pF=pF) # dual bound D(s) in Equation (12) [= inf]
    } else {
        ## if s = 0, then t in [0, s/d] requires t to be 0 *and* f(0) = 0, so
        ## 0 is a root (as s/d). Furthermore, at t=0 (and with s=0),
        ## dual_bound_2_deriv_term(...) = d
        d
    }
}


### Wang's methods #############################################################

##' @title Conditional expectation (\bar{I}(c)) for computing the worst VaR as
##'        in Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014, Prop. 3.1)
##' @param c evaluation point
##' @param alpha confidence level alpha
##' @param d dimension d
##' @param method character string giving the method
##'        generic = numerical integration; Wang.Par = Pareto distibution
##' @param ... ellipsis argument passed to integrate()
##' @return \bar{I}(c)
##' @author Marius Hofert
Wang_Ibar <- function(c, alpha, d, method=c("generic", "Wang.Par"), ...)
{
    a <- alpha + (d-1)*c
    b <- 1-c
    ddd <- list(...)
    method <- match.arg(method)
    switch(method,
           "generic" = {
               stopifnot(length(c) == 1) # not vectorized (due to integrate())
               qF <- ddd$qF # use provided 'qF()'
               ## properly deal with the limit c == (1-alpha)/d (empty integration interval)
               if(a == b) { # or c == (1-alpha)/d
                   qF((d-1+alpha)/d)
               } else { # now the non-empty case
                   ddd$qF <- NULL
                   h <- function(...)
                       integrate(qF, lower=a, upper=b, ...)$value / (b-a)
                   do.call(h, ddd) # call integrate() on the remaining arguments in '...'
               }
           },
           "Wang.Par" = {
               th <- ddd$theta # use provided 'theta'
               ## properly deal with the limit c == (1-alpha)/d (empty integration interval)
               ## (vectorized)
               res <- rep(NA, length(c))
               i <- a == b # or c == (1-alpha)/d
               if(any(i))   res[i] <- rep(qPar((d-1+alpha)/d, theta=th), sum(i))
               if(any(!i)) res[!i] <- if(th == 1) log((1-a[!i])/(1-b[!i]))/(b[!i]-a[!i]) - 1
                   else (th/(1-th))*((1-b[!i])^(1-1/th)-(1-a[!i])^(1-1/th))/(b[!i]-a[!i]) - 1
               res
           },
           stop("Wrong method"))
}

##' @title Right-hand side term in the objective function for computing the worst VaR
##'        as in Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014, Prop. 3.1)
##' @param c evaluation point
##' @param alpha confidence level alpha
##' @param d dimension d
##' @param method character string giving the method
##'        generic = numerical integration; Wang.Par = Pareto distibution
##' @param ... ellipsis argument containing theta (for method="Wang.Par")
##'        or qF (for method="generic")
##' @return Right-hand side term in Prop. 3.1
##' @author Marius Hofert
##' @note for the correct 'c', this is the conditional expectation
Wang_h_aux <- function(c, alpha, d, method=c("generic", "Wang.Par"), ...)
{
    ddd <- list(...)
    method <- match.arg(method)
    qF <- if(method=="Wang.Par") function(y) qPar(y, theta=ddd$theta) else ddd$qF
    a <- alpha + (d-1)*c
    b <- 1-c
    qF(a)*(d-1)/d + qF(b)/d
}

##' @title Objective function for computing the worst VaR as in
##'        Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014, Prop. 3.1)
##' @param c evaluation point
##' @param alpha confidence level alpha
##' @param d dimension d
##' @param method character string giving the method
##' @param ... ellipsis argument passed to Wang_h_aux() and Wang_Ibar()
##' @return objective function for computing the worst VaR
##' @author Marius Hofert
Wang_h <- function(c, alpha, d, method=c("generic", "Wang.Par"), ...)
{
    stopifnot(0 <= c, c <= (1-alpha)/d) # sanity check (otherwise b > a)
    Wang_Ibar(c, alpha=alpha, d=d, method=method, ...) -
    Wang_h_aux(c, alpha=alpha, d=d, method=method, ...)
}


### Main wrapper function for computing the worst VaR in the homogeneous case ##

## Assumptions:
## - d=2: ultimately decreasing density (for x >= x0), alpha >= F(x0)
## - "Wang": F needs to live on [0, Inf), admitting a positive density which is
##           ultimately decreasing (for x >= x0), alpha >= F(x0)
## - "dual": F needs to be continuous with unbounded support and and ultimately
##           decreasing density, F(0) = 0 (otherwise, 0 as a lower bound for
##           uniroot() in dual_bound() is not valid)

##' @title Compute the worst VaR_\alpha in the homogeneous case with:
##'        1) d=2: Embrechts, Puccetti, Rueschendorf (2013, Proposition 2)
##'        2) d>=3:
##'           "Wang": Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014, Prop. 3.1)
##'           "Wang.Par": the same, just with explicit formula for the integral
##'                       in the Pareto case
##'           "dual": Embrechts, Puccetti, Rueschendorf (2013, Proposition 4)
##' @param alpha confidence level
##' @param d dimension
##' @param method character string giving the method
##' @param interval initial interval
##' @param tol uniroot() x-tolerance
##' @param ... ellipsis arguments passed to Wang_h() (for d>=3), dual_bound()
##' @return worst VaR in the homogeneous case
##' @author Marius Hofert
worst_VaR_hom <- function(alpha, d, method=c("Wang", "Wang.Par", "dual"),
                          interval=NULL, tol=NULL, ...)
{
    stopifnot(0<alpha, alpha<1, d>=2)
    method <- match.arg(method)
    if(is.null(tol)) # use higher tolerance (MATLAB default) for Wang's approach (matters!)
        tol <- if(method!="dual") 2.2204e-16 else .Machine$double.eps^0.25 # uniroot()'s default

    ## d == 2
    if(d == 2) { # see Embrechts, Puccetti, Rueschendorf (2013, Proposition 2)
        qF <- NULL # make CRAN check happy
        if(!hasArg(qF))
            stop("The bivariate case requires the quantile function qF of F")
        return(2*list(...)$qF((1+alpha)/2))
    }

    ## d >= 3
    switch(method,
           "Wang" = {
               ## check qF()
               qF <- NULL # make CRAN check happy
               if(!hasArg(qF))
                   stop("Method 'Wang' requires the quantile function qF of F")
               ## check 'interval'
               if(is.null(interval)) interval <- c(0, (1-alpha)/d)
               else {
                   if(interval[1] < 0) stop("interval[1] needs to be >= 0")
                   if(interval[1] > (1-alpha)/d) stop("interval[2] needs to be <= (1-alpha)/d")
                   if(interval[1] >= interval[2]) stop("interval[1] needs to be smaller than interval[2]")
               }
               ## c_u: guarantee that uniroot() doesn't fail due to root at (1-alpha)/d
               h.up <- .Machine$double.xmin
               ## c_l: compute and check validity
               h.low <- Wang_h(interval[1], alpha=alpha, d=d, ...)
               if(is.na(h.low))
                   stop("Objective function at interval[1] is NA/NaN. Provide a larger interval[1].")
               if(h.up * h.low >= 0)
                   stop("Objective function at end points of 'interval' not of opposite sign. Provide a smaller interval[1].")
               ## root-finding on 'interval'
               c. <- uniroot(function(c) Wang_h(c, alpha=alpha, d=d, ...),
                             interval=interval, f.lower=h.low, f.upper=h.up, tol=tol)$root
               d * Wang_h_aux(c., alpha=alpha, d=d, qF=list(...)$qF)
           },
           "Wang.Par" = { # only kept for ruling out problems due to numerical integration
               ## check 'theta'
               theta <- NULL # make CRAN check happy
               if(!hasArg(theta))
                   stop("Method 'Wang.Par' requires the parameter theta")
               th <- list(...)$theta
               stopifnot(length(th) == 1, th > 0) # check theta here
               ## check 'interval'
               if(is.null(interval)) # note: we can be more specific here
                   interval <- c(if(th <= 1) .Machine$double.eps else 0, (1-alpha)/d)
               else {
                   if(interval[1] < 0) stop("interval[1] needs to be >= 0")
                   if(interval[1] > (1-alpha)/d) stop("interval[2] needs to be <= (1-alpha)/d")
                   if(interval[1] >= interval[2]) stop("interval[1] needs to be smaller than interval[2]")
               }
               ## c_u: guarantee that uniroot() doesn't fail due to root at (1-alpha)/d
               h.up <- .Machine$double.xmin
               ## c_l: compute and check validity (note: we can be more specific here)
               if(th <= 1 && interval[1] == 0)
                   stop("If theta <=1, interval[1] has to be > 0 as otherwise the internal Wang_h() is NaN")
               h.low <- Wang_h(interval[1], alpha=alpha, d=d, method="Wang.Par", ...)
               ## former check (not necessary anymore due to check theta <= 1 && interval[1] == 0)
               ## if(is.na(h.low))
               ##     stop("Objective function at interval[1] is NA/NaN. Provide a larger interval[1].")
               if(h.up * h.low >= 0)
                   stop("Objective function at end points of 'interval' not of opposite sign. Provide a smaller interval[1].")
               ## root-finding on 'interval'
               c. <- uniroot(function(c)
                             Wang_h(c, alpha=alpha, d=d, method="Wang.Par", ...),
                             interval=interval, f.lower=h.low, f.upper=h.up, tol=tol)$root
               d * Wang_h_aux(c., alpha=alpha, d=d, method="Wang.Par", theta=th)
           },
           "dual" = {
               pF <- NULL # make CRAN check happy
               if(!hasArg(pF))
                   stop("Method 'dual' requires the distribution function pF")
               if(!hasArg(interval))
                   stop("Method 'dual' requires an initial interval")
               uniroot(function(s) dual_bound(s, d=d, ...)-(1-alpha),
                       interval=interval, tol=tol)$root # s interval
               ## note: we can't pass arguments to the inner root-finding (yet)
           },
           stop("Wrong method"))
}


### 3) Worst VaR in the inhomogeneous case #####################################

### Rearrangement algorithm ####################################################

##' @title Auxiliary Function for Computing Steps 4 and 5 for the RA
##' @param X (N, d)-matrix (either \underline{X}^\alpha or \overline{X}^\alpha)
##' @param method character indicating which VaR is approximated (worst/best)
##'        determines optimizing function (min for worst VaR; max
##'        for best VaR)
##' @param err character string indicating the error function used
##'        ("absolute" or "relative")
##' @param maxiter maximal number of iterations
##' @param eps epsilon error to determine convergence; if NULL, then the
##'        iteration is done until the matrix doesn't change
##' @param impl string indicating the implementation (either "C" or "R")
##' @return 5-list containing the
##'         1) computed (lower or upper [depending on X]) bound for (worst or
##'            best [depending on method]) VaR
##'         2) (individual) error reached
##'         3) minimal [for worst VaR] or maximal [for best VaR] row sum
##'         4) number of oppositely ordered columns
##'         5) number of iterations through the matrix columns used for the last N
##' @author Marius Hofert
##' @note - We use "<= abs.err" to determine convergence instead of "< abs.err" as
##'         this then also nicely works with "= 0" (if abs.err=0) which stops in
##'         case the matrices are identical (no change at all).
##'       - No checking here due to speed
RA_aux <- function(X, method, err, maxiter, eps, impl)
{
    d <- ncol(X)
    if(impl=="C") {
        stop("Still in the works... Not available yet.") # TODO
        if(is.null(eps)) eps <- -1 # for C code
        if(is.infinite(maxiter)) maxiter <- -1 # for C code
        RA_aux_ <- NULL # to avoid "RA_aux: no visible binding for global variable 'RA_aux_'"
        .Call("RA_aux_", X, method, err, maxiter, eps)
    } else { # R implementation
          ## Define helper functions
          optim.fun <- if(method=="worst") min else max
          err.fun <- if(err=="absolute") function(x, y) abs(x-y) else function(x, y) abs((x-y)/y)
          ## Loop through the columns
          m.row.sums <- c()
          count <- 0
          while (TRUE) {
              ## Counter related quantities
              count <- count + 1 # increase counter
              mrs.old <- if(count == 1) optim.fun(rowSums(X)) else mrs.new # old min/max row sum
              ## Oppositely order X (=> Y)
              Y <- X
              for(j in 1:d)
                  Y[,j] <- sort(Y[,j], decreasing=TRUE)[rank(rowSums(Y[,-j]))]
              ## Compute minimal/maximal row sums
              mrs.new <- optim.fun(rowSums(Y)) # new min/max row sum
              m.row.sums <- c(m.row.sums, mrs.new) # append min/max row sum
              ## Check convergence (we use "<= eps" as it entails eps=0)
              stp <- (count == maxiter) || if(is.null(eps)) all(Y == X) else
              err.fun(mrs.new, mrs.old) <= eps
              if(stp) {
                  num.opp.ordered <- sum(sapply(seq_len(d), function(j)
                      all(sort(Y[,j], decreasing=TRUE)[rank(rowSums(Y[,-j]))] == Y[,j]))) # count number of oppositely ordered columns
                  individual.err <- err.fun(mrs.new, mrs.old) # compute the (individual) error
                  break
              } else X <- Y
          }
          ## Return
          list(bound=mrs.new, individual.err=individual.err, m.row.sums=m.row.sums,
               num.opp.ordered=num.opp.ordered, num.iter=count)
      }
}

##' @title Computing Lower/Upper Bounds for the Worst VaR with the RA
##' @param alpha confidence level
##' @param d dimension
##' @param qF a marginal quantile function (homogeneous case) or a d-list of such;
##'        note that each marginal quantile function has to be vectorized in p
##' @param N number of discretization points
##' @param abs.err absolute error to determine convergence; if NULL, then the
##'        iteration is done until the matrix doesn't change anymore, i.e., each
##'        column is oppositely ordered with respect to the sum of all others
##' @param maxiter maximal number of iterations
##' @param method character indicating which VaR is approximated (worst/best)
##' @param sample logical indicating whether each column of the two working
##'        matrices are sampled before the iteration begins
##' @param impl string indicating the implementation (either "C" or "R")
##' @return 5-list containing the
##'         1) computed lower and upper bounds for (worst or best) VaR
##'         2) the relative dependence uncertainty spread
##'            "|(upper bound - lower bound) / upper bound|"
##'         3) individual (for each bound) errors reached
##'         4) minimal [for worst VaR] or maximal [for best VaR] row sums for
##'            lower and upper bounds
##'         5) number of oppositely ordered columns for lower and upper bounds
##'         6) number of iterations through the matrix columns used
##' @author Marius Hofert
##' @note Notation is from p. 2757 in Embrechts, Puccetti, Rueschendorf (2013);
##'       variables are named according to the 'worst' VaR case.
RA <- function(alpha, d, qF, N, abs.err=NULL, maxiter=Inf,
               method=c("worst", "best"), sample=TRUE, impl=c("R", "C"))
{
    ## Checks and Step 1 (get N, abs.err)
    stopifnot(0 < alpha, alpha < 1, is.null(abs.err) || abs.err >= 0,
              length(N) >= 1, N >= 2, maxiter >= 1, is.logical(sample))
    method <- match.arg(method)
    impl <- match.arg(impl)
    if(missing(d))
        stopifnot(is.list(qF), sapply(qF, is.function), (d <- length(qF)) >= 2)
    else { # extend the given quantile function to a list of functions
        stopifnot(is.function(qF))
        qF <- rep(list(qF), d)
    }

    ## Compute lower bound

    ## Step 2 (build \underline{X}^\alpha)
    p <- if(method=="worst") alpha + (1-alpha)*0:(N-1)/N else alpha*0:(N-1)/N # N-vector of prob.
    X <- sapply(qF, function(qF) qF(p))
    ## adjust those that are -Inf (for method="best")
    ## use alpha*((0+1)/2 / N) = alpha/(2N) instead of 0 quantile
    if(method == "best")
        X[1,] <- sapply(1:d, function(j) if(is.infinite(X[1,j])) qF[[j]](alpha/(2*N)) else X[1,j])

    ## Step 3 (randomly permute each column of \underline{X}^\alpha)
    if(sample) X <- apply(X, 2, sample)

    ## Steps 4, 5 (determine \underline{X}^*)
    ## repeat oppositely ordering \underline{X}^\alpha until there is only an
    ## abs.err change in the min (method="worst") or max (method="best") row sum
    ## or until we reached maxiter number of iterations
    res.low <- RA_aux(X, method=method, err="absolute", maxiter=maxiter,
                      eps=abs.err, impl=impl)

    ## Compute upper bound

    ## Step 2 (build \overline{X}^\alpha)
    p <- if(method=="worst") alpha + (1-alpha)*1:N/N else alpha*1:N/N # N-vector of prob.
    X <- sapply(qF, function(qF) qF(p))
    ## adjust those that are Inf (for method="worst")
    ## use alpha+(1-alpha)*(N-1+N)/(2*N) = alpha+(1-alpha)*(1-1/(2*N)) instead of 1 quantile
    if(method == "worst") X[N,] <- sapply(1:d, function(j)
        if(is.infinite(X[N,j])) qF[[j]](alpha+(1-alpha)*(1-1/(2*N))) else X[N,j])

    ## Step 3 (randomly permute each column of \overline{X}^\alpha)
    if(sample) X <- apply(X, 2, sample)

    ## Step 6 (determine \overline{X}^*)
    ## repeat oppositely ordering \overline{X}^\alpha until there is only an
    ## abs.err change in the min (method="worst") or max (method="best") row sum
    ## or until we reached maxiter number of iterations
    res.up <- RA_aux(X, method=method, err="absolute", maxiter=maxiter,
                     eps=abs.err, impl=impl)

    ## Step 7 (return \underline{s}_N, \overline{s}_N and other info)
    list(bounds=c(res.low$bound, res.up$bound),
         rel.DU.spread=abs((res.up$bound-res.low$bound)/res.up$bound),
         individual.err=c(res.low$individual.err, res.up$individual.err),
	 num.iter=c(res.low$num.iter, res.up$num.iter),
         m.row.sums=list(res.low$m.row.sums, res.up$m.row.sums),
         num.opp.ordered=c(res.low$num.opp.ordered, res.up$num.opp.ordered))
}

##' @title Computing Lower/Upper Bounds for the Worst VaR with the ARA
##' @param alpha confidence level
##' @param d dimension
##' @param qF a marginal quantile function (homogeneous case) or a d-list of such;
##'        note that each marginal quantile function has to be vectorized in p
##' @param N vector of numbers of discretization points
##' @param rel.err 2-vector of relative errors for determining the individual
##'        relative error (i.e., the relative error in the minimal/maximal row sum
##'        for each of the bounds) and the joint relative error (i.e., the relative
##'        error between the computed lower and upper bounds). rel.err[1] can be
##'        NULL (=> iteration is done until matrices don't change anymore, i.e.,
##'        until each column is oppositely ordered to the sum of all others,
##'        but only if maxiter hasn't been reached)
##' @param maxiter maximal number of iterations per N
##' @param method character indicating which VaR is approximated (worst/best)
##' @param sample logical indicating whether each column of the two working
##'        matrices are sampled before the iteration begins
##' @param impl string indicating the implementation (either "C" or "R")
##' @return 5-list containing the
##'         1) computed lower and upper bounds for (worst or best) VaR
##'         2) the relative dependence uncertainty spread
##'            "|(upper bound - lower bound) / upper bound|"
##'         3) individual (for each bound) errors reached
##'         4) minimal [for worst VaR] or maximal [for best VaR] row sums for
##'            lower and upper bounds
##'         5) number of oppositely ordered columns for lower and upper bounds
##'         6) number of iterations through the matrix columns used
##' @author Marius Hofert
ARA <- function(alpha, d, qF, N=2^seq(8, 20, by=1), rel.err=c(0.001, 0.01),
                maxiter=10, method=c("worst", "best"), sample=TRUE,
                impl=c("R", "C"))
{
    ## Checks and Step 1 (get N, rel.err)
    stopifnot(0 < alpha, alpha < 1, length(rel.err) == 2,
              is.null(rel.err[1]) || rel.err[1] >=0, rel.err[2] >= 0,
              length(N) >= 1, N >= 2, maxiter >= 1, is.logical(sample))
    method <- match.arg(method)
    impl <- match.arg(impl)
    if(missing(d))
        stopifnot(is.list(qF), sapply(qF, is.function), (d <- length(qF)) >= 2)
    else { # extend the given quantile function to a list of functions
        stopifnot(is.function(qF))
        qF <- rep(list(qF), d)
    }

    ## Loop over N
    for(N. in N) {

        ## Compute lower bound

        ## Step 2 (build \underline{X}^\alpha)
        p <- if(method=="worst") alpha + (1-alpha)*0:(N.-1)/N. else alpha*0:(N.-1)/N. # N.-vector of prob.
        X <- sapply(qF, function(qF) qF(p))
        ## adjust those that are -Inf (for method="best")
        ## use alpha*((0+1)/2 / N.) = alpha/(2N.) instead of 0 quantile
        if(method == "best")
            X[1,] <- sapply(1:d, function(j) if(is.infinite(X[1,j])) qF[[j]](alpha/(2*N.)) else X[1,j])

        ## Step 3 (randomly permute each column of \underline{X}^\alpha)
        if(sample) X <- apply(X, 2, sample)

        ## Steps 4, 5 (determine \underline{X}^*)
        ## repeat oppositely ordering \underline{X}^\alpha until there is only an
        ## rel.err[1] change in the min (method="worst") or max (method="best") row sum
        ## or until we reached maxiter number of iterations
        res.low <- RA_aux(X, method=method, err="relative", maxiter=maxiter,
                          eps=rel.err[1], impl=impl)

        ## Compute upper bound

        ## Step 2 (build \overline{X}^\alpha)
        p <- if(method=="worst") alpha + (1-alpha)*1:N./N. else alpha*1:N./N. # N.-vector of prob.
        X <- sapply(qF, function(qF) qF(p))
        ## adjust those that are Inf (for method="worst")
        ## use alpha+(1-alpha)*(N.-1+N.)/(2*N.) = alpha+(1-alpha)*(1-1/(2*N.)) instead of 1 quantile
        if(method == "worst") X[N.,] <- sapply(1:d, function(j)
            if(is.infinite(X[N.,j])) qF[[j]](alpha+(1-alpha)*(1-1/(2*N.))) else X[N.,j])

        ## Step 3 (randomly permute each column of \overline{X}^\alpha)
        if(sample) X <- apply(X, 2, sample)

        ## Step 6 (determine \overline{X}^*)
        ## repeat oppositely ordering \overline{X}^\alpha until there is only an
        ## rel.err[1] change in the min (method="worst") or max (method="best") row sum
        ## or until we reached maxiter number of iterations
        res.up <- RA_aux(X, method=method, err="relative", maxiter=maxiter,
                         eps=rel.err[1], impl=impl)

        ## Determine convergence (individual + joint)
        joint.err <- abs((res.low$bound-res.up$bound)/res.up$bound)
        if( (max(res.low$individual.err, res.up$individual.err) <= rel.err[1]) &&
            joint.err <= rel.err[2] ) break

    }

    ## Step 7 (return \underline{s}_N, \overline{s}_N and other info)
    list(bounds=c(res.low$bound, res.up$bound),
         rel.DU.spread=abs((res.up$bound-res.low$bound)/res.up$bound),
         individual.err=c(res.low$individual.err, res.up$individual.err),
         joint.err=joint.err, N.used=N.,
         num.iter=c(res.low$num.iter, res.up$num.iter),
         m.row.sums=list(res.low$m.row.sums, res.up$m.row.sums),
         num.opp.ordered=c(res.low$num.opp.ordered, res.up$num.opp.ordered))
}

