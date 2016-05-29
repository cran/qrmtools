\name{VaR_bounds}
\alias{crude_VaR_bounds}
\alias{VaR_bounds_hom}
\alias{dual_bound}
\alias{rearrange}
\alias{RA}
\alias{ARA}
\title{Best and Worst Value-at-Risk for Given Margins}
\description{
  Compute the best and worst Value-at-Risk for given marginal distributions.
}
\usage{
crude_VaR_bounds(alpha, qF, ...)
VaR_bounds_hom(alpha, d, method=c("Wang", "Wang.Par", "dual"),
               interval=NULL, tol=NULL, ...)
dual_bound(s, d, pF, tol=.Machine$double.eps^0.25, ...)

rearrange(X, tol=0, tol.type=c("relative", "absolute"), max.ra=Inf,
          method=c("worst", "best"), sample=TRUE, is.sorted=FALSE, trace=FALSE)
RA(alpha, qF, N, abstol=0, max.ra=Inf, method=c("worst", "best"),
   sample=TRUE)
ARA(alpha, qF, N.exp=seq(8, 19, by=1), reltol=c(0, 0.01),
    max.ra=10*length(qF), method=c("worst", "best"), sample=TRUE)
}
\arguments{
  \item{alpha}{Value-at-Risk confidence level (e.g., 0.99).}
  \item{d}{Dimension (number of risk factors; \eqn{\ge 2}{>=2}).}
  \item{qF}{The marginal quantile function (in the homogeneous case) or
    a \code{d}-list containing the marginal quantile functions (in the
          general case).}
  \item{method}{A \code{\link{character}} string. For
    \describe{
      \item{\code{VaR_bounds_hom()}:}{\code{method="Wang"} and
	\code{method="Wang.Par"}
         apply the approach of McNeil et al. (2015,
	 Proposition 8.32) for computing best (i.e., smallest) and
	 worst (i.e., largest) Value-at-Risk. The
         latter method assumes Pareto margins and thus does
	 not require numerical integration.
	\code{method="dual"} applies the dual bound approach as
	 in Embrechts et al. (2013, Proposition 4)
         for computing worst Value-at-Risk (no value for the best
         Value-at-Risk can be obtained with this approach and thus
          \code{\link{NA}} is returned for the best Value-at-Risk).}
      \item{\code{rearrange()}, \code{RA()}, \code{ARA()}:}{
	\code{method} indicates whether bounds
	for the best or for the worst Value-at-Risk should be computed.
	These bounds are termed
        \eqn{\underline{s}_N} and \eqn{\overline{s}_N} in the
        literature (and below) and are theoretically not guaranteed
        bounds of (best or worst) Value-at-Risk;
	however, they are treated as such in practice and are typically in line
	with results from \code{VaR_bounds_hom()} in the homogeneous
	case, for example.}
  }}
  \item{interval}{Initial interval (a \code{\link{numeric}(2)}) for
    computing worst Value-at-Risk. If not provided, these are the defaults chosen:
    \describe{
      \item{\code{method="Wang"}:}{The initial
	interval is \eqn{[0,(1-\alpha)/d]}{[0,(1-alpha)/d]}.}
      \item{\code{method="Wang.Par"}:}{The initial
	interval is \eqn{[c_l,c_u]}, where \eqn{c_l} and \eqn{c_u}
	are chosen as in Hofert et al. (2015).}
      \item{\code{method="dual"}:}{In this case, no good defaults are known.
	Note that the lower endpoint of the initial interval has to be
	sufficiently large in order for the the inner root-finding algorithm
	to find a root; see Details.}
  }}
  \item{tol}{
    \describe{
      \item{\code{VaR_bounds_hom()}:}{The tolerance for \code{\link{uniroot}()}
	for computing worst Value-at-Risk. This defaults (for \code{tol=NULL}) to
	\eqn{2.2204*10^{-16}}{2.2204*10^{-16}}
	for \code{method="Wang"} or \code{method="Wang.Par"} (where a
	smaller tolerance is crucial) and to \code{\link{uniroot}()}'s
	default \code{.Machine$double.eps^0.25} otherwise. Note that for
	\code{method="dual"}, \code{tol} is used for both the outer
	and the inner root-finding procedure.}
      \item{\code{rearrange()}:}{(Absolute or relative) tolerance to
	determine (the individual) convergence. This should normally be
        a number greater than or equal to 0, but we also allow \code{tol=NULL}
        which rearranges the columns until all columns are oppositely
        ordered to the sum of all other columns (this should only be used by
        experts).}
  }}
  \item{tol.type}{\code{\link{character}} string indicating the
    type of convergence tolerance function to be used (\code{"relative"}
    for relative tolerance and \code{"absolute"} for absolute tolerance).}
  \item{s}{Dual bound evaluation point.}
  \item{pF}{The marginal loss distribution function.}
  \item{X}{An (\code{N}, \code{d})-matrix of quantiles (to be
    rearranged). If \code{is.sorted} it is assumed that the columns of
    \code{X} are sorted in \emph{increasing} order.}
  \item{max.ra}{Maximal number of (considered) column rearrangements
    of the underlying matrix of quantiles (can be set to \code{Inf}).}
  \item{N}{The number of discretization points.}
  \item{N.exp}{The exponents of the number of discretization points
    (a \code{\link{vector}}) over which the algorithm iterates to find
    the smallest number of discretization points for which the desired
    accuracy (specified by \code{reltol}) is attained; for each number
    of discretization points, at most \code{max.ra}-many
    column rearrangements are of the underlying matrix of quantiles
    are considered.}
  \item{abstol}{Absolute convergence tolerance \eqn{\epsilon}{epsilon}
    to determine the individual convergence, i.e., the change in the computed minimal
    (for \code{method="worst"}) or maximal (for \code{method="best"})
    row sums for the lower bound \eqn{\underline{s}_N} and the upper
    bound \eqn{\overline{s}_N}. \code{abstol} is typically
    \eqn{\ge0}{>=0}; it can also be \code{\link{NULL}}, see \code{tol} above.}
  \item{reltol}{A \code{\link{vector}} of length one or two containing the
    relative convergence tolerances. If \code{reltol} is of length two,
    the first component is the the individual relative tolerance
    (used to determine
    convergence of the minimal (for \code{method="worst"}) or maximal
    (for \code{method="best"}) row sums for
    \eqn{\underline{s}_N} and \eqn{\overline{s}_N})
    and the second component is the joint relative tolerance (i.e., the
    relative tolerance between the computed \eqn{\underline{s}_N} and
    \eqn{\overline{s}_N} with respect to
    \eqn{\overline{s}_N}). If \code{reltol} is of length one, it denotes
    the joint relative tolerance; the individual relative tolerance is
    then taken as \code{NULL} (see \code{tol} above).}
  \item{sample}{A \code{\link{logical}} indicating whether each column of
    the two underlying matrices of quantiles (see Step 3 of the Rearrangement
    Algorithm in Embrechts et al. (2013))
    are randomly permuted before the rearrangements begin. This typically has
    quite a positive effect on run time (as most of the time is spent (oppositely)
    ordering columns).}
  \item{is.sorted}{A \code{\link{logical}} indicating whether the columns of \code{X}
    are sorted in increasing order.}
  \item{trace}{A \code{\link{logical}} indicating whether the underlying matrix is
    printed after each rearrangement step. See \code{vignette("VaR_bounds",
    package="qrmtools")} for how to interpret the output.}
  \item{\dots}{
    \describe{
      \item{\code{crude_VaR_bounds()}:}{Ellipsis argument passed
        to (all provided) quantile functions.}
      \item{\code{VaR_bounds_hom()}:}{The case \code{method="Wang"}
        requires the quantile function \code{qF()} to be provided
        and additional arguments passed via the ellipsis argument are passed on to
	% the underlying objective function Wang_h() from which they
	% are passed on to
	the underlying \code{\link{integrate}()}. For
        \code{method="Wang.Par"}
	the ellipsis argument must contain the parameter \eqn{\theta>0}{theta>0} of
	the Pareto distribution.
	For \code{method="dual"}, the ellipsis argument must contain the distribution
	function \code{pF()} and the initial interval \code{interval} for
	the outer root finding procedure (not for \code{d=2});
	additional arguments are passed on to the underlying
	\code{\link{integrate}()} for computing the dual bound
	\eqn{D(s)}.}
      \item{\code{dual_bound()}:}{The ellipsis argument is passed to the
	underlying \code{\link{integrate}()}.}
    }}
}
\value{
  \code{crude_VaR_bounds()} returns crude lower and upper bounds for
  Value-at-Risk at confidence level \eqn{\alpha}{alpha} for any
  \eqn{d}-dimensional model with marginal quantile functions
  specified by \code{qF}.

  \code{VaR_bounds_hom()} returns the best and worst Value-at-Risk at
  confidence level \eqn{\alpha}{alpha} for \eqn{d} risks with equal
  distribution function specified by the ellipsis \code{...}.

  \code{dual_bound()} returns the value of the dual bound \eqn{D(s)} as
  given in Embrechts, Puccetti, \enc{Rüschendorf}{Rueschendorf}
  (2013, Eq. (12)).

  \code{rearrange()} returns a \code{\link{list}} containing
  \describe{
    \item{\code{bound}:}{The computed \eqn{\underline{s}_N}
      or \eqn{\overline{s}_N}.}
    \item{\code{tol}:}{The reached tolerance (i.e., the (absolute or
      relative) change of the minimal (for \code{method="worst"}) or
      maximal (for \code{method="best"}) row sum after the last rearranged
      column).}
    \item{\code{converged}:}{A \code{\link{logical}} indicating whether
      the desired (absolute or relative) tolerance \code{tol} has been reached.}
    \item{\code{m.row.sums}:}{A \code{\link{vector}} containing the
      computed minimal (for \code{method="worst"}) or maximal
      (for \code{method="best"}) row sums after each (considered)
      column rearrangement.}
    \item{\code{X.rearranged}:}{An (\code{N},
      \code{d})-\code{\link{matrix}}
      containing the rearranged \code{X}.}
  }

  \code{RA()} returns a \code{\link{list}} containing
  \describe{
    \item{\code{bounds}:}{A bivariate vector containing the computed
      \eqn{\underline{s}_N} and \eqn{\overline{s}_N} (the so-called
      rearrangement range) which are typically treated as bounds for
      (worst or best) Value-at-Risk; see also above.}
    \item{\code{rel.ra.gap}:}{The reached relative tolerance (also known as
      relative rearrangement gap) between
      \eqn{\underline{s}_N} and \eqn{\overline{s}_N} computed with
      respect to \eqn{\overline{s}_N}.}
    \item{\code{ind.abs.tol}:}{A bivariate \code{\link{vector}} containing the
      reached individual absolute tolerances (i.e., the absolute change
      of the minimal
      (for \code{method="worst"}) or maximal (for \code{method="best"})
      row sums for computing \eqn{\underline{s}_N} and \eqn{\overline{s}_N};
      see also \code{tol} returned by \code{rearrange()} above).}
    \item{\code{converged}:}{A bivariate \code{\link{logical}} vector
      indicating convergence of the computed \eqn{\underline{s}_N} and
      \eqn{\overline{s}_N} (i.e., whether the desired tolerances were
      reached).}
    \item{\code{num.ra}:}{A bivariate vector containing the number
      of column rearrangments of the underlying matrices
      of quantiles for
      \eqn{\underline{s}_N} and \eqn{\overline{s}_N}.}
    \item{\code{m.row.sums}:}{A \code{\link{list}} of length two containing minimal
      (for \code{method="worst"}) or maximal (for \code{method="best"})
      row sums (after each (considered) column rearrangement)
      for the computed \eqn{\underline{s}_N} and \eqn{\overline{s}_N}.}
    \item{\code{X}:}{The initially constructed (\code{N}, \code{d})-matrices
      of quantiles for computing
      \eqn{\underline{s}_N} and \eqn{\overline{s}_N}.}
    \item{\code{X.rearranged}:}{The rearranged matrices \code{X} (for
      \eqn{\underline{s}_N} and \eqn{\overline{s}_N}).}
  }

  \code{ARA()} returns a \code{\link{list}} containing
  \describe{
    \item{\code{bounds}:}{See \code{RA()}.}
    \item{\code{rel.ra.gap}:}{See \code{RA()}.}
    \item{\code{rel.tol}:}{A trivariate \code{\link{vector}} containing
      the reached individual relative tolerances and the reached joint
      relative tolerance (computed with respect to \eqn{\overline{s}_N}).}
    \item{\code{converged}:}{A trivariate \code{\link{logical}} \code{\link{vector}}
      indicating individual convergence of the computed
      \eqn{\underline{s}_N} (first entry) and \eqn{\overline{s}_N} (second entry) and
      indicating joint convergence of the two bounds according to the attained
      joint relative tolerance (third entry).}
    \item{\code{N.used}:}{The actual \code{N} used for computing
      the (final) \eqn{\underline{s}_N} and \eqn{\overline{s}_N}.}
    \item{\code{num.ra}:}{See \code{RA()}; computed for \code{N.used}.}
    \item{\code{m.row.sums}:}{See \code{RA()}; computed for \code{N.used}.}
    \item{\code{X}:}{See \code{RA()}; computed for for \code{N.used}.}
    \item{\code{X.rearranged}:}{See \code{RA()}; computed for for \code{N.used}.}
  }
}
\details{
  For \code{d=2}, \code{VaR_bounds_hom()} uses the method of
  Embrechts et al. (2013,
  Proposition 2). For \code{method="Wang"} and \code{method="Wang.Par"}
  the method presented in McNeil et al. (2015, Prop. 8.32) is
  implemented; this goes back to Embrechts et al. (2014, Prop. 3.1; note that
  the published version of this paper contains typos for both bounds).
  This requires
  one \code{\link{uniroot}()} and, for the generic \code{method="Wang"},
  one \code{\link{integrate}()}. The critical part for the
  generic \code{method="Wang"} is the lower endpoint of the initial
  interval for \code{\link{uniroot}()}. If the (marginal)
  distribution function has finite first moment, this can be taken as
  0. However, if it has infinite first moment, the lower endpoint has to
  be positive (but must lie below the unknown root). Note that the upper
  endpoint \eqn{(1-\alpha)/d}{(1-alpha)/d} also happens to be a
  root and thus one needs a proper initional interval containing the
  root and being stricticly contained in
  \eqn{(0,(1-\alpha)/d}{(1-alpha)/d)}.
  In the case of Pareto margins, Hofert et al. (2015) have
  derived such an initial (which is used by
  \code{method="Wang.Par"}).
  Also note that the chosen smaller default tolerances for
  \code{\link{uniroot}()} in case of \code{method="Wang"} and
  \code{method="Wang.Par"} are crucial for obtaining reliable
  Value-at-Risk values; see Hofert et al. (2015).

  For \code{method="dual"} for computing worst Value-at-Risk, the method
  presented of Embrechts et al. (2013, Proposition 4) is implemented.
  This requires two (nested) \code{\link{uniroot}()}, and an
  \code{\link{integrate}()}. For the inner root-finding procedure to
  find a root, the lower endpoint of the provided initial
  \code{interval} has to be \dQuote{sufficiently large}.

  Note that these approaches for computing the
  Value-at-Risk bounds in the homogeneous case are numerically non-trivial;
  see the source code and \code{vignette("VaR_bounds",
    package="qrmtools")}
  for more details. As a
  rule of thumb, use \code{method="Wang"} if you have to (i.e., if the
  margins are not Pareto) and \code{method="Wang.Par"} if you can (i.e.,
  if the margins are Pareto). It is not recommended to use
  (the numerically even more challenging) \code{method="dual"}.


  Concerning the inhomogeneous case, \code{rearrange()} is an auxiliary
  function which is called by \code{RA()} and \code{ARA()}.
  After a column of \code{X} has been rearranged, the tolerance between
  the minimal (for the worst Value-at-Risk) or maximal (for the best
  Value-at-Risk) row sum after this rearrangement and the
  one of \eqn{d} steps before (so typically when that column
  was rearranged the last time) is computed and convergence determined.
  For performance reasons, no checking is done and \code{rearrange()}
  can change in future versions to (futher) improve run time. Overall
  it should only be used by experts.

  For the Rearrangement Algorithm \code{RA()}, convergence of
  \eqn{\underline{s}_N} and \eqn{\overline{s}_N} is determined if the
  minimal (for the worst Value-at-Risk) or maximal (for the best
  Value-at-Risk) row sum satisfies the specified \code{abstol}
  (so \eqn{\le\epsilon}{<=eps})
  after at most \code{max.ra}-many column rearrangements. This is different
  from Embrechts et al. (2013) who use \eqn{<\epsilon}{< eps} and
  only check for convergence after an iteration through all
  columns of the underlying matrix of quantiles has been completed.

  For the Adaptive Rearrangement Algorithm \code{ARA()},
  convergence of \eqn{\underline{s}_N} and \eqn{\overline{s}_N}
  is determined if, after at most \code{max.ra}-many column
  rearrangements, the (the individual relative tolerance)
  \code{reltol[1]} is satisfied \emph{and} the
  relative (joint) tolerance between both bounds is at most \code{reltol[2]}.

  Note that both \code{RA()} and \code{ARA()} need to evalute the
  0-quantile (for the lower bound for the best Value-at-Risk) and
  the 1-quantile (for the upper bound for the
  worst Value-at-Risk). As the algorithms can only handle finite values, the
  0-quantile and the 1-quantile need to be adjusted if infinite. Instead
  of the 0-quantile, the \eqn{\alpha/(2N)}{alpha/(2N)}-quantile is
  computed and instead of the 1-quantile the
  \eqn{\alpha+(1-\alpha)(1-1/(2N))}{alpha+(1-alpha)(1-1/(2N))}-quantile
  is computed for such margins (if the 0-quantile or the 1-quantile is
  finite, no adjustment is made).

  As a rule of thumb (see the examples below,
  \code{vignette("VaR_bounds", package="qrmtools")}
  and Hofert et al. (2015) for the reasons),
  it is recommended to use \code{ARA()} instead of \code{RA()}.

  On the theoretical side, let us again stress the following.
  \code{rearrange()}, \code{RA()} and \code{ARA()} compute
  \eqn{\underline{s}_N} and \eqn{\overline{s}_N} which are, from a
  practical point of view, treated as bounds for the worst
  (i.e., largest) or the best (i.e., smallest) Value-at-Risk
  (whatever is chosen with \code{method}), but which are not
  known to be such bounds from a theoretical point of view; see also
  above. Calling them \dQuote{bounds} for worst or best Value-at-Risk is
  thus theoretically not correct (unless proven) but \dQuote{practical}.
  The literature thus speaks of \eqn{(\underline{s}_N, \overline{s}_N)}
  as the rearrangement range (rather than an interval containing worst
  or best Value-at-Risk).
}
\author{Marius Hofert}
\references{
  Embrechts, P., Puccetti, G., \enc{Rüschendorf}{Rueschendorf}, L.,
  Wang, R., Beleraj, A. (2014).  An Academic Response to Basel
  3.5. \emph{Risks} \bold{2}(1), 25--48.

  Embrechts, P., Puccetti, G., \enc{Rüschendorf}{Rueschendorf}, L. (2013).
  Model uncertainty and VaR aggregation. \emph{Journal of Banking \&
    Finance} \bold{37}, 2750--2764.

  McNeil, A. J., Frey, R., and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.

  Hofert, M., Memartoluie, A., Saunders, D., Wirjanto, T. (2015).
  Improved Algorithms for Computing Worst Value-at-Risk:
  Numerical Challenges and the Adaptive Rearrangement Algorithm.
  See http://arxiv.org/abs/1505.02281.
}
\seealso{
  \code{vignette("VaR_bounds", package="qrmtools")}
  for more example calls, numerical challenges
  encoutered and a comparison of the different methods for computing
  the worst (i.e., largest) Value-at-Risk.
}
\examples{
## Pareto setup
alpha <- 0.99 # VaR confidence level
th <- 2 # Pareto parameter theta
qF <- function(p, theta=th) qPar(p, theta=theta) # Pareto quantile function
pF <- function(q, theta=th) pPar(q, theta=theta) # Pareto distribution function


## d=2: Compute best/worst VaR explicitly (hom. case) and compare with (A)RA ###

d <- 2 # dimension

## Explicit
VaRbounds <- VaR_bounds_hom(alpha, d=d, qF=qF) # (best VaR, worst VaR)

## Adaptive Rearrangement Algorithm (ARA)
set.seed(271) # set seed (for reproducibility)
ARAbest  <- ARA(alpha, qF=rep(list(qF), d), method="best")
ARAworst <- ARA(alpha, qF=rep(list(qF), d))

## Rearrangement Algorithm (RA) with N as in ARA()
RAbest  <- RA(alpha, qF=rep(list(qF), d), N=ARAbest$N.used, method="best")
RAworst <- RA(alpha, qF=rep(list(qF), d), N=ARAworst$N.used)

## Compare
stopifnot(all.equal(c(ARAbest$bounds[1], ARAbest$bounds[2],
                      RAbest$bounds[1],  RAbest$bounds[2]),
                    rep(VaRbounds[1], 4), tolerance=0.004, check.names=FALSE))
stopifnot(all.equal(c(ARAworst$bounds[1], ARAworst$bounds[2],
                      RAworst$bounds[1],  RAworst$bounds[2]),
                    rep(VaRbounds[2], 4), tolerance=0.003, check.names=FALSE))


## d=8: Compute best/worst VaR (hom. case) and compare with (A)RA ##############

d <- 8 # dimension

## Compute VaR bounds with various methods
I <- crude_VaR_bounds(alpha, qF=rep(list(qF), d)) # crude bound
VaR.W     <- VaR_bounds_hom(alpha, d=d, method="Wang", qF=qF)
VaR.W.Par <- VaR_bounds_hom(alpha, d=d, method="Wang.Par", theta=th)
VaR.dual  <- VaR_bounds_hom(alpha, d=d, method="dual", interval=I, pF=pF)

## Adaptive Rearrangement Algorithm (ARA) (with different relative tolerances)
set.seed(271) # set seed (for reproducibility)
ARAbest  <- ARA(alpha, qF=rep(list(qF), d), reltol=c(0.001, 0.01), method="best")
ARAworst <- ARA(alpha, qF=rep(list(qF), d), reltol=c(0.001, 0.01))

## Rearrangement Algorithm (RA) with N as in ARA and abstol (roughly) chosen as in ARA
RAbest  <- RA(alpha, qF=rep(list(qF), d), N=ARAbest$N.used,
              abstol=mean(tail(abs(diff(ARAbest$m.row.sums$low)), n=1),
                          tail(abs(diff(ARAbest$m.row.sums$up)), n=1)),
              method="best")
RAworst <- RA(alpha, qF=rep(list(qF), d), N=ARAworst$N.used,
              abstol=mean(tail(abs(diff(ARAworst$m.row.sums$low)), n=1),
                          tail(abs(diff(ARAworst$m.row.sums$up)), n=1)))

## Compare
stopifnot(all.equal(c(VaR.W[1], ARAbest$bounds, RAbest$bounds),
                    rep(VaR.W.Par[1],5), tolerance=0.004, check.names=FALSE))
stopifnot(all.equal(c(VaR.W[2], VaR.dual[2], ARAworst$bounds, RAworst$bounds),
                    rep(VaR.W.Par[2],6), tolerance=0.003, check.names=FALSE))


## Using (some of) the additional results computed by (A)RA() ##################

xlim <- c(1, max(sapply(RAworst$m.row.sums, length)))
ylim <- range(RAworst$m.row.sums)
plot(RAworst$m.row.sums[[2]], type="l", xlim=xlim, ylim=ylim,
     xlab="Number or rearranged columns",
     ylab=paste0("Minimal row sum per rearranged column"),
     main=substitute("Worst VaR minimal row sums ("*alpha==a.*","~d==d.*" and Par("*
                     th.*"))", list(a.=alpha, d.=d, th.=th)))
lines(1:length(RAworst$m.row.sums[[1]]), RAworst$m.row.sums[[1]], col="royalblue3")
legend("bottomright", bty="n", lty=rep(1,2),
       col=c("black", "royalblue3"), legend=c("upper bound", "lower bound"))
## => One should use ARA() instead of RA()


## "Reproducing" examples from Embrechts et al. (2013) #########################

## "Reproducing" Table 1 (but seed and eps are unknown)
## Left-hand side of Table 1
N <- 50
qPar <- rep(list(qF), 3)
p <- alpha + (1-alpha)*(0:(N-1))/N # for 'worst' (= largest) VaR
X <- sapply(qPar, function(qF) qF(p))
cbind(X, rowSums(X))
## Right-hand side of Table 1
set.seed(271)
res <- RA(alpha, qF=qPar, N=N)
row.sum <- rowSums(res$X.rearranged$low)
cbind(res$X.rearranged$low, row.sum)[order(row.sum),]

## "Reproducing" Table 3 for alpha=0.99 (but seed is unknown)
N <- 2e4 # we use a smaller N here to save run time
eps <- 0.1 # absolute tolerance
xi <- c(1.19, 1.17, 1.01, 1.39, 1.23, 1.22, 0.85, 0.98)
beta <- c(774, 254, 233, 412, 107, 243, 314, 124)
qF.lst <- lapply(1:8, function(j){ function(p) qGPD(p, xi=xi[j], beta=beta[j])})
set.seed(271)
res.best <- RA(0.99, qF=qF.lst, N=N, abstol=eps, method="best")
print(format(res.best$bounds, scientific=TRUE), quote=FALSE) # close to first value of 1st row
res.worst <- RA(0.99, qF=qF.lst, N=N, abstol=eps)
print(format(res.worst$bounds, scientific=TRUE), quote=FALSE) # close to last value of 1st row
}
\keyword{programming}