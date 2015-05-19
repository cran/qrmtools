\name{worst_VaR}
\alias{crude_VaR_bounds}
\alias{worst_VaR_hom}
\alias{dual_bound}
\alias{RA}
\alias{ARA}
\title{Worst Value-at-Risk for Given Margins}
\description{
  Compute the worst Value-at-Risk for given marginal distributions under
  various setups.
}
\usage{
crude_VaR_bounds(alpha, d, qF, ...)
worst_VaR_hom(alpha, d, method=c("Wang", "Wang.Par", "dual"),
              interval=NULL, tol=NULL, ...)
dual_bound(s, d, pF, ...)
RA(alpha, d, qF, N, abs.err=NULL, maxiter=Inf, method=c("worst", "best"),
   sample=TRUE, impl=c("R", "C"))
ARA(alpha, d, qF, N=2^seq(8, 20, by=1), rel.err=c(0.001, 0.01), maxiter=10,
    method=c("worst", "best"), sample=TRUE, impl=c("R", "C"))
}
\arguments{
  \item{alpha}{Value-at-Risk confidence level (e.g., 0.99).}
  \item{d}{Dimension (number of risk factors).}
  \item{qF}{The marginal quantile function in the homogeneous case or,
    for \code{crude_VaR_bounds()} and \code{RA()}, a list (of length
    \code{d}) with the marginal quantile functions.}
  \item{method}{A \code{\link{character}}. For
    \describe{
      \item{\code{worst_VaR_hom()}:}{\code{method="Wang"}
	and \code{method="Wang.Par"} apply the approach of Embrechts,
	Puccetti, \enc{Rüschendorf}{Rueschendorf}, Wang and Beleraj (2014,
	Proposition 3.1); the latter assumes Pareto margins and thus does
	not require numerical integration. \code{method="dual"} applies the
	dual bound approach as in Embrechts, Puccetti,
	\enc{Rüschendorf}{Rueschendorf} (2013, Proposition 4).}
      \item{\code{RA()}, \code{ARA()}:}{\code{method} indicates whether bounds
	for the worst (i.e., largest) or the best (i.e., smallest)
	Value-at-Risk are computed.}
  }}
  \item{interval}{Initial interval (a \code{\link{numeric}(2)}). For
    \code{method="Wang"} and \code{method="Wang.Par"}, if not provided,
    the initivial interval is
    \eqn{[0,(1-\alpha)/d]}{[0,(1-alpha)/d]}. For \code{method="dual"},
    note that the lower endpoint of the initial interval has to be
    sufficiently large in order for the the inner root-finding algorithm
    to find a root; see Details.}
  \item{tol}{Tolerance for the involved root-finding algorithm (passed
    to \code{\link{uniroot}()}; defaults to \code{2.2204e–16}).
    For \code{method="dual"}, \code{tol} is
    used for the outer root-finding algorithm only (defaults to
    \code{.Machine$double.eps^0.25} as \code{\link{uniroot}()}).}
  \item{s}{Dual bound evaluation point.}
  \item{pF}{The marginal loss distribution function.}
  \item{N}{Number of discretization points. For
    \describe{
      \item{\code{RA()}:}{A positive integer.}
      \item{\code{ARA()}:}{A vector of positive integers. For each
	such number, at most \code{maxiter}-many iterations through the columns
	of the underlying matrix of quantiles are conducted, checking after
	each of these iterations whether the desired relative errors
	(both individually and jointly; see below) has been
	achieved.}
  }}
  \item{abs.err}{Absolute error \eqn{\epsilon}{epsilon} to determine
    the individual convergence (change in the computed minimal
    (for worst Value-at-Risk) or maximal (for best Value-at-Risk)
    row sums for the lower and upper bound). If \code{NULL} (the
    default) and \code{maxiter} has not been reached, the iterations
    through the columns of the underlying matrix of quantiles is done
    until the matrix does not change anymore, i.e., until each column
    is oppositely ordered to the sum of all others.}
  \item{rel.err}{A vector of length containing two relative errors,
    the individual relative error (first component; used to determine
    convergence of the minimal (for worst Value-at-Risk) or maximal
    (for best Value-at-Risk) row sums for the lower and upper bound)
    and the joint relative error (second component; relative error
    between the computed lower and upper bound computed with respect
    to the upper bound). Note that \code{rel.err[1]} is allowed to be
    \code{NULL}; see \code{abs.err} above for how to interpret this case.}
  \item{maxiter}{maximal number of iterations over the columns of the
    underlying matrix of quantiles (can be set to \code{Inf}).}
  \item{sample}{\code{\link{logical}} indicating whether each column of
    the two internal matrices of quantiles (see Step 3 of the Rearrangement
    Algorithm in Embrechts, Puccetti, \enc{Rüschendorf}{Rueschendorf} (2013))
    are sampled before the iteration begins.}
  \item{impl}{\code{\link{character}} providing the type of
    implementation used (C vs R).}
  \item{\dots}{For
    \describe{
      \item{\code{worst_VaR_hom()}:}{For \code{d=2} or
	\code{method="Wang"} the quantile function \code{qF()} (for the latter,
	additional arguments passed to the underlying \code{\link{integrate}()}
	have to be provided); for \code{method="dual"} the distribution
	function \code{pF()} and additionalarguments passed to the underlying
	\code{\link{integrate}()} for computing the dual bound
	\eqn{D(s)} have to be provided; and for \code{method="Wang.Par"}
	the Pareto parameter \eqn{\theta>0}{theta>0} has to be provided.}
      \item{\code{dual_bound()}:}{Additional arugments passed to the
	underlying \code{\link{integrate}()}.}
      \item{\code{crude_VaR_bounds()}:}{Additional arguments passed to
	(all provided) quantile functions.}
    }}
}
\value{
  \code{crude_VaR_bounds()} returns crude lower and upper bounds for
  Value-at-Risk at confidence level \eqn{\alpha}{alpha} for any
  \eqn{d}-dimensional model with marginal quantile functions
  specified by \code{qF}.

  \code{worst_VaR_hom()} returns the worst (i.e., largest) Value-at-Risk at
  confidence level \eqn{\alpha}{alpha} for \eqn{d} risks with equal
  distribution function specified by the ellipsis \code{...}.

  \code{dual_bound()} returns the value of the dual bound \eqn{D(s)} as
  given in Embrechts, Puccetti, \enc{Rüschendorf}{Rueschendorf}
  (2013, Eq. (12)).

  \code{RA()} returns a \code{\link{list}} containing
  \describe{
    \item{\code{bounds}:}{A bivariate vector containing the computed lower
      and upper bound for (worst or best) Value-at-Risk.}
    \item{\code{rel.DU.spread}:}{relative error between lower and upper
      bound (computed with respect to the upper bound).}
    \item{\code{individual.err}:}{A bivariate vector containing the
      individual errors (i.e., the absolute change of the minimal
      (for worst Value-at-Risk) or maximal (for best Value-at-Risk)
      row sums for the lower (first component) and upper bound (second
      component) in the last iteration).}
    \item{\code{num.iter}:}{A bivariate vector containing, for the lower
      and upper bound, the actual used number of iterations over the
      columns of the underlying matrix of quantiles.}
    \item{\code{m.row.sums}:}{A list of length two containing minimal
      (for worst Value-at-Risk) or maximal (for best Value-at-Risk)
      row sums for the lower and upper bound.}
    \item{\code{num.opp.ordered}:}{A bivariate vector containing
      the number of oppositely ordered columns for the lower and upper bound.}
  }

  \code{ARA()} returns a \code{\link{list}} containing
  \describe{
    \item{\code{bounds}:}{See \code{RA()}.}
    \item{\code{rel.DU.spread}:}{See \code{RA()}.}
    \item{\code{individual.err}:}{See \code{RA()}.}
    \item{\code{joint.err}:}{The joint relative error reached (i.e.,
      the relative error between the computed lower and upper bound
      computed with respect to the upper bound).}
    \item{\code{N.used}:}{The actual \code{N} used for computing
      the (final) lower and upper bound.}
    \item{\code{num.iter}:}{See \code{RA()}.}
    \item{\code{m.row.sums}:}{See \code{RA()}.}
    \item{\code{num.opp.ordered}:}{See \code{RA()}.}
  }
}
\details{
  For \code{d=2}, \code{worst_VaR_hom()} implements the method of
  Embrechts, Puccetti, \enc{Rüschendorf}{Rueschendorf} (2013,
  Proposition 2). For \code{method="Wang"} and \code{method="Wang.Par"},
  the method presented in Embrechts, Puccetti, \enc{Rüschendorf}{Rueschendorf},
  Wang, Beleraj (2014, Proposition 3.1) is implemented. This requires
  one \code{\link{uniroot}()} and, for the generic \code{method="Wang"},
  one \code{\link{integrate}()}. The critical part for these methods is
  the lower endpoint of the initial interval. If the (marginal)
  distribution function has finite first moment, this can be taken as
  0. However, if it has infinite first moment, the lower endpoint has to
  be positive (below the unknown root and smaller than the upper
  endpoint \eqn{(1-\alpha)/d}{(1-alpha)/d}). Also note that our default
  uses a higer tolerance than that of \code{\link{uniroot}()}, which is
  critical here. For \code{method="dual"}, the method presented of
  Embrechts, Puccetti, \enc{Rüschendorf}{Rueschendorf}
  (2013, Proposition 4) is implemented. This
  requires two (nested) \code{\link{uniroot}()}, and an
  \code{\link{integrate}()}. For the inner root-finding
  procedure to find a root, the lower endpoint of the provided initial
  \code{interval} has to be \dQuote{sufficiently large}.

  Note that both of these approaches for computing the worst (i.e., largest)
  Value-at-Risk in the homogeneous case are numerically non-trivial.
  See the source code and \code{demo(worst_VaR)} for more details.

  For the Rearrangement Algorithm \code{RA()}, unless
  \code{maxiter}-many iterations have been reached, convergence is
  determined for the lower and the upper bound if the minimal
  (for the worst Value-at-Risk) or maximal (for the best
  Value-at-Risk) row sums after an iteration did not change by more than
  (the absolute error) \code{abs.err} (so \eqn{\le\epsilon}{<= eps}). This
  is slightly different from Embrechts, Puccetti, \
  enc{Rüschendorf}{Rueschendorf} (2013) who use \eqn{<\epsilon}{< eps}
  but this slight modification has the advantage of allowing
  \eqn{\epsilon=0}{eps=0}. For the Adaptive Rearrangement
  Algorithm \code{ARA()}, again unless \code{maxiter}-many iterations
  have been reached, convergence is determined if the minimal
  (for the worst Value-at-Risk) or maximal (for the best
  Value-at-Risk) row sums after an iteration did not change by more than
  (the individual relative error) \code{rel.err[1]} \emph{and} the
  relative (joint) error between both bounds is at most \code{rel.err[2]}.

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
}
\author{Marius Hofert}
\references{
  Embrechts, P., Puccetti, G., \enc{Rüschendorf}{Rueschendorf}, L.,
  Wang, R., Beleraj, A. (2014).  An Academic Response to Basel
  3.5. \emph{Risks} \bold{2}(1), 25--48.

  Embrechts, P., Puccetti, G., \enc{Rüschendorf}{Rueschendorf}, L. (2013).
  Model uncertainty and VaR aggregation. \emph{Journal of Banking \&
    Finance} \bold{37}, 2750--2764.
}
\seealso{
  \code{demo(worst_VaR)} for more example calls, numerical challenges
  encoutered and a comparison of the different methods for computing
  the worst (i.e., largest) Value-at-Risk.
}
\examples{
## Setup
alpha <- 0.99 # VaR confidence level
th <- 2 # Pareto parameter
qF <- function(p) qPar(p, theta=th) # Pareto quantile function
pF <- function(q) pPar(q, theta=th) # Pareto distribution function

## d=2: Compute worst VaR ('method' is ignored)
stopifnot(all.equal(worst_VaR_hom(alpha, d=2, qF=qF), 26.28427, tolerance=1e-7))

## d=8: Compare the various methods in the homogeneous case (equal margins)
d <- 8 # dimension
I <- crude_VaR_bounds(alpha, d=d, qF=qF) # crude bound
## or: I <- crude_VaR_bounds(alpha, d=d, qF=qPar, theta=th)
(worst.VaR.Wang     <- worst_VaR_hom(alpha, d=d, method="Wang", qF=qF))
(worst.VaR.Wang.Par <- worst_VaR_hom(alpha, d=d, method="Wang.Par", theta=th))
(worst.VaR.dual     <- worst_VaR_hom(alpha, d=d, method="dual", interval=I, pF=pF))
stopifnot(all.equal(c(worst.VaR.Wang, worst.VaR.Wang.Par), rep(worst.VaR.dual, 2),
                    tolerance=1e-6))
## Compare with the Rearrangement Algorithm
set.seed(271) # for (the sampling) Step 3 in the Rearrangement Algorithm
RA1 <- RA(alpha, d=d, qF=qF, abs.err=0.01, N=1000)
(worst.VaR.RA <- RA1$bounds)
set.seed(271)
RA2 <- RA(alpha, d=d, qF=qF, N=1000)
(worst.VaR.RA <- RA2$bounds)
## Adaptive Rearrangement Algorithm (with defaults)
set.seed(271)
ARA1 <- ARA(alpha, d=d, qF=qF)
(worst.VaR.ARA <- ARA1$bounds)

## Using additional output of RA()
ran <- range(RA2$m.row.sums)
plot(RA2$m.row.sums[[2]], type="b", ylim=ran, xlab="Iteration",
     ylab=paste0("Minimal row sum (# of oppositely ordered columns: lower=",
     RA2$num.opp.ordered[1],", upper=",
     RA2$num.opp.ordered[2],")"),
     main=substitute("Rearrangement Algorithm for"~alpha==a.*","~d==d.*" and Par("*
                     th.*")", list(a.=alpha, d.=d, th.=th)))
lines(1:length(RA2$m.row.sums[[1]]), RA2$m.row.sums[[1]], type="b",
      col="blue")
legend("bottomright", inset=0.02, y.intersp=1.2, bty="n", lty=rep(1,2),
       col=c("black", "blue"), legend=c("upper bound", "lower bound"))

## "Reproducing" an example from Embrechts, Puccetti, Rueschendorf (2013)
## (seed unknown, eps unknown)
set.seed(271)
(RA3 <- RA(alpha, d=3, qF=qF, N=50))
stopifnot(all.equal(RA3$bounds[1], 45.11593, tolerance=1e-7),
          all.equal(RA3$bounds[2], 46.40413, tolerance=1e-7),
          RA3$num.opp.ordered == c(3,3))

## An example from Embrechts, Puccetti, Rueschendorf (2013) with larger N
## (first two decimals are equal)
if(FALSE) {
    RA4 <- RA(alpha, d=3, qF=qF, N=1e5) # ~ 15s
    stopifnot(all.equal(RA4$bounds, rep(45.98, 2), tolerance=3e-4))
}
}
\keyword{distribution}