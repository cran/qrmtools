\name{returns}
\alias{log_returns}
\title{Compute Log-Returns or the Inverse Transformation}
\description{
  Compute log-returns or data from given
  log-returns (the inverse transformation).
}
\usage{
log_returns(x, inverse = FALSE, start.value, drop = TRUE)
}
\arguments{
  \item{x}{A matrix of values to be turned into log-returns (if \code{inverse = FALSE})
    or log-returns (if \code{inverse = TRUE}).}
  \item{inverse}{A \code{\link{logical}} indicating whether the inverse transformation
    (data from given log-returns) is to be computed (if TRUE, this
    requires \code{start.value} to be specified).}
  \item{start.value}{If \code{inverse = TRUE}, the last available value
    of the time series to be constructed from the given log-returns \code{x}.}
  \item{drop}{A \code{\link{logical}} indicating whether 1-column matrices (or vectors)
   are returned as vectors.}
}
\value{
  A \code{\link{matrix}} with the same number of columns as \code{x}
  (if one, then a \code{\link{vector}}), just one row less (if \code{inverse = FALSE})
  or one row more (if \code{inverse = TRUE}).
}
\details{
  For negative log-returns, use -\code{log_returns(x)} or
  \code{log_returns(-x, inverse = TRUE, start.value = ...)}.
}
\author{Marius Hofert}
\examples{
## Generate some data
n <- 1000
set.seed(271)
S <- cumsum(sample(c(-1,1), size = n, replace = TRUE))
S <- if(any(S <= 0)) S <- 1 - min(S) + S

## Log-returns
X <- log_returns(S) # build log-returns
Y <- log_returns(X, inverse = TRUE, start.value = S[1]) # transform back (first date is lost)
stopifnot(all.equal(S, Y))

## -Log-returns
X <- -log_returns(S) # build -log-returns
Y <- log_returns(-X, inverse = TRUE, start.value = S[1]) # transform back
stopifnot(all.equal(S, Y))
}
\keyword{utilities}