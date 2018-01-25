\name{nested_matrix}
\alias{nested_matrix}
\title{Construction of Nested Block Matrices}
\description{
  Constructing nested block matrices, used, for example, for
  hierarchical dependence models.
}
\usage{
nested_matrix(x, diagonal = rep(1, d))
}
\arguments{
  \item{x}{\code{\link{list}} of length 2 or 3 containing the
    homogeneous \code{\link{numeric}}
    entry of the current block of the nested block matrix, the
    \code{\link{integer}} components belongning to the current block (or
    \code{NULL}) and, possibly, another (nested) \code{\link{list}} of
    the same type.}
  \item{diagonal}{diagonal elements of the nested block matrix.}
}
\value{
  A nested block \code{\link{matrix}} of the structure
  as specified in \code{x} with off-diagonal entries as specified
  in \code{x} and diagonal entries as specified in \code{diagona}.
}
\details{
  See the examples for how to use.
}
\author{Marius Hofert}
\examples{
rho <- c(0.2, 0.3, 0.5, 0.8) # some entries (e.g., correlations)

## Test homogeneous case
x <- list(rho[1], 1:6)
nested_matrix(x)

## Two-level case with one block of size 2
x <- list(rho[1], 1, list(rho[2], 2:3))
nested_matrix(x)

## Two-level case with one block of size 2 and a larger nonnested block
x <- list(rho[1], 1:3, list(rho[2], 4:5))
nested_matrix(x)

## Test two-level case with three blocks of size 2
x <- list(rho[1], NULL, list(list(rho[2], 1:2),
                             list(rho[3], 3:4),
                             list(rho[4], 5:6)))
nested_matrix(x)

## Test three-level case
x <- list(rho[1], 1:3, list(rho[2], NULL, list(list(rho[3], 4:5),
                                               list(rho[4], 6:8))))
nested_matrix(x)

## Test another three-level case
x <- list(rho[1], c(3, 6, 1), list(rho[2], c(9, 2, 7, 5), list(rho[3], c(8, 4))))
nested_matrix(x)
}
\keyword{utilities}