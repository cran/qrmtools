\name{fire}
\alias{fire}
\docType{data}
\title{Danish Fire Insurance Claims}
\description{
  Fire insurance claims in 1M DKK in Denmark from 1980-01-03 to 1990-12-31.
}
\usage{
data("fire")
}
\format{
  \code{\link{xts}} object containing 2167 observations.
}
\source{
  originally Mette Rytgaard (Copenhagen Re); see \pkg{evir} and
  \pkg{QRM}.
%% Created via (after comparing with evir):
%% library(QRM)
%% data(danish)
%% fire <- danish
%% names(fire) <- "Loss"
%% library(xts)
%% fire <- as.xts(fire)
%% save(fire, file = "fire.rda", compress = "xz")
}
\author{Marius Hofert}
\examples{
library(xts)
data("fire")
str(fire)
stopifnot(inherits(fire, "xts"), length(fire) == 2167)
plot.zoo(fire, ylab = "Fire insurance claim")
}
\keyword{datasets}
