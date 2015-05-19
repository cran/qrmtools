\name{Black_Scholes}
\alias{Black_Scholes}
\alias{Black_Scholes_Greeks}
\title{Black--Scholes formula and the Greeks}
\description{
  Compute the Black--Scholes formula and the Greeks.
}
\usage{
Black_Scholes(t, S, r, sigma, K, T, type=c("call", "put"))
Black_Scholes_Greeks(t, S, r, sigma, K, T)
}
\arguments{
  \item{t}{initial or current time \eqn{t} (in years).}
  \item{S}{stock price at time \eqn{t}.}
  \item{r}{risk-free annual interest rate.}
  \item{sigma}{annual volatility (standard deviation).}
  \item{K}{strike.}
  \item{T}{maturity (in years).}
  \item{type}{\code{\link{character}} string indicating whether the price
    of a call (the default) or of put option is to be computed.}
}
\value{
  \code{Black_Scholes()} returns the value of a European call or put
  option (depending on the chosen \code{type}).

  \code{Black_Scholes_Greeks()} returns the first-order derivatives
  delta, theta, rho, vega and the second-order derivatives gamma, vanna
  and vomma (in this order).
}
%% \details{
%% }
\author{Marius Hofert}
\references{
   McNeil, A. J., Frey, R., and Embrechts, P. (2015).
   \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
   Princeton University Press.
}
%% \examples{
%% }
\keyword{distribution}