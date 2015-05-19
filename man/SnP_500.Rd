\name{SnP_500}
\alias{SnP_500}
\docType{data}
\title{S&P 500 Data}
\description{
  Data set consisting of daily open, high, low and close prices, trading volume
  (number of shares traded that day) and adjusted closing prices of the
  S&P 500 from 1960-01-04 to 2013-12-31.
}
\usage{
data("SnP_500")
}
\format{
  A 6-column \code{xts} object (open, high, low, close, volume,
  adjusted) consisting of the corresponding \code{\link{numeric}} values
  of the S&P 500 from 1960-01-04 to 2013-12-31.
}
\source{
  The data was drawn from Yahoo! Finance (\url{http://finance.yahoo.com/}).
}
\examples{
data("SnP_500") # use as.data.frame(SnP_500) for converting it to a data.frame
head(SnP_500)
plot(SnP_500[,"close"], main="S&P 500 daily closing prices")
head(rownames(as.data.frame(SnP_500))) # row names; or use zoo's index()
}
\keyword{datasets}