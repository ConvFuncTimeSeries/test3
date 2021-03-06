% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tvAR.r
\name{rcAR}
\alias{rcAR}
\title{Estimating of Random-Coefficient AR Models}
\usage{
rcAR(x, lags = c(1), include.mean = TRUE)
}
\arguments{
\item{x}{a time series of data.}

\item{lags}{the lag of AR models. This is more flexible than using order. It can skip unnecessary lags.}

\item{include.mean}{a logical value indicating whether the constant terms are included.}
}
\value{
\code{rcAR} function returns a list with following components:
\item{par}{estimated parameters.}
\item{se.est}{standard errors.}
\item{residuals}{residuals.}
\item{sresiduals}{standardized residuals.}
}
\description{
Estimate random-coefficient AR models.
}
\examples{
x1=rnorm(100)
x2=arima.sim(n = 63, list(ar = c(0.8897, -0.4858)),sd = sqrt(0.1796))
x=c(x1,x2)
est=rcAR(x,2,TRUE)
}
