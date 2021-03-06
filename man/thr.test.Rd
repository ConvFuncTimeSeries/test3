% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Hpo_test.r
\name{thr.test}
\alias{thr.test}
\title{Threshold Nonlinearity Test}
\usage{
thr.test(y, p = 1, d = 1, thrV = NULL, ini = 40, include.mean = T)
}
\arguments{
\item{y}{a time seris.}

\item{p}{AR order.}

\item{d}{delay for the threhosld variable.}

\item{thrV}{threshold variable.}

\item{ini}{inital number of data to start RLS estimation.}

\item{include.mean}{a logical value for including constant terms.}
}
\value{
\code{thr.test} returns a list with components: 
\item{F-ratio}{F statistic.}
\item{df}{the numerator and denominator degrees of freedom.}
\item{ini}{initial number of data to start RLS estimation.}
}
\description{
Threshold nonlinearity test.
}
\examples{
arorder=rep(1,2)
ar.coef=matrix(c(0.7,-0.8),2,1)
y=uTAR.sim(100,arorder,ar.coef,1,0)
thr.test(y$series,1,1,y$series,40,TRUE)
}
\references{
Tsay, R. (1989) Testing and Modeling Threshold Autoregressive Processes. \emph{Journal of the American Statistical Associations} \strong{84}(405), 231-240.
}
