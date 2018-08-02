#' Estimate Time-Varying Coefficient AR Models
#'
#' Estimate time-varying coefficient AR models.
#' @param x a time series of data.
#' @param lags the lagged variables used, e.g. lags=c(1,3) means lag-1 and lag-3 are used as regressors. 
#' It is more flexible than specifying an order.
#' @param include.mean a logical value indicating whether the constant terms are included. 
#' @return \code{trAR} function returns the value from function \code{dlmMLE}.
#' @examples
#' x1=rnorm(100)
#' x2=arima.sim(n = 63, list(ar = c(0.8897, -0.4858)),sd = sqrt(0.1796))
#' x=c(x1,x2)
#' est=tvAR(x,2)
#' @export 
"tvAR" <- function(x,lags=c(1),include.mean=TRUE){
if(is.matrix(x))x <- c(x[,1])
nlag <-  length(lags)
if(nlag > 0){p <- max(lags)
 }else{
  p=0; inlcude.mean=TRUE}
ist <- p+1
nT <- length(x)
nobe <- nT-p
X <- NULL
if(include.mean)X <- rep(1,nobe)
if(nlag > 0){
for (i in 1:nlag){
 ii = lags[i]
 X <- cbind(X,x[(ist-ii):(nT-ii)])
 }
}
X <- as.matrix(X)

if(p > 0){c1 <- paste("lag",lags,sep="")
 }else{c1 <- NULL}
if(include.mean) c1 <- c("cnt",c1)
colnames(X) <- c1
k <- ncol(X)
y <- x[ist:nT]
m1 <- lm(y~-1+X)
coef <- m1$coefficients
#### Estimation
build <- function(parm){
 if(k == 1){
  dlm(FF=matrix(rep(1,k),nrow=1),V=exp(parm[1]), W=exp(parm[2]),
  GG=diag(k),m0=coef, C0 = diag(k),JFF=matrix(c(1:k),nrow=1),X=X)
  }else{
  dlm(FF=matrix(rep(1,k),nrow=1),V=exp(parm[1]), W=diag(exp(parm[-1])),
  GG=diag(k),m0=coef, C0 = diag(k),JFF=matrix(c(1:k),nrow=1),X=X)
  }
 }
mm <- dlmMLE(y,rep(-0.5,k+1),build)
par <- mm$par; value <- mm$value; counts <- mm$counts

cat("par:", par,"\n")
tvAR <-list(par=par,value=value,counts=counts,convergence=mm$convergence,message=mm$message)
}

#' Filtering and Smoothing for Time-Varying AR Models
#' 
#' This function performs forward filtering and backward smoothing for a fitted time-varying AR model with parameters in 'par'.
#' @param x a time series of data.
#' @param lags the lag of AR order.
#' @param par the fitted time-varying AR models. It can be an object returned by function. \code{tvAR}.
#' @param include.mean a logical value indicating whether the constant terms are included.
#' @examples
#' x1=rnorm(100)
#' x2=arima.sim(n = 63, list(ar = c(0.8897, -0.4858)),sd = sqrt(0.1796))
#' x=c(x1,x2)
#' est=tvAR(x,2)
#' tvARFiSm(x,2,TRUE,est$par)
#' @return \code{trARFiSm} function return values returned by function \code{dlmFilter} and \code{dlmSmooth}.
#' @export
"tvARFiSm" <- function(x,lags=c(1),include.mean=TRUE,par){ 
if(is.matrix(x))x <- c(x[,1])
nlag <-  length(lags)
if(nlag > 0){p <- max(lags)
 }else{
  p=0; inlcude.mean=TRUE}
ist <- p+1
nT <- length(x)
nobe <- nT-p
X <- NULL
if(include.mean)X <- rep(1,nobe)
if(nlag > 0){
for (i in 1:nlag){
 ii = lags[i]
 X <- cbind(X,x[(ist-ii):(nT-ii)])
 }
}
X <- as.matrix(X)
k <- ncol(X)
y <- x[ist:nT]
m1 <- lm(y~-1+X)
coef <- m1$coefficients
### Model specification
if(k == 1){
  tvAR <- dlm(FF=matrix(rep(1,k),nrow=1),V=exp(par[1]), W=exp(par[2]),
  GG=diag(k),m0=coef, C0 = diag(k),JFF=matrix(c(1:k),nrow=1),X=X)
  }else{
  tvAR <- dlm(FF=matrix(rep(1,k),nrow=1),V=exp(par[1]), W=diag(exp(par[-1])),
  GG=diag(k),m0=coef, C0 = diag(k),JFF=matrix(c(1:k),nrow=1),X=X)
  }
mF <- dlmFilter(y,tvAR)
mS <- dlmSmooth(mF)

tvARFiSm <- list(filter=mF,smooth=mS)
}

#' Estimating of Random-Coefficient AR Models
#'
#' Estimate random-coefficient AR models.
#' @param x a time series of data.
#' @param lags the lag of AR models. This is more flexible than using order. It can skip unnecessary lags.
#' @param include.mean a logical value indicating whether the constant terms are included.
#' @return \code{rcAR} function returns a list with following components:
#' \item{par}{estimated parameters.}
#' \item{se.est}{standard errors.}
#' \item{residuals}{residuals.}
#' \item{sresiduals}{standardized residuals.}
#' @examples
#' x1=rnorm(100)
#' x2=arima.sim(n = 63, list(ar = c(0.8897, -0.4858)),sd = sqrt(0.1796))
#' x=c(x1,x2)
#' est=rcAR(x,2,TRUE)
"rcAR" <- function(x,lags=c(1),include.mean=TRUE){
if(is.matrix(x))x <- c(x[,1])
if(include.mean){
 mu <- mean(x)
 x <- x-mu
 cat("Sample mean: ",mu,"\n")
}
nlag <-  length(lags)
if(nlag  < 1){lags <- c(1); nlag <- 1}
p <- max(lags)
ist <- p+1
nT <- length(x)
nobe <- nT-p
X <- NULL
for (i in 1:nlag){
 ii = lags[i]
 X <- cbind(X,x[(ist-ii):(nT-ii)])
 }
X <- as.matrix(X)
k <- ncol(X)
y <- x[ist:nT]
m1 <- lm(y~-1+X)
par <- m1$coefficients
par <- c(par,rep(-3.0,k),-0.3)
##cat("initial estimates: ",par,"\n")
##
rcARlike <- function(par,y=y,X=X){
k <- ncol(X)
nobe <- nrow(X)
sigma2 <- exp(par[length(par)])
if(k > 1){
 beta <- matrix(par[1:k],k,1)
 yhat <- X%*%beta
 gamma <- matrix(exp(par[(k+1):(2*k)]),k,1)
 sig <- X%*%gamma
 }else{
  beta <- par[1]; gamma <- exp(par[2])
  yhat <- X*beta
  sig <- X*gamma
  }
 at <- y-yhat
 v1 <- sig+sigma2
 ll <- sum(dnorm(at,mean=rep(0,nobe),sd=sqrt(v1),log=TRUE))
rcARlike <- -ll
}
#
mm <- optim(par,rcARlike,y=y,X=X,hessian=TRUE)
est <- mm$par
H <- mm$hessian
Hi <- solve(H)
se.est <- sqrt(diag(Hi))
tratio <- est/se.est
tmp <- cbind(est,se.est,tratio)
cat("Estimates:","\n")
print(round(tmp,4))
### Compute residuals and standardized residuals
sigma2 <- exp(est[length(est)])
if(k > 1){
 beta <- matrix(est[1:k],k,1)
 yhat <- X%*%beta
 gamma <- matrix(exp(est[(k+1):(2*k)]),k,1)
 sig <- X%*%gamma
 }else{
  beta <- est[1]; gamma <- exp(est[2])
  yhat <- X*beta
  sig <- X*gamma
  }
 at <- y-yhat
 v1 <- sig+sigma2
 sre <- at/sqrt(v1)
#
rcAR <- list(par=est,se.est=se.est,residuals=at,sresiduals=sre)
}


