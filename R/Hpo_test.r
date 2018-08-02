#' Threshold Nonlinearity Test 
#'
#' Threshold nonlinearity test.
#' @references
#' Tsay, R. (1989) Testing and Modeling Threshold Autoregressive Processes. \emph{Journal of the American Statistical Associations} \strong{84}(405), 231-240.
#'
#' @param y a time seris.
#' @param p AR order.
#' @param d delay for the threhosld variable.
#' @param thrV threshold variable.
#' @param ini inital number of data to start RLS estimation.
#' @param include.mean a logical value for including constant terms.
#' @return \code{thr.test} returns a list with components: 
#' \item{F-ratio}{F statistic.}
#' \item{df}{the numerator and denominator degrees of freedom.}
#' \item{ini}{initial number of data to start RLS estimation.} 
#' @examples
#' arorder=rep(1,2)
#' ar.coef=matrix(c(0.7,-0.8),2,1)
#' y=uTAR.sim(100,arorder,ar.coef,1,0)
#' thr.test(y$series,1,1,y$series,40,TRUE)
#' @export
"thr.test" <- function(y,p=1,d=1,thrV=NULL,ini=40,include.mean=T){
if(is.matrix(y))y <- y[,1]
nT <- length(y)
if(p < 1)p <-1
if(d < 1) d <-1
ist <- max(p,d)+1
nobe <- nT-ist+1
if(length(thrV) < nobe){
 cat("SETAR model is entertained","\n")
 thrV <- y[(ist-d):(nT-d)]
 }else{thrV <- thrV[1:nobe]}
#
Y <- y[ist:nT]
X <- NULL
for (i in 1:p){
  X <- cbind(X,y[(ist-i):(nT-i)])
}
if(include.mean){X <- cbind(rep(1,nobe),X)}
X <- as.matrix(X)
kx <- ncol(X)
m1 <- sort(thrV,index.return=TRUE)
idx <- m1$ix
Y <- Y[idx]
if(kx == 1){X <- X[idx]
       }else{ X <- X[idx,]}
while(kx > ini){ini <- ini+10}
if(kx == 1){
   XpX = sum(X[1:ini]^2)
   XpY <- sum(X[1:ini]*Y[1:ini])
   betak <- XpY/XpX
   Pk <- 1/XpX
   resi <- Y[1:ini] - X[1:ini]*betak
   }
   else{ XpX <- t(X[1:ini,])%*%X[1:ini,]
         XpY <- t(X[1:ini,])%*%as.matrix(Y[1:ini])
	 Pk <- solve(XpX)
	 betak <- Pk%*%XpY
	 resi <- Y[1:ini]-X[1:ini,]%*%betak
	}
### RLS to obtain standardized residuals
sresi <- resi
if(ini < nobe){
 for (i in (ini+1):nobe){
   if(kx == 1){xk <- X[i]
               er <- xk*betak-Y[i]
	       deno <- 1+xk^2/Pk
	       betak <- betak -Pk*er/deno
	       tmp <- Pk*xk
	       Pk <- tmp*tmp/deno
	       sresi <- c(sresi,er/sqrt(deno))
	       }
	       else{
	       xk <- matrix(X[i,],kx,1)
	       tmp <- Pk%*%xk
	       er <- t(xk)%*%betak - Y[i]
	       deno <- 1  + t(tmp)%*%xk
	       betak <- betak - tmp*c(er)/c(deno)
	       Pk <- Pk -tmp%*%t(tmp)/c(deno)
	       sresi <- c(sresi,c(er)/sqrt(c(deno)))
	       }
	}
###cat("final esimates: ",betak,"\n")
### testing
Ys <- sresi[(ini+1):nobe]
if(kx == 1){x <- data.frame(X[(ini+1):nobe])
           }else{ x <- data.frame(X[(ini+1):nobe,])}
m2 <- lm(Ys~.,data=x)
Num <- sum(Ys^2)-sum(m2$residuals^2)
Deno <- sum(m2$residuals^2)
h <- max(1,p+1-d)
df1 <- p+1
df2 <- nT-d-ini-p-h
Fratio <- (Num/df1)/(Deno/df2)

cat("Threshold nonlinearity test for (p,d): ",c(p,d),"\n")
pv <- pf(Fratio,df1,df2,lower.tail=FALSE)
cat("F-ratio and p-value: ",c(Fratio,pv),"\n")
}
else{cat("Insufficient sample size to perform the test!","\n")}

thr.test <- list(F.ratio = Fratio,df=c(df1,df2),ini=ini)
}



#' Tsay Test for Nonlinearity
#'
#' Perform Tsay (1986) nonlinearity test.
#' @references
#' Tsay, R. (1986) Nonlinearity tests for time series. \emph{Biometrika} \strong{73}(2), 461-466.
#' @param y time series.
#' @param p AR order.
#' @return The function outputs the F statistic, p value, and the degrees of freedom. The null hypothsis is there is no nonlinearity.
#' @examples
#' y=MSM.sim(100,c(1,1),0.7,-0.5,c(0.5,0.6),c(1,1),c(0,0),500)
#' Tsay(y$series,1)
#' @export
"Tsay" <- function(y,p=1){
if(is.matrix(y))y=y[,1]
nT <- length(y)
if(p < 1)p <- 1
ist <- p+1
Y <- y[ist:nT]
ym <- scale(y,center=TRUE,scale=FALSE)
X <- rep(1,nT-p)
for (i in 1:p){
 X <- cbind(X,ym[(ist-i):(nT-i)])
 }
colnames(X) <- c("cnst",paste("lag",1:p))
m1 <- lm(Y~-1+.,data=data.frame(X))
resi <- m1$residuals
XX <- NULL
for (i in 1:p){
 for (j in 1:i){
  XX <- cbind(XX,ym[(ist-i):(nT-i)]*ym[(ist-j):(nT-j)])
  }
 }
XR <- NULL
for (i in 1:ncol(XX)){
 mm <- lm(XX[,i]~.,data=data.frame(X))
 XR <- cbind(XR,mm$residuals)
 }
colnames(XR) <- paste("crossp",1:ncol(XX))
m2 <- lm(resi~-1+.,data=data.frame(XR))
resi2 <- m2$residuals
SS0 <- sum(resi^2)
SS1 <- sum(resi2^2)
df1 <- ncol(XX)
df2 <- nT-p-df1-1
Num <- (SS0-SS1)/df1
Deno <- SS1/df2
F <- Num/Deno
pv <- 1-pf(F,df1,df2)
cat("Non-linearity test & its p-value: ",c(F,pv),"\n")
cat("Degrees of freedom of the test: ",c(df1,df2),"\n")
}


#' Backtest
#'
#' @param m1 a time-series model object.
#' @param rt the time series.
#' @param orig forecast origin.
#' @param h forecast horizon.
#' @param xre the independent variables.
#' @param fixed parameter constriant.
#' @param inc.mean a logicial value for constant term of the model. Default is TRUE.
#' @return The function returns a list with following components:
#' \item{orig}{the starting forecast origin.}
#' \item{err}{observed value minus fitted value.}
#' \item{rmse}{RMSE of out-of-sample forecasts.}
#' \item{mabso}{mean absolute error of out-of-sample forecasts.}
#' \item{bias}{bias of out-of-sample foecasts.}
#' @export
"backtest" <- function(m1,rt,orig,h,xre=NULL,fixed=NULL,include.mean=TRUE){
# m1: is a time-series model object
# orig: is the starting forecast origin
# rt: the time series
# xre: the independent variables
# h: forecast horizon
# fixed: parameter constriant
# inc.mean: flag for constant term of the model.
#
regor=c(m1$arma[1],m1$arma[6],m1$arma[2])
seaor=list(order=c(m1$arma[3],m1$arma[7],m1$arma[4]),period=m1$arma[5])
nT=length(rt)
if(orig > nT)orig=nT
if(h < 1) h=1
rmse=rep(0,h)
mabso=rep(0,h)
bias <- rep(0,h)
nori=nT-orig
err=matrix(0,nori,h)
jlast=nT-1
for (n in orig:jlast){
 jcnt=n-orig+1
 x=rt[1:n]
 if (is.null(xre))
  pretor=NULL else pretor=xre[1:n,]
 mm=arima(x,order=regor,seasonal=seaor,xreg=pretor,fixed=fixed,include.mean=include.mean)
 if (is.null(xre)){nx=NULL} 
   else {nx=matrix(xre[(n+1):(n+h),],h,ncol(xre))}
 fore=predict(mm,h,newxreg=nx)
 kk=min(nT,(n+h))
# nof is the effective number of forecats at the forecast origin n.
 nof=kk-n
 pred=fore$pred[1:nof]
 obsd=rt[(n+1):kk]
 err[jcnt,1:nof]=obsd-pred
}
#
for (i in 1:h){
iend=nori-i+1
tmp=err[1:iend,i]
mabso[i]=sum(abs(tmp))/iend
rmse[i]=sqrt(sum(tmp^2)/iend)
bias[i]=mean(tmp)
}
print("RMSE of out-of-sample forecasts")
print(rmse)
print("Mean absolute error of out-of-sample forecasts")
print(mabso)
print("Bias of out-of-sample foecasts")
print(bias)
backtest <- list(origin=orig,error=err,rmse=rmse,mabso=mabso,bias=bias)
}



#' Backtest for Univariate TAR Models
#'
#' Perform back-test of a univariate SETAR model.
#' @param model SETAR model.
#' @param orig forecast origin.
#' @param h forecast horizon.
#' @param iter number of iterations.
#' @examples
#' arorder=rep(1,2)
#' ar.coef=matrix(c(0.7,-0.8),k,1)
#' y=uTAR.sim(100,arorder,ar.coef,1,0)
#' est=uTAR.est(y$series,arorder,0,1)
#' backTAR(est,50,1,3000)
#' @return \code{backTAR} returns a list of components:
#' \item{model}{SETAR model.}
#' \item{error}{prediction errors.}
#' \item{State}{predicted states.}
#' @export
"backTAR" <- function(model,orig,h=1,iter=3000){
y <- model$data
order <- model$arorder
thr1 <- c(model$thr)
thr2 <- c(min(y)-1,thr1,max(y)+1)
cnst <- model$cnst
d1 <- model$delay
nT <- length(y)
if(orig < 1)orig <- nT-1
if(orig >= nT) orig <- nT-1
Err <- NULL
Y <- c(y,rep(mean(y),h))
nregime <- length(thr1)+1
State <- rep(nregime,(nT-orig))
for (n in orig:(nT-1)){
 y1 <- y[1:n]
 v1 <- y[n-d1]
 jj = nregime
 for (j in 1:nregime){
  if((v1 > thr2[j]) && (v1 <= thr2[j+1]))jj=j
  }
 State[n-orig+1] = jj
 m1 <- uTAR.est(y1,arorder=order,thr=thr1,d=d1,include.mean=cnst,output=FALSE)
 m2 <- uTAR.pred(m1,n,h=h,iterations=iter,output=FALSE)
 er <- Y[(n+1):(n+h)]-m2$pred
 if(h == 1) {Err=c(Err,er)
            }
	    else{Err <- rbind(Err,matrix(er,1,h))}
 }
## summary statistics
MSE <- NULL
MAE <- NULL
Bias <- NULL
if(h == 1){
 MSE = mean(Err^2)
 MAE = mean(abs(Err))
 Bias = mean(Err)
 nf <- length(Err)
 }else{
  nf <- nrow(Err)
  for (j in 1:h){
   err <- Err[1:(nf-j+1),j]
   MSE <- c(MSE,mean(err^2))
   MAE <- c(MAE,mean(abs(err)))
   Bias <- c(Bias,mean(err))
   }
  }
cat("Starting forecast origin: ",orig,"\n")
cat("1-step to ",h,"-step out-sample forecasts","\n")
cat("RMSE: ",sqrt(MSE),"\n")
cat(" MAE: ",MAE,"\n")
cat("Bias: ",Bias,"\n")
cat("Performance based on the regime of forecast origins: ","\n")
for (j in 1:nregime){
 idx=c(1:nf)[State==j]
 if(h == 1){
  Errj=Err[idx]
  MSEj = mean(Errj^2)
  MAEj = mean(abs(Errj))
  Biasj = mean(Errj)
  }else{
   MSEj <- NULL
   MAEj <- NULL
   Biasj <- NULL
   for (i in 1:h){
    idx=c(1:(nf-i+1))[State[1:(nf-i+1)]==j]
    err = Err[idx,i]
    MSEj <- c(MSEj,mean(err^2))
    MAEj <- c(MAEj,mean(abs(err)))
    Biasj <- c(Biasj,mean(err))
    }
   }
  cat("Summary Statistics when forecast origins are in State: ",j,"\n")
  cat("Number of forecasts used: ",length(idx),"\n")
  cat("RMSEj: ",sqrt(MSEj),"\n")
  cat(" MAEj: ",MAEj,"\n")
  cat("Biasj: ",Biasj,"\n")
  }
backTAR <- list(model=model,error=Err,State=State)
}


#' Rank-Based Portmanteau Tests
#'
#' Performs rank-based portmanteau statistics.
#' @param zt time series.
#' @param lag the maximum lag to calculate the test statistic.
#' @param ouput a logical value for output. Default is TRUE.
#' @return \code{rankQ} function outputs the test statistics and p-values for Portmanteau tests, and returns a list with components:
#' \item{Qstat}{test statistics.}
#' \item{pv}{p-values.}
#' @examples
#' y=rnorm(1000)
#' rankQ(y,10,output=TRUE)
#' @export
"rankQ" <- function(zt,lag=10,output=TRUE){
nT <- length(zt)
Rt <- rank(zt)
erho <- vrho <- rho <- NULL
rbar <- (nT+1)/2
sse <- nT*(nT^2-1)/12
Qstat <- NULL
pv <- NULL
rmRt <- Rt - rbar
deno <- 5*(nT-1)^2*nT^2*(nT+1)
n1 <- 5*nT^4
Q <- 0
for (i in 1:lag){
 tmp <- crossprod(rmRt[1:(nT-i)],rmRt[(i+1):nT])
 tmp <- tmp/sse
 rho <- c(rho,tmp)
 er <- -(nT-i)/(nT*(nT-1))
 erho <- c(erho,er)
 vr <- (n1-(5*i+9)*nT^3+9*(i-2)*nT^2+2*i*(5*i+8)*nT+16*i^2)/deno
 vrho <- c(vrho,vr)
 Q <- Q+(tmp-er)^2/vr
 Qstat <- c(Qstat,Q)
 pv <- c(pv,1-pchisq(Q,i))
 }
if(output){
 cat("Results of rank-based Q(m) statistics: ","\n")
 Out <- cbind(c(1:lag),rho,Qstat,pv)
 colnames(Out) <- c("Lag","ACF","Qstat","p-value")
 print(round(Out,3))
 }

rankQ <- list(rho = rho, erho = erho, vrho=vrho,Qstat=Qstat,pv=pv)
}


#' F Test for Nonlinearity
#'
#' Compute the F-test statistic for nonlinearity
#' @param x time series.
#' @param order AR order.
#' @param thres threshold value.
#' @return The function outputs the test statistic and its p-value, and return a list with components:
#' \item{test.stat}{test statistic.}
#' \item{p.value}{p-value.}
#' \item{order}{AR order.}
#' @examples
#' y=rnorm(100)
#' F.test(y,2,0)
#' @export
"F.test" <- function(x,order,thres=0.0){
if (missing(order)) order=ar(x)$order
 m <- order
 ist <- m+1
 nT <- length(x)
 y <- x[ist:nT]
 X <- NULL
 for (i in 1:m) X <- cbind(X,x[(ist-i):(nT-i)])
 lm1 <- lm(y~X)
 a1 <- summary(lm1)
 coef <- a1$coefficient[-1,]
 idx <- c(1:m)[abs(coef[,3]) > thres]
 jj <- length(idx)
 if(jj==0){
        idx <- c(1:m)
	jj <- m
        }
 for (j in 1:jj){
   for (k in 1:j){
    X <- cbind(X,x[(ist-idx[j]):(nT-idx[j])]*x[(ist-idx[k]):(nT-idx[k])])
    }
   }
 lm2 <- lm(y~X)
 a2 <- anova(lm1,lm2)
 list(test.stat = signif(a2[[5]][2],4),p.value=signif(a2[[6]][2],4),order=order)
}

#' ND Test
#'
#' Compute the ND test statistic of Pena and Rodriguez (2006, JSPI).
#' @param x time series.
#' @param m the maximum number of lag of correlation to test.
#' @param p AR order.
#' @param q MA order.
#' @references
#' Pena, D., and Rodriguez, J. (2006) A powerful Portmanteau test of lack of fit for time series. series. \emph{Journal of American Statistical Association}, 97, 601-610.
#' @return \code{PRnd} function outputs the ND test statistic and its p-value.
#' @examples
#' y=arima.sim(n=500,list(ar=c(0.8,-0.6,0.7)))
#' PRnd(y,10,3,0)
#' @export
"PRnd" <- function(x,m=10,p=0,q=0){
pq <- p+q
nu <- 3*(m+1)*(m-2*pq)^2
de <- (2*m*(2*m+1)-12*(m+1)*pq)*2
alpha <- nu/de
#cat("alpha: ",alpha,"\n")
beta <- 3*(m+1)*(m-2*pq)
beta <- beta/(2*m*(m+m+1)-12*(m+1)*pq)
#
n1 <- 2*(m/2-pq)*(m*m/(4*(m+1))-pq)
d1 <- 3*(m*(2*m+1)/(6*(m+1))-pq)^2
r1 <- 1-(n1/d1)
lambda <- 1/r1
cat("lambda: ",lambda,"\n")
if(is.matrix(x))x <- c(x[,1])
nT <- length(x)
adjacf <- c(1)
xwm <- scale(x,center=T,scale=F)
deno <- crossprod(xwm,xwm)
if(nT > m){
 for (i in 1:m){
   nn <- crossprod(xwm[1:(nT-i)],xwm[(i+1):nT])
   adjacf <- c(adjacf,(nT+2)*nn/((nT-i)*deno))
   }
## cat("adjacf: ",adjacf,"\n")
 }
 else{
  adjacf <- c(adjacf,rep(0,m))
   }
Rm <- adjacf
tmp <- adjacf
for (i in 1:m){
 tmp <- c(adjacf[i+1],tmp[-(m+1)])
 Rm <- rbind(Rm,tmp)
 }
#cat("Rm: ","\n")
#print(Rm)

tst <- -(nT/(m+1))*log(det(Rm))
a1 <- alpha/beta
nd <- a1^(-r1)*(lambda/sqrt(alpha))*(tst^r1-a1^r1*(1-(lambda-1)/(2*alpha*lambda^2)))
pv <- 2*(1-pnorm(abs(nd)))
cat("ND-stat & p-value ",c(nd,pv),"\n")
}


